#include <itensor/all.h>
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <math.h>       /* exp */
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>
#include "mytn.h"
#include <filesystem>

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// Dynamics dissipative (photon losses) matter interacting via a single-mode cavity (bosonic)
// It is possible to select either Dicke or Tavis-Cummings coupling via the string variable
// photon_matter_coupling concerning the unitary part.
// photon_matter_coupling = "dicke"
//      H = omega0 a^dag a + h \sum_j \sigma_j^z + g/sqrt{N} (a+a^\dag)\sum_j \sigma_j^x
// photon_matter_coupling = "tavis"
//      H = omega0 a^dag a + h \sum_j \sigma_j^z + g/sqrt{N} \sum_j (\sigma_j^+ a + h.c.)

// The whole dynamics is obtained by integrating the Lindbland master equation
// \dot{rho} = -i[H,\rho] + k (a \rho a^dag - a^\dag a \rho - \rho a^\dag a),

// Used purification (doubling space) approach


//------------------------------------------------------------------MAIN--------------------------------------------------------------------

int main(int argc , char* argv[]){
	

    if(argc != 2)
    {
        cerr << "Passed wrong number of arguments.\n";
        exit(-1);
    }

    InputGroup input = InputGroup(argv[1],"input");
    int N           = input.getInt("N",10); // length
    int max_occ     = input.getInt("max_occ",25); // number of fermions up spin
    double omega0   = 1.; // photon frequency (energy scale)
    double h        = input.getReal("h",0.2);
    double g_ratio  = input.getReal("g",1.6);
    double kappa    = input.getReal("kappa",0.4);
    double T        = input.getReal("T",100);
    double dt       = input.getReal("dt",0.01);
    double cut_off  = input.getReal("cut_off",1E-14);
    int maxDim      = input.getInt("maxDim",1024);
    string initial_state = input.getString("initial_state","polarized_x");
    string photon_matter_coupling =  "dicke"; // choose either "tavis" g\sum_j(a s_j^+ + h.c.); "dicke" g(a+a^\dag)\sum_j s_j^x

    // mean field critical photon-matter coupling
    double gc = sqrt(0.5 * abs(h)     * (omega0*omega0 + kappa*kappa/4.) / omega0 );
    double g   = g_ratio * gc;

    bool dissipative = (kappa > 1E-10);

    // observables to measure
    vector<string> name_obs;
    name_obs.push_back("norm"); // written for every measurement (column 2)
    name_obs.push_back("Sx");
    name_obs.push_back("Sz");
    name_obs.push_back("Na");
    name_obs.push_back("MaxD");

    // spin coherent state (theta is the polar angle, phi is the azimuthal angle)

    double theta = 0.9 * M_PI;
    double phi   = 0.;


    int steps_measure = 1;
    if (dt < 0.01) steps_measure = int(0.01/dt);
    int total_steps = int(T / dt);

    SiteSet sites_single = custom_spin_boson(N+1,max_occ);
    MPS psi = initialize_spin_boson_state(sites_single , 0 , theta, phi);
    
    cerr << "Magnetization along x" << endl;
    vector<double> mj = measure_magnetization(&psi,sites_single,"x");
    for(double m : mj) cerr << " " << m ;
    cerr << "\nMagnetization along z" << endl;
    mj = measure_magnetization(&psi,sites_single,"z");
    for(double m : mj) cerr << " " << m ;
    

    // doubling space

    SiteSet sites = custom_spin_boson_doubling(N+1,max_occ);
    MPS psi_t = randomMPS(sites);    
    
    // inserting bra (it gets inverted and dag) between [1,N]
    insert_state(&psi_t, psi, 1,   true , true);
    // inserting ket bewteen [N+1,2*N]
    insert_state(&psi_t, psi, N+2, false, false);

    double  norm = compute_norm_purified_impurity(&psi_t);
    psi_t /= norm;
    // operators of interest
    Index sph   = sites(N+2);
    Index sph_p = prime(sites(N+2));
    ITensor Nb  = ITensor(sph,sph_p);
    
    for(int d=1; d <= dim(sph) ; d++) Nb.set( sph(d)  ,sph_p(d),d-1);
    
    sph   = sites(1);
    sph_p = prime(sites(1));

    ITensor Sx = ITensor(sph,sph_p);
    ITensor Sz = ITensor(sph,sph_p);

    Sx.set(sph(1),sph_p(2),1.);
    Sx.set(sph(2),sph_p(1),1.);	
              
    Sz.set(sph(1),sph_p(1),1.);
    Sz.set(sph(2),sph_p(2),-1.);

    vector<complex<double> > number_photons = measure_local_obs_impurity_first_site(&psi_t , Nb, false, 1);
    vector<complex<double> > mj_x  = measure_local_obs_impurity_first_site(&psi_t , Sx, false, 2);
    vector<complex<double> > mj_z  = measure_local_obs_impurity_first_site(&psi_t , Sz, false, 2);
    
    vector<ITensor> Lj;
    int impurity_sites[2] = {N+1,N+2};
    for(int j : impurity_sites)
    {
        sph   = sites(j);
	    sph_p = prime(sites(j));
        ITensor A  = ITensor(sph,sph_p);

        for(int d=1; d < dim(sph) ; d++)
        {
            A.set( sph(d+1)  ,sph_p(d),sqrt(d));
        }
        Lj.push_back( A ); 
    }


    // build gates for ket
    vector<BondGate> gates_matter_single;
    vector<BondGate> gates_pm_single;
    vector<MyBondGateDiss> gates_D;

    if(dissipative)
    {
        cerr << "Dissipative is true.\n" ;
        gates_matter_single = gates_photon_matter(sites_single , omega0 ,h , -g/sqrt(N), dt/2., "short-range",photon_matter_coupling);
        gates_pm_single     = gates_photon_matter(sites_single , omega0 ,h , -g/sqrt(N), dt/2., "long-range",photon_matter_coupling);
        gates_D = gates_dissipative_impurity(sites , Lj, kappa , dt);
    }              
    else
    {
        cerr << "Dissipative is false.\n" ;
        gates_matter_single = gates_photon_matter(sites_single , omega0 ,h , -g/sqrt(N), dt, "short-range",photon_matter_coupling);
        gates_pm_single     = gates_photon_matter(sites_single , omega0 ,h , -g/sqrt(N), dt, "long-range",photon_matter_coupling);
    }

    // using the gates above to define the gates of our interest
    vector<MyBondGate> gates_matter;
    vector<MyBondGate> gates_pm;
    
    // ket dynamics
    for( BondGate g : gates_matter_single)
    {
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
        // position along the doubled space
        int inew_1 = N + 1 + i1;
        int inew_2 = N + 1 + i2;
        
        ITensor gate = g.gate();
        Index si1 = sites_single(i1);
        Index si2 = sites_single(i2);
        Index sinew1 = sites(inew_1);
        Index sinew2 = sites(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));
        
        MyBondGate gnew = MyBondGate(sites,{inew_1,inew_2},0,gate);
        gnew.modify_gate(gate);
        gates_matter.push_back(gnew);
    }

    // bra dynamics
    for( BondGate g : gates_matter_single)
    {
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
        // position along the doubled space
        int inew_1 = N + 2 - i1;
        int inew_2 = N + 2 - i2;
        
        ITensor gate = g.gate();
        Index si1 = sites_single(i1);
        Index si2 = sites_single(i2);
        Index sinew1 = sites(inew_1);
        Index sinew2 = sites(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));
        // gate = swapPrime(gate,0,1); // NOT SURE

        
        MyBondGate gnew = MyBondGate(sites,{inew_1,inew_2},0,gate);
        gnew.modify_gate(dag(gate));
        gates_matter.push_back(gnew);
    }

    // ket dynamics
    for( BondGate g : gates_pm_single)
    {
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
        // position along the doubled space
        int inew_1 = N + 1 + i1;
        int inew_2 = N + 1 + i2;
        
        ITensor gate = g.gate();
        Index si1 = sites_single(i1);
        Index si2 = sites_single(i2);
        Index sinew1 = sites(inew_1);
        Index sinew2 = sites(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));
        
        MyBondGate gnew = MyBondGate(sites,{inew_1,inew_2},0,gate);
        gnew.modify_gate(gate);
        gates_pm.push_back(gnew);
    }

    // bra dynamics
    for( BondGate g : gates_pm_single)
    {
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
        // position along the doubled space
        int inew_1 = N + 2 - i1;
        int inew_2 = N + 2 - i2;
        
        ITensor gate = g.gate();
        Index si1 = sites_single(i1);
        Index si2 = sites_single(i2);
        Index sinew1 = sites(inew_1);
        Index sinew2 = sites(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));
        // gate = swapPrime(gate,0,1); // NOT SURE

        MyBondGate gnew = MyBondGate(sites,{inew_1,inew_2},0,gate);
        gnew.modify_gate(dag(gate));
        gates_pm.push_back(gnew);

    }


    cerr << "Initialized all gates.\n";
  
    MPS psi_t0 = psi_t;
    vector<double> overlap_t;

    cerr << setprecision(10);

    fs::create_directories("data");
    string file_obs    = tinyformat::format("data/%s_obs_TN_N%d_maxocc%d_omega%.2f_h%.2f_gratio%.2f_kappa%.2f_maxDim%d.txt", photon_matter_coupling,N, max_occ, omega0, h, g_ratio, kappa, maxDim);
    ofstream save_file( file_obs) ;

    save_file << "# t";
    for(string name : name_obs) save_file << " . " + name;
    save_file << endl;

    save_file << setprecision(4);

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = (k+1)*dt;

        // short-range interaction part Hamiltonian
        // no need of swap gates
        // cerr << "Applying single-site terms.\n";
        for (MyBondGate g : gates_matter)
        {
            vector<int> jn = g.jn();
            int j = *min_element(jn.begin(), jn.end());
            // change orthogonality center to minimize error
            psi_t.position(j);
            // apply gate
            ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
            auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
            
            psi_t.set(j,U);
            psi_t.set(j+1,S*V);
        }

        // long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
        // cerr << "Applying long-range photon matter interactions.\n";

        for (MyBondGate g : gates_pm)
        {   
            vector<int> jn = g.jn();
            int j1 = jn[0];
            int j2 = jn[1];

            // acting on ket part
            if(j1 < j2)
            {
                // cerr << "Case 1.\n";
                int j = j2;
                // change orthogonality center to minimize error
                psi_t.position(j);
                // cerr << "Acting on sites " << j-1 << " " << j << "\n";
                // cerr << g.gate() << "\n";
                // apply gate
                ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
                auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                
                psi_t.set(j-1,U);
                psi_t.set(j,S*V);
                
                swap_gate(&psi_t,j-1,j,cut_off,maxDim);

                // cerr << "Swapping " << j-1 << " - " << j << "\n";
            }

            // acting on bra
            else
            {
                // cerr << "Case 2.\n";
                int j = j2;
                // cerr << j << "\n";
                // change orthogonality center to minimize error
                psi_t.position(j);
                // cerr << "Acting on sites " << j << " " << j+1 << "\n";
                // cerr << g.gate() << "\n";

                // apply gate
                ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
                auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);


                swap_gate(&psi_t,j,j+1,cut_off,maxDim);
                // cerr << "Swapping " << j << " - " << j+1 << "\n";
            }

                        
        }


        if(dissipative)
        {

            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                int j = jket[0];
                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);
                AA = AA + dpsi;
                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);
            }

            // short-range interaction part Hamiltonian
            // no need of swap gates
            // cerr << "Applying single-site terms.\n";
            for (MyBondGate g : gates_matter)
            {
                vector<int> jn = g.jn();
                int j = *min_element(jn.begin(), jn.end());
                psi_t.position(j);
                ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
                auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);
            }

            // no need of swap gates
            for (MyBondGate g : gates_pm)
            {   
                vector<int> jn = g.jn();
                int j1 = jn[0];
                int j2 = jn[1];

                // acting on ket part
                if(j1 < j2)
                {
                    int j = j2;
                    // change orthogonality center to minimize error
                    psi_t.position(j);
                    // apply gate
                    ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
                    auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                    psi_t.set(j-1,U);
                    psi_t.set(j,S*V);
                    swap_gate(&psi_t,j-1,j,cut_off,maxDim);
                }

                // acting on bra
                else
                {
                    int j = j2;
                    psi_t.position(j);
                    ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
                    auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
                    psi_t.set(j,U);
                    psi_t.set(j+1,S*V);
                    swap_gate(&psi_t,j,j+1,cut_off,maxDim);
                }

                            
            }


        }

        if (k % steps_measure == 0)
        {
            norm = compute_norm_purified_impurity(&psi_t);
            cerr << t << " " << norm << "\n";
            save_file << t << " " << norm;
            for(string name : name_obs)
            {
                
                if(name=="Sx")
                {
                    mj_x  = measure_local_obs_impurity_first_site(&psi_t , Sx, false, 2);
                    save_file << " " << real(mj_x[0])/norm; 
                }

                if(name=="Sz")
                {
                    mj_z  = measure_local_obs_impurity_first_site(&psi_t , Sz, false, 2);
                    save_file << " " << real(mj_z[0])/norm; 
                }

                if(name=="Na")
                {
                    number_photons = measure_local_obs_impurity_first_site(&psi_t , Nb, false, 1);
                    save_file << " " << real(number_photons[0])/norm; 
                }
                if(name=="MaxD")
                {
                    save_file << " " << maxLinkDim(psi_t) << endl;
                }
            }
            
        }
        
    }       

    save_file.close();
    return 0;

}