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
#include "../../spin_boson.h"
#include <filesystem>

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// Dynamics of isolated (no dissipation) matter interacting via a single-mode cavity (bosonic)
// It is possible to select either Dicke or Tavis-Cummings coupling via the string variable
// photon_matter_coupling.
// photon_matter_coupling = "dicke"
//      H = omega0 a^dag a + h \sum_j \sigma_j^z + g/sqrt{N} (a+a^\dag)\sum_j \sigma_j^x
// photon_matter_coupling = "tavis"
//      H = omega0 a^dag a + h \sum_j \sigma_j^z + g/sqrt{N} \sum_j (\sigma_j^+ a + h.c.)

// Benchmarked with ED: successful

//------------------------------------------------------------------MAIN--------------------------------------------------------------------

int main(int argc , char* argv[]){
	
    int N = 3;
    int max_occ = 6;
    double omega0 = 1.;
    double h = 0.5;
    double g = 1.5;
    string photon_matter_coupling =  "dicke"; // choose either "tavis" g\sum_j(a s_j^+ + h.c.); "dicke" g(a+a^\dag)\sum_j s_j^x

    // observables to measure
    vector<string> name_obs;
    name_obs.push_back("fidelity");
    name_obs.push_back("Sx");
    name_obs.push_back("Sz");
    name_obs.push_back("Na");
    name_obs.push_back("MaxD");

    // spin coherent state (theta is the polar angle, phi is the azimuthal angle)
    double theta = 0.5 * M_PI;
    double phi   = 0.;

    double T  = 50.;   // total time
    double dt = 0.005;   // time step
    double cut_off = 1E-8; // cut_off TEBD
    int maxDim = 50;       // max bond dimension

    int steps_measure = 10;
    int total_steps = int(T / dt);

    SiteSet sites = custom_spin_boson(N+1,max_occ);
    MPS psi = initialize_spin_boson_state(sites , 0 , theta, phi);
    
    cerr << "Magnetization along x" << endl;
    vector<double> mj = measure_magnetization(&psi,sites,"x");
    for(double m : mj) cerr << " " << m ;
    cerr << "\nMagnetization along z" << endl;
    mj = measure_magnetization(&psi,sites,"z");
    for(double m : mj) cerr << " " << m ;
    cerr << "\n";

    vector<BondGate> gates_matter = gates_photon_matter(sites , omega0 ,h , g/sqrt(N), dt, "short-range",photon_matter_coupling);
    vector<BondGate> gates_pm = gates_photon_matter(sites , omega0 ,h , g/sqrt(N), dt, "long-range",photon_matter_coupling);

    MPS psi_t0 = psi;
    vector<double> overlap_t;

    cerr << setprecision(10);

    fs::create_directories("data");
    string file_obs    = tinyformat::format("data/%s_obs_TN_N%d_maxocc%d_omega%.2f_h%.2f_g%.2f.txt", photon_matter_coupling, N, max_occ, omega0, h, g);
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
        for (BondGate g : gates_matter)
        {
            int j = g.i1();
            // change orthogonality center to minimize error
            psi.position(j);
            psi.normalize();
            // apply gate
            ITensor AA = psi(j)*psi(j+1)*g.gate();
            auto [U,S,V] = svd(noPrime(AA),inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
            
            psi.set(j,U);
            psi.set(j+1,S*V);
        }

        // long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
        for (BondGate g : gates_pm)
        {   
            int j = g.i2();
            // change orthogonality center to minimize error
            psi.position(j);
            psi.normalize();
            // apply gate
            ITensor AA = psi(j-1)*psi(j)*g.gate();
            auto [U,S,V] = svd(noPrime(AA),inds(psi(j-1)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
            
            psi.set(j-1,U);
            psi.set(j,S*V);
            swap_gate(&psi,j-1,j,cut_off,maxDim);
        }


        if (k % steps_measure == 0)
        {
            save_file << t ;
            for(string name : name_obs)
            {
                if(name=="fidelity")
                {
                    double fidelity = abs(innerC(psi_t0,psi));
                    fidelity *= fidelity;
                    save_file << " " << fidelity;
                }
                if(name=="Sx")
                {
                    mj = measure_magnetization(&psi,sites,"x");
                    double Sx = 0.;
                    for(int idx=1; idx < mj.size() ; idx++) Sx += mj[idx];
                    save_file << " " << Sx/N ;
                }

                if(name=="Sz")
                {
                    mj = measure_magnetization(&psi,sites,"z");
                    double Sz = 0.;
                    for(int idx=1; idx < mj.size() ; idx++) Sz += mj[idx];
                    save_file << " " << Sz/N ;
                }

                if(name=="Na")
                {
                    mj = measure_magnetization(&psi,sites,"x");
                    save_file << " " << mj[0]/N;
                }
                if(name=="MaxD")
                {
                    save_file << " " << maxLinkDim(psi) << endl;
                }
            }
            // save_file << endl;
            // double overlap = abs(innerC(psi_t0,psi));
            // save_file << t  overlap.real() << endl; 
            // save_imag_F << overlap.imag() << endl;
            // save_time   << t << endl;
            cerr << t << " " << maxLinkDim(psi) << endl;
        }
        
    }       

    // save_real_F.close();
    // save_imag_F.close();
    // save_time.close();
    save_file.close();
    return 0;

}