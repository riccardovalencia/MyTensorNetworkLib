#include "TEBD_long_range.h"
#include "observables.h"
#include "MyClasses.h"
#include "/home/ricval/Documenti/MyLibrary/TDVP/tdvp.h"
#include "/home/ricval/Documenti/MyLibrary/TDVP/basisextension.h"
#include "state_manipulation.h"
#include <itensor/all.h>
#include <iostream>
#include <math.h>       
#include <complex>
using namespace std;
using namespace itensor;


MPS
TDVP_time_evolve(MPS psi_t, SiteSet sites, MPO H , Args TDVP_args, double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int N = length(sites);
	
 
    
    string file_obs    = tinyformat::format("%s.txt",file_root);
    ofstream save_file( file_obs) ;
    save_file << "# t . norm \n";

    file_obs    = tinyformat::format("%s_nj_up.txt",file_root);
    ofstream save_file_nj_up( file_obs) ;

    file_obs    = tinyformat::format("%s_nj_dn.txt",file_root);
    ofstream save_file_nj_dn( file_obs) ;

    int total_steps = int(T/dt);
	int MaxDim = TDVP_args.getInt("MaxDim");
	double cut_off = TDVP_args.getReal("Cutoff");
	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 20;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = (k+1)*dt;

        if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",normalize,
                                      "Silent",true});
            }
        
        // TDVP sweep

        double energy = tdvp(psi_t,H,-1_i*dt,sweeps,{"Truncate",true,
                                        "DoNormalize",normalize,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});


        double normalization = norm(psi_t);
        
        save_file << t << " " << normalization << "\n";
        save_file_nj_up << t ;
		save_file_nj_dn << t ;
		for(int j : range1(N))
		{
			psi_t.position(j);
			ITensor ket = psi_t(j);
			ITensor bra = dag(prime(ket,"Site"));
			ITensor Njup = op(sites,"Nup",j);
			ITensor Njdn = op(sites,"Ndn",j);

			save_file_nj_up << " " << real(eltC(ket * Njup * bra))/normalization ;
        	save_file_nj_dn << " " << real(eltC(ket * Njdn * bra))/normalization ; 
		}
		
		save_file_nj_up << "\n";
		save_file_nj_dn << "\n";
        
        save_file_nj_dn.flush();
        save_file_nj_up.flush();


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
                exit(0);
			}
            // writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       


	return psi_t;
}



MPS
TDVP_split_lindbland_time_evolve(MPS psi_t, MPO Hbra , MPO Hket, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;

	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 10;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        // if(k < 3)
        //     {
        //     // Global subspace expansion
        //     std::vector<Real> epsilonK = {1E-12, 1E-12};
        //     addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
        //                               "Method","DensityMatrix",
        //                               "KrylovOrd",3,
        //                               "DoNormalize",false,
        //                               "Silent",true});
        //     }
        
        // TDVP sweep

        double energy = tdvp(psi_t,Hbra,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

	
		energy = tdvp(psi_t,Hket,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

    
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

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			// if(k < 3)
            // {
            // // Global subspace expansion
            // std::vector<Real> epsilonK = {1E-12, 1E-12};
            // addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
            //                           "Method","DensityMatrix",
            //                           "KrylovOrd",3,
            //                           "DoNormalize",false,
            //                           "Silent",true});
            // }
            // TDVP sweep
			energy = tdvp(psi_t,Hbra,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

	
			energy = tdvp(psi_t,Hket,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}


MPS
TDVP_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;

	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 10;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",false,
                                      "Silent",true});
            }
        
        // TDVP sweep

        double energy = tdvp(psi_t,H,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

    
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

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",false,
                                      "Silent",true});
            }
            // TDVP sweep
			energy = tdvp(psi_t,H,-1_i*dt_step,sweeps,{"Truncate",true,
										"DoNormalize",false,
										"Silent",true,
										"NumCenter",2,
										"ErrGoal",1E-10});
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}
