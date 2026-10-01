/**
 * @file time_evolution.cc
 * @brief Implementation of time_evolution.h (the functions are documented in the header).
 */
#include "time_evolution.h"
#include "../dynamics/purified_state.h"
#include "../mps/gates.h"
#include "../mps/mps_tools.h"
#include <itensor/all.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;
using namespace itensor;


MPS
TEBD_lindblad_time_evolve(MPS psi_t, vector<BondGate> gates , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");
    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

    	gateTEvol( gates , dt , dt , psi_t , TEBD_args); 

        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();


                int j = jket[0];

                cerr <<  j << " "; 

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }

            gateTEvol( gates , dt , dt , psi_t , TEBD_args); 
        }


        if(normalize)
        {
            double norm = compute_norm_purified_impurity(&psi_t);
            psi_t /= norm;
        }


        if ( (k+1) % steps_save_state == 0)
        {
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


// -----------------------------------------------------------------
// Time evolution via TEBD of a a density matrix unfolded as an MPS.
// The structure of the unfolded density matrix is 
// | | | | |
// o-o-o-o-o-   (ket)
// |
// o-o-o-o-o-   (bra)
// | | | | |
// The algorithm handles long range interactions + local dissipative channel acting on the first physical site (corresponding
// to the N and N+1 site in the unfolded MPS)

MPS
TEBD_long_range_int_lindblad_time_evolve(MPS psi_t, vector<BondGate> gates_H, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim     = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");
	int N2 = length(psi_t);
	int N = int(N2/2);


    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;
		cerr << "Time : " << t << "\n";
		// long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
        for (BondGate g : gates_H)
        {   
            int j1 = g.i1();
			int j2 = g.i2();
			int j ;

			// cerr << "Applying swap gate between " << j1 << " " << j2 << "\n";

			// the swap gates move around the impurity site, to make it near 
			// at then end, it will be back to the original position
			// it acts on bra
			if(j1 <= N && j2 <= N)
			{
				j = min(j1,j2);
				ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
				auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
				
				psi_t.set(j,U);
				psi_t.set(j+1,S*V);
				swap_gate(&psi_t,j,j+1,cut_off,MaxDim);

			}

			// it acts on ket

			else
			{
				j = max(j1,j2);
				ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
				auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
				psi_t.set(j-1,U);
				psi_t.set(j,S*V);
				swap_gate(&psi_t,j-1,j,cut_off,MaxDim);
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

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			// long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
			for (BondGate g : gates_H)
			{   

				int j1 = g.i1();
				int j2 = g.i2();
				int j ;
				// cerr << "Applying swap gate between " << j1 << " " << j2 << "\n";

				// the swap gates move around the impurity site, to make it near 
				// at then end, it will be back to the original position
				// it acts on bra
				if(j1 <= N && j2 <= N)
				{
					j = min(j1,j2);
					ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
					auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
					
					psi_t.set(j,U);
					psi_t.set(j+1,S*V);
					swap_gate(&psi_t,j,j+1,cut_off,MaxDim);

				}

				// it acts on ket

				else
				{
					j = max(j1,j2);
					ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
					auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
					psi_t.set(j-1,U);
					psi_t.set(j,S*V);
					swap_gate(&psi_t,j-1,j,cut_off,MaxDim);
				}
			}
            
        }


        if(normalize)
        {
            double norm = compute_norm_purified_impurity(&psi_t);
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


// -----------------------------------------------------------------
// Time evolution via application of MPO of a a density matrix unfolded as an MPS.
// The structure of the unfolded density matrix is 
// | | | | |
// o-o-o-o-o-   (ket)
// |
// o-o-o-o-o-   (bra)
// | | | | |
// The algorithm handles long range interactions + local dissipative channel acting on the first physical site (corresponding
// to the N and N+1 site in the unfolded MPS).
// The coherent part is applied via a first-order approximation of exp(-iH t) = 1 -i H t.

MPS
MPO_lindblad_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;
	
	MPO Ht =  -1_i * dt_step * H; 


    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        MPS dpsi = applyMPO(Ht,psi_t,{"Method=","DensityMatrix","MaxDim=",MaxDim,"Cutoff=",cut_off});
        psi_t = sum(psi_t, dpsi.noPrime());
	
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


			dpsi = applyMPO(Ht,psi_t,{"Method=","DensityMatrix","MaxDim=",MaxDim,"Cutoff=",cut_off});
        	psi_t = sum(psi_t, dpsi.noPrime());
        }


        if(normalize)
        {
            double norm = compute_norm_purified_impurity(&psi_t);
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
