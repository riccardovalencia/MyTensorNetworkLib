/**
 * @file bosonic_east_model_states.cc
 * @brief Implementation of bosonic_east_model_states.h (the functions are documented in the header).
 */
#include "bosonic_east_model_states.h"
#include "../dof/boson.h"
#include "../dynamics/adiabatic.h"
#include "../models/bosonic_east_model.h"
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


// Initial state |0>^k \otimes |n_0> \otimes |GS(n_0)_L> \otimes |0>^(N-L-k-1)
// It's a state with k 0's, n_0, the ground state of size L in this symmetry sector, and the other 0's.
// The total size of the system is L

MPS
make_super_bosonic_state(MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int L, const int N, const int n0, const int k)
{
	MPS psi_t0 = randomMPS(sites);
	set_fock_excitation( &psi_t0, sites, N , n0, k+1);
	insert_state(&psi_t0, state_to_insert, sites, sites_state_to_insert, k+2, L, N);
	return psi_t0;
}


// return the super-bosonic coherent state

MPS
make_super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const complex<double> alpha, const double s, const double c, double dt, double T)
{
	int L = length(*psi_coherent);
	set_coherent_state_on_site( psi_coherent, sites_coherent, L , 1, alpha );
	evolve_adiabatic_linear_ramp( psi_coherent, sites_coherent, s, c, dt, T);

	return *psi_coherent;
}


// return the super-bosonic coherent state


// Return a super-bosonic squeezed state
MPS
make_super_bosonic_squeezed_state( MPS *psi_squeezed, const SiteSet sites_squeezed, const double alpha, const double s, const double c, double dt, double T)
{
	int L = length(*psi_squeezed);
	set_squeezed_state_on_site( psi_squeezed, sites_squeezed, L , 1, alpha );
	evolve_adiabatic_linear_ramp( psi_squeezed, sites_squeezed, s, c, dt, T);

	return *psi_squeezed;

}


MPO
make_dressed_operator(MPO O,  const SiteSet sites, double s, double c, double dt, double T)
{

	int size = length(sites);
	double J_target = exp(-s);
	int number_step = int(T/dt);

	for(int step=1; step<=int(number_step); step++)
	{
		double J = J_target * step * dt / T;
		cerr << J << endl;
		// make_bosonic_east_model_evolution_mpo( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt)

		MPO expH  = make_bosonic_east_model_evolution_mpo( sites, size , 1, 0.0 , J,  c, dt);
		MPO expHd = make_bosonic_east_model_evolution_mpo( sites, size , 1, 0.0 , J,  c, -1*dt);
		
		O = nmultMPO( O , prime(expH) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);
		O = nmultMPO( expHd , prime(O) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);

	}	
	return O;
}


tuple<double, double>
compute_overlap_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir)
{
    // if you want to compute the variance over the Hamiltonian with smaller n_0 over the state with greater n_0
    int max_n0 = max(n0_1, n0_2);
    int min_n0 = min(n0_1, n0_2);


    vector<ITensor> expanded_state; 

    MPS psi_max_n0; 
    MPS psi_min_n0;
    Boson sites_max_n0;
    Boson sites_min_n0;

    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
    Index physical_previous ;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

    if (min_n0 == n0_1)
    {
        psi_min_n0 = *(psi1);
        psi_max_n0 = *(psi2);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_min_n0 = Boson(ind1);
        sites_max_n0 = Boson(ind2);
    }
    else
    {
        psi_min_n0 = *(psi2);
        psi_max_n0 = *(psi1);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_min_n0 = Boson(ind2);
        sites_max_n0 = Boson(ind1);
    }


    // site 1

    rightindexj = rightLinkIndex( psi_max_n0, 1);
    physical = sites_min_n0(1);
    physical_previous = sites_max_n0(1);
    psi_tocopy_j = psi_max_n0(1);
    
    Tj = ITensor(physical, rightindexj);

    for(int r=1; r<= dim(rightindexj); r++)
    {
        for(int d=1; d<= dim(physical); d++)
        {
            if( d <= dim(physical_previous) ) Tj.set(physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_previous=d, rightindexj=r) );
            else Tj.set(physical=d,rightindexj=r , 0);
        }
    }
    
    expanded_state.push_back(Tj);

    // in the middle of the chain

    for (int j =2 ; j<= size-1 ; j ++)
    {  
        leftindexj  = leftLinkIndex( psi_max_n0, j);
        rightindexj = rightLinkIndex( psi_max_n0, j);
        physical = sites_min_n0(j);
        physical_previous = sites_max_n0(j);
        psi_tocopy_j = psi_max_n0(j);
    
        Tj = ITensor(rightindexj, physical, leftindexj);


        for(int l=1; l<= dim(leftindexj) ; l++)
        {
            for(int r=1; r<= dim(rightindexj); r++)
            {
                for(int d=1; d<= dim(physical); d++)
                {
                    if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d, rightindexj=r) );
                    else Tj.set(leftindexj=l,physical=d, rightindexj=r, 0);
                
                }
            }
        }
 
        expanded_state.push_back(Tj);

    }

    // last site

    leftindexj  = leftLinkIndex( psi_max_n0, size);
    physical = sites_min_n0(size);
    physical_previous = sites_max_n0(size);
    psi_tocopy_j = psi_max_n0(size);
    
    Tj = ITensor(leftindexj, physical);
    
    for(int l=1; l<= dim(leftindexj) ; l++)
    {
            for(int d=1; d<= dim(physical); d++)
            {
                if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d));
                else Tj.set(leftindexj=l,physical=d , 0);
            
            }
    }

    expanded_state.push_back(Tj);


    // define a new MPS of the same kind of sites_max_n0 whose elements are set to be equal to the expanded_state
    // in this way I can use built-in functions in ITensor

    MPS psi_constrained = MPS(sites_min_n0);

    for(int j=1 ; j<=size ; j++) psi_constrained.set(j, expanded_state[j-1]);

    complex<double> overlap_amplitude = innerC(psi_min_n0,psi_constrained);
    double overlap_absolute = overlap_amplitude.real() * overlap_amplitude.real() + overlap_amplitude.imag() * overlap_amplitude.imag();


    cerr << "Computing variance of n0: " << max_n0 << " over Hamiltonian with n0=" << min_n0 << endl;
    double variance_constraned_space = compute_bosonic_east_model_energy_variance(&psi_constrained, sites_min_n0 , size , lambda, min_n0, symmetry_sector, s, c, symmetry_sector_dir);

    if( overlap_absolute < 1E-10 ) overlap_absolute = 0.;

    return {overlap_absolute,variance_constraned_space};
}


tuple<double, double>
compute_overlap_different_cutoffs( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir)
{

    int max_cut_off = max(cut_off_fock_space1, cut_off_fock_space2);

    vector<ITensor> expanded_state; 

    MPS psi_maxcutoff; 
    MPS psi_mincutoff;
    Boson sites_maxcutoff;
    Boson sites_mincutoff;

    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
    Index physical_previous ;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

    if (max_cut_off == cut_off_fock_space1)
    {
        psi_maxcutoff = *(psi1);
        psi_mincutoff = *(psi2);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_maxcutoff = Boson(ind1);
        sites_mincutoff = Boson(ind2);
    }
    else
    {
        psi_maxcutoff = *(psi2);
        psi_mincutoff = *(psi1);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_maxcutoff = Boson(ind2);
        sites_mincutoff = Boson(ind1);
    }

    // site 1

    rightindexj = rightLinkIndex( psi_mincutoff, 1);
    physical = sites_maxcutoff(1);
    physical_previous = sites_mincutoff(1);
    psi_tocopy_j = psi_mincutoff(1);
    
    Tj = ITensor(physical, rightindexj);

    for(int r=1; r<= dim(rightindexj); r++)
    {
        for(int d=1; d<= dim(physical); d++)
        {
            if( d <= dim(physical_previous) ) Tj.set(physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_previous=d, rightindexj=r) );
            else Tj.set(physical=d,rightindexj=r , 0);
        }
    }
    
    expanded_state.push_back(Tj);

    // in the middle of the chain

    for (int j =2 ; j<= size-1 ; j ++)
    {  
        leftindexj  = leftLinkIndex( psi_mincutoff, j);
        rightindexj = rightLinkIndex( psi_mincutoff, j);
        physical = sites_maxcutoff(j);
        physical_previous = sites_mincutoff(j);
        psi_tocopy_j = psi_mincutoff(j);
    
        Tj = ITensor(rightindexj, physical, leftindexj);


        for(int l=1; l<= dim(leftindexj) ; l++)
        {
            for(int r=1; r<= dim(rightindexj); r++)
            {
                for(int d=1; d<= dim(physical); d++)
                {
                    if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d, rightindexj=r) );
                    else Tj.set(leftindexj=l,physical=d, rightindexj=r, 0);
                
                }
            }
        }
 
        expanded_state.push_back(Tj);

    }

    // last site

    leftindexj  = leftLinkIndex( psi_mincutoff, size);
    physical = sites_maxcutoff(size);
    physical_previous = sites_mincutoff(size);
    psi_tocopy_j = psi_mincutoff(size);
    
    Tj = ITensor(leftindexj, physical);
    
    for(int l=1; l<= dim(leftindexj) ; l++)
    {
            for(int d=1; d<= dim(physical); d++)
            {
                if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d));
                else Tj.set(leftindexj=l,physical=d , 0);
            
            }
    }

    expanded_state.push_back(Tj);


    // define a new MPS of the same kind of sites_maxcutoff whose elements are set to be equal to the expanded_state
    // in this way I can use built-in functions in ITensor

    MPS psi_expanded = MPS(sites_maxcutoff);

    for(int j=1 ; j<=size ; j++) psi_expanded.set(j, expanded_state[j-1]);

    complex<double> overlap_amplitude = innerC(psi_maxcutoff,psi_expanded);
    double overlap_absolute = overlap_amplitude.real() * overlap_amplitude.real() + overlap_amplitude.imag() * overlap_amplitude.imag();


    double variance_expanded_space = compute_bosonic_east_model_energy_variance(&psi_expanded, sites_maxcutoff , size , max_cut_off, n0, symmetry_sector, s, c, symmetry_sector_dir);


    return {overlap_absolute,variance_expanded_space};
}
