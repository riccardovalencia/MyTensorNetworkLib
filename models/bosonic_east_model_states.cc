/**
 * @file bosonic_east_model_states.cc
 * @brief Implementation of bosonic_east_model_states.h (interfaces documented in the header, logic commented here).
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


// |0>^k (x) |n0> (x) |ground state of L sites> (x) |0>^(N-L-k-1): Fock excitation n0 on site k+1,
// then the ground state copied onto the sites k+2..k+L+1 (insert_state)

MPS
make_super_bosonic_state(MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int L, const int N, const int n0, const int k)
{
	MPS psi_t0 = randomMPS(sites);
	set_fock_excitation( &psi_t0, sites, n0, k+1);
	insert_state(&psi_t0, state_to_insert, sites, sites_state_to_insert, k+2, L, N);
	return psi_t0;
}


// coherent state on site 1, then the linear adiabatic ramp of the hopping

MPS
make_super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const complex<double> alpha, const double J_target, const double U, double dt, double T)
{
	set_coherent_state_on_site( psi_coherent, sites_coherent, 1, alpha );
	evolve_adiabatic_linear_ramp( psi_coherent, sites_coherent, J_target, U, dt, T);

	return *psi_coherent;
}


// squeezed vacuum on site 1, then the linear adiabatic ramp of the hopping
MPS
make_super_bosonic_squeezed_state( MPS *psi_squeezed, const SiteSet sites_squeezed, const double alpha, const double J_target, const double U, double dt, double T)
{
	set_squeezed_state_on_site( psi_squeezed, sites_squeezed, 1, alpha );
	evolve_adiabatic_linear_ramp( psi_squeezed, sites_squeezed, J_target, U, dt, T);

	return *psi_squeezed;

}


// Heisenberg picture along the ramp J(t) = J_target t / T: at every step O -> e^{i dt H} O e^{-i dt H}
// with first-order MPOs of the evolution (make_bosonic_east_model_evolution_mpo), compressed
// (cutoff 1E-14, at most 1000 states).
MPO
make_dressed_operator(MPO O,  const SiteSet sites, double J_target, double U, double dt, double T)
{
	int size = length(sites);
	int number_step = int(T/dt);

	for(int step=1; step<=int(number_step); step++)
	{
		double J = J_target * step * dt / T;
		MPO expH  = make_bosonic_east_model_evolution_mpo( sites, size , 0.0 , J,  U, dt);
		MPO expHd = make_bosonic_east_model_evolution_mpo( sites, size , 0.0 , J,  U, -1*dt);
		
		O = nmultMPO( O , prime(expH) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);
		O = nmultMPO( expHd , prime(O) ,{"MaxDim",1000,"Cutoff",1E-14}); 
		O.mapPrime(2,1);

	}	
	return O;
}


// The state with the larger n0 is moved onto the site set (local dimensions) of the state with the smaller n0;
// returns its squared overlap with that state (0 below 1E-10) and its energy variance with the smaller n0.

tuple<double, double>
compute_overlap_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir)
{
    bool first_is_min = (n0_1 <= n0_2);
    int min_n0          = min(n0_1, n0_2);
    MPS psi_min_n0      = first_is_min ? *psi1 : *psi2;
    MPS psi_max_n0      = first_is_min ? *psi2 : *psi1;
    Boson sites_min_n0  = Boson(inds(first_is_min ? sites1 : sites2));

    MPS psi_constrained = make_resized_state(psi_max_n0, sites_min_n0);

    double overlap = std::norm(innerC(psi_min_n0, psi_constrained));   // |<psi_min|psi_constrained>|^2
    if(overlap < 1E-10) overlap = 0.;

    double variance = compute_bosonic_east_model_energy_variance(&psi_constrained, sites_min_n0, size, lambda, min_n0, symmetry_sector, s, c, symmetry_sector_dir);
    return {overlap, variance};
}


// The state with the smaller cutoff is moved onto the site set of the state with the larger cutoff;
// returns the squared overlap of the two and the energy variance of the moved state.

tuple<double, double>
compute_overlap_different_cutoffs( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir)
{
    bool first_is_max     = (cut_off_fock_space1 >= cut_off_fock_space2);
    int max_cut_off       = max(cut_off_fock_space1, cut_off_fock_space2);
    MPS psi_maxcutoff     = first_is_max ? *psi1 : *psi2;
    MPS psi_mincutoff     = first_is_max ? *psi2 : *psi1;
    Boson sites_maxcutoff = Boson(inds(first_is_max ? sites1 : sites2));

    MPS psi_expanded = make_resized_state(psi_mincutoff, sites_maxcutoff);

    double overlap = std::norm(innerC(psi_maxcutoff, psi_expanded));

    double variance = compute_bosonic_east_model_energy_variance(&psi_expanded, sites_maxcutoff, size, max_cut_off, n0, symmetry_sector, s, c, symmetry_sector_dir);
    return {overlap, variance};
}
