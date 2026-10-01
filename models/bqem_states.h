/**
 * @file bqem_states.h
 * @brief "Super-bosonic" states of the bosonic quantum east model (bQEM), dressed operators,
 *        and overlaps between states computed in different truncations or sectors.
 *
 * See models/bqem.h for the conventions (s, c, n0, symmetry).
 */
#ifndef MYTN_MODELS_BQEM_STATES_H
#define MYTN_MODELS_BQEM_STATES_H

#include <itensor/all.h>
#include <complex>
#include <tuple>

using namespace std;
using namespace itensor;

/**
 * @brief State |0>^k (x) |n0> (x) |ground_state> (x) |0>...: a ground state of L sites inserted
 *        after the site with n0 bosons, in a chain of N sites.
 * @param ground_state State of L sites (site set ground_sites).
 * @param sites        Site set of the N-site chain.
 * @param ground_sites Site set of ground_state.
 * @param L            Number of sites of ground_state.
 * @param N            Number of sites of the chain.
 * @param n0           Occupation of site k+1.
 * @param k            Number of empty sites before the excited site.
 */
MPS
super_bosonic_state(MPS ground_state, const SiteSet sites, const SiteSet ground_sites, const int L, const int N, const int n0, const int k);

/**
 * @brief Super-bosonic coherent state: coherent state |alpha> on site 1, dressed by switching
 *        on the bQEM hopping adiabatically (adiabatic_transformation_linear_protocol).
 * @param psi_coherent   State to overwrite (also returned).
 * @param sites_coherent Boson site set.
 * @param alpha          Coherent amplitude.
 * @param s, c           Target bQEM parameters.
 * @param dt             Time step of the ramp.
 * @param T              Duration of the ramp.
 */
MPS
super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const complex<double> alpha, const double s, const double c, double dt, double T);

/** @brief As super_bosonic_coherent_state, starting from a squeezed vacuum with parameter alpha on site 1. */
MPS
super_bosonic_squeezed_state( MPS *psi_coherent, const SiteSet sites_coherent, const double alpha, const double s, const double c, double dt, double T);

/**
 * @brief Dress an operator with a linear ramp of the bQEM hopping from 0 to e^{-s}:
 *        O -> U^dag O U with U the time-ordered evolution of duration T.
 * @param O     Operator (MPO) to dress.
 * @param sites Boson site set.
 * @param dt    Time step of the ramp.
 * @param T     Duration of the ramp.
 */
MPO
bqem_dressed_operator(MPO O, const SiteSet sites, const double s, const double c, double dt, double T);

/**
 * @brief Overlap between ground states in sectors with different n0, and energy variance of
 *        the embedded state.
 * @param psi1, psi2     States (sites1, sites2) with n0 = n0_1 and n0_2.
 * @param size           Number of sites.
 * @param lambda         Fock-space cutoff.
 * @param symmetry_sector, s, c bQEM parameters (see compute_variance_hamiltonian_bqem).
 * @param symmetry_sector_dir   Folder of the symmetry-sector data files (see compute_variance_hamiltonian_bqem).
 * @return {|<psi1|psi2>|, variance of H on the state with smaller n0 embedded in the larger space}.
 */
tuple<double, double>
scalar_product_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir);

/**
 * @brief Overlap between states computed with different Fock-space cutoffs; the state with the
 *        smaller cutoff is embedded in the larger space.
 * @param symmetry_sector_dir Folder of the symmetry-sector data files (see compute_variance_hamiltonian_bqem).
 * @return {|<psi1|psi2>|, variance of H on the embedded state}.
 */
tuple<double, double>
scalar_product_different_cutoff( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir);

#endif
