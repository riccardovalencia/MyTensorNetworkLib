/**
 * @file bosonic_east_model_states.h
 * @brief "Super-bosonic" states of the bosonic quantum east model (bosonic east model), dressed operators,
 *        and overlaps between states computed in different truncations or sectors.
 *
 * See models/bosonic_east_model.h for the conventions (J = e^{-s}, U = 1 - 2c, n0, symmetry).
 */
#ifndef MYTN_MODELS_BOSONIC_EAST_MODEL_STATES_H
#define MYTN_MODELS_BOSONIC_EAST_MODEL_STATES_H

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
make_super_bosonic_state(MPS ground_state, const SiteSet sites, const SiteSet ground_sites, const int L, const int N, const int n0, const int k);

/**
 * @brief Super-bosonic coherent state: coherent state |alpha> on site 1, dressed by switching
 *        on the bosonic east model hopping adiabatically (evolve_adiabatic_linear_ramp).
 * @param psi_coherent   State to overwrite (also returned).
 * @param sites_coherent Boson site set.
 * @param alpha          Coherent amplitude.
 * @param J_target       Final facilitated hopping amplitude (e^{-s}).
 * @param U              Density-density coefficient (1 - 2c).
 * @param dt             Time step of the ramp.
 * @param T              Duration of the ramp.
 */
MPS
make_super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const complex<double> alpha, const double J_target, const double U, double dt, double T);

/**
 * @brief As make_super_bosonic_coherent_state, starting from a squeezed vacuum on site 1.
 * @param psi_coherent   State to overwrite (also returned).
 * @param sites_coherent Boson site set.
 * @param alpha          Squeezing parameter r of the initial state.
 * @param J_target       Final facilitated hopping amplitude (e^{-s}).
 * @param U              Density-density coefficient (1 - 2c).
 * @param dt             Time step of the ramp.
 * @param T              Duration of the ramp.
 */
MPS
make_super_bosonic_squeezed_state( MPS *psi_coherent, const SiteSet sites_coherent, const double alpha, const double J_target, const double U, double dt, double T);

/**
 * @brief Dress an operator with a linear ramp of the bosonic east model hopping from 0 to J_target:
 *        O -> U O U^dag with U the time-ordered evolution of duration T (the dressed O acts on the
 *        dressed states U|psi> as O acts on |psi>). Built with heisenberg_step (dynamics/heisenberg.h).
 * @param O        Operator (MPO) to dress.
 * @param sites    Boson site set.
 * @param J_target Final facilitated hopping amplitude (e^{-s}).
 * @param U        Density-density coefficient (1 - 2c).
 * @param dt       Time step of the ramp.
 * @param T        Duration of the ramp.
 */
MPO
make_dressed_operator(MPO O, const SiteSet sites, const double J_target, const double U, double dt, double T);

/**
 * @brief Overlap between ground states in sectors with different n0: the state with the larger n0 is
 *        projected onto the site set of the other one (make_resized_state).
 * @param psi1, psi2     States, with n0 = n0_1 and n0_2.
 * @param sites1, sites2 Their site sets.
 * @param size           Number of sites.
 * @param n0_1, n0_2     Occupations of the virtual site 0 of the two states.
 * @param lambda         Fock-space cutoff (selects the symmetry eigenvalue in the data file).
 * @param symmetry_sector, s, c bosonic east model parameters (see compute_bosonic_east_model_energy_variance).
 * @param symmetry_sector_dir   Folder of the symmetry-sector data files.
 * @return {|<psi_min|psi_projected>|^2 (0 below 1e-10), variance of H (smaller n0) on the projected state}.
 */
tuple<double, double>
compute_overlap_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir);

/**
 * @brief Overlap between states computed with different Fock-space cutoffs; the state with the
 *        smaller cutoff is embedded in the larger space (make_resized_state).
 * @param psi1, psi2     States.
 * @param sites1, sites2 Their site sets.
 * @param size           Number of sites.
 * @param cut_off_fock_space1, cut_off_fock_space2 Fock-space cutoffs of the two states.
 * @param n0, symmetry_sector, s, c bosonic east model parameters (see compute_bosonic_east_model_energy_variance).
 * @param symmetry_sector_dir Folder of the symmetry-sector data files.
 * @return {|<psi_max|psi_embedded>|^2, variance of H on the embedded state}.
 */
tuple<double, double>
compute_overlap_different_cutoffs( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c, const string symmetry_sector_dir);

#endif
