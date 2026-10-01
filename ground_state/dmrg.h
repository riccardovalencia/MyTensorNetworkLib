/**
 * @file dmrg.h
 * @brief DMRG drivers: ground states with increasing bond dimension, excited states, Fermi sea.
 *
 * The perform_dmrg* drivers read their parameters from two ITensor Args:
 * - physical_args : "size", "cut_off_fock_space", "n0" (Int), "s", "c" (Real) - used in file names;
 * - numerical_args: "bond_dimension", "scaling_bond_dimension", "max_bond_dimension",
 *                   "number_sweep_fixed_bond_dimension" (Int), "precision_dmrg",
 *                   "lower_bound_singular_values", "original_noise" (Real).
 * They write energies and convergence data to files in the working directory.
 */
#ifndef MYTN_GROUND_STATE_DMRG_H
#define MYTN_GROUND_STATE_DMRG_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Ground state of H with a bond dimension increased by scaling_bond_dimension until the
 *        relative energy change is below precision_dmrg or max_bond_dimension is reached.
 *
 * Output files: "energy_size..._n0<n0>.dat" (energies and variances), "delta_energy_size...dat"
 * (convergence), and the intermediate states "ground_state_file_n0<n0>_chi<chi>".
 * @throws ITError if max_bond_dimension is exceeded.
 *
 * @param ground_state         Initial state, overwritten with the result.
 * @param H                    Hamiltonian.
 * @param sites                Site set.
 * @param set_output_precision Digits written to the output files.
 * @return The final bond dimension.
 */
int
perform_dmrg(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

/** @brief As perform_dmrg, allowing up to 3 sweeps per bond dimension (instead of 2) before increasing it. */
int
perform_dmrg_soft(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

/** @brief As perform_dmrg for the mean-field problem; writes "meanfield_*energy_size...dat". */
void
perform_dmrg_meanfield(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

/**
 * @brief Fixed schedule of 8 DMRG rounds (bond dimension 10 -> 50, decreasing noise) on H, meant to
 *        target excited states close to energy_target by passing e.g. H = (H_0 - energy_target)^2.
 *
 * Relative energy changes are appended to "excited_states_delta_energy_size...dat", energies to
 * "energy_size...dat".
 *
 * @param energy_target Target energy (written to the output file).
 * @param ground_state  Initial state, overwritten with the result.
 * @return Variance of H on the final state.
 */
double
perform_dmrg_variance(double energy_target, MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

/**
 * @brief Initial guess for perform_dmrg_variance: product state with int(energy_target) bosons on
 *        random sites (at most cut_off_fock_space per site).
 * @param ground_state_variance State to overwrite; must have link indices (e.g. randomMPS(sites)).
 */
void
set_excited_state_guess( MPS *ground_state_variance, const SiteSet sites, const int size, const double energy_target, const int cut_off_fock_space);

/**
 * @brief Ground state of H at fixed filling, starting from make_electron_product_state.
 * @param H        Hamiltonian (e.g. make_tight_binding_mpo).
 * @param sites    Electron site set.
 * @param Nupfill  Number of up electrons.
 * @param Ndnfill  Number of down electrons.
 * @param sweeps   DMRG sweeps.
 * @param min_varH Maximum accepted energy variance.
 * @throws ITError if the variance exceeds min_varH or the filling is not reproduced.
 */
MPS
find_fermi_sea(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, const Sweeps sweeps, double min_varH);

#endif
