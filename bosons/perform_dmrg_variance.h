#ifndef MYTN_BOSONS_PERFORM_DMRG_VARIANCE_H
#define MYTN_BOSONS_PERFORM_DMRG_VARIANCE_H

#include <itensor/all.h>
#include <iostream>
#include <fstream>	//output file
#include <sstream> // for ostringstream
#include <string>
#include <iomanip>

using namespace std;
using namespace itensor;

// Fixed schedule of 8 DMRG rounds (bond dimension 10 -> 50, decreasing noise) on H, starting from *ground_state
// (overwritten). Meant to target excited states close to energy_target, passing e.g. H = (H_bQEM - energy_target)^2.
// Relative energy changes are appended to "excited_states_delta_energy_size...dat", energies to "energy_size...dat".
// physical_args : "size", "cut_off_fock_space", "n0" (Int), "s", "c" (Real) - used for file names.
// numerical_args: "bond_dimension", "scaling_bond_dimension", "max_bond_dimension",
//                 "number_sweep_fixed_bond_dimension" (Int), "precision_dmrg",
//                 "lower_bound_singular_values", "original_noise" (Real).
// Returns the variance of H on the final state.
double
perform_DMRG_variance(double energy_target, MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

// Initial guess for perform_DMRG_variance: product state with int(energy_target) bosons on random sites.
// Note: sets the site tensors without link indices (*ground_state_variance should not be used for
// operations that require links before DMRG).
void
initialize_excited_state( MPS *ground_state_variance, const SiteSet sites, const int size, const double energy_target, const int cut_off_fock_space);


#endif
