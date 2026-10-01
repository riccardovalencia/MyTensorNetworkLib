#ifndef MYTN_BOSONS_PERFORM_DMRG_H
#define MYTN_BOSONS_PERFORM_DMRG_H

// DMRG drivers for the bosonic quantum east model (bQEM).

#include <itensor/all.h>
#include <iostream>
#include <fstream>	//output file
#include <sstream> // for ostringstream
#include <string>
#include <iomanip>

using namespace std;
using namespace itensor;

// DMRG for the ground state of H, starting from *ground_state (overwritten with the result).
// The bond dimension is increased by scaling_bond_dimension until the energy changes by less than
// precision_dmrg or max_bond_dimension is reached. Energies/variances are written to
// "energy_size..._n0<n0>.dat" and "delta_energy_size...dat"; intermediate states to "ground_state_file_n0<n0>_chi<chi>".
// physical_args : "size", "cut_off_fock_space", "n0" (Int), "s", "c" (Real) - used for file names.
// numerical_args: "bond_dimension", "scaling_bond_dimension", "max_bond_dimension",
//                 "number_sweep_fixed_bond_dimension" (Int), "precision_dmrg",
//                 "lower_bound_singular_values", "original_noise" (Real).
// Returns the final bond dimension.
int
perform_DMRG(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

// Same as perform_DMRG, but allows up to 3 sweeps per bond dimension (instead of 2) before increasing it.
int
perform_DMRG_soft(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

#endif
