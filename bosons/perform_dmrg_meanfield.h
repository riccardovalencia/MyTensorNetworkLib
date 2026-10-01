#ifndef MYTN_BOSONS_PERFORM_DMRG_MEANFIELD_H
#define MYTN_BOSONS_PERFORM_DMRG_MEANFIELD_H


#include <itensor/all.h>
#include <iostream>
#include <fstream>	//output file
#include <sstream> // for ostringstream
#include <string>
#include <iomanip>

using namespace std;
using namespace itensor;

// As perform_DMRG for the mean-field (bond dimension 1) problem; writes "meanfield_*energy_size...dat".
// physical_args : "size", "cut_off_fock_space", "n0" (Int), "s", "c" (Real) - used for file names.
// numerical_args: "bond_dimension", "scaling_bond_dimension", "max_bond_dimension",
//                 "number_sweep_fixed_bond_dimension" (Int), "precision_dmrg",
//                 "lower_bound_singular_values", "original_noise" (Real).
void
perform_DMRG_meanfield(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);


#endif
