#ifndef MYTN_BOSONS_GET_DATA_H
#define MYTN_BOSONS_GET_DATA_H

// Command-line parsing for the bQEM drivers.
//get_data.h

using namespace std;

//----------------------------------------------------------------------
// input of the physical information of the bosonic quantum east model chain

void 
get_data_system(  char* argv[] , int *size , double *s , double *c, int *n0 , double *simmetry_sector , int *cut_off_fock_space );

//----------------------------------------------------------------------
// input about numerical quantities for DMRG in the bosonic quantum east model chain
void
get_data_DMRG(  char* argv[]  , int *bond_dimension, double *lower_bound_singular_values, double *scaling_bond_dimension, double *precision_dmrg, int *max_bond_dimension);

//----------------------------------------------------------------------
// input about numerical quantities for TEBD in the bosonic quantum east model chain
void
get_data_TEBD(  char* argv[] , int *max_bond_dimension, double *lower_bound_singular_values, double *total_time, double *delta_t, int *number_steps_samplig);

#endif
