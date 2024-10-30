#include "get_data.h"
#include <iostream>

using namespace std;
	
//----------------------------------------------------------------------
// input of the physical information of the bosonic quantum east model chain
void 
get_data_system(  char* argv[] , int *size , double *s , double *c , int *n0 , double *simmetry_sector , int *cut_off_fock_space )
{
	*size = atoi( argv[1] );
	*s = atof( argv[2] );
	*c = atof( argv[3] );
	*n0 = atoi( argv[4] );
	*simmetry_sector = atof( argv[5] );
	*cut_off_fock_space = atoi( argv[6] );
}

//----------------------------------------------------------------------
// input about numerical quantities for DMRG in the bosonic quantum east model chain
void
get_data_DMRG(  char* argv[]  , int *bond_dimension, double *lower_bound_singular_values, double *scaling_bond_dimension, double *precision_dmrg, int *max_bond_dimension)
{
	*bond_dimension = atoi( argv[7] );
	*lower_bound_singular_values = atof( argv[8] );
	*scaling_bond_dimension = atof(argv[9] );
	*precision_dmrg = atof( argv[10] );
	*max_bond_dimension = atoi( argv[11] );
}

//----------------------------------------------------------------------
// input about numerical quantities for TEBD in the bosonic quantum east model chain
void
get_data_TEBD(  char* argv[] , int *max_bond_dimension, double *lower_bound_singular_values, double *total_time, double *delta_t, int *number_steps_samplig)
{
	*max_bond_dimension = atoi( argv[6] );
	*lower_bound_singular_values = atof( argv[7] );
	*total_time = atof( argv[10] );
	*delta_t = atof( argv[9] );
	*number_steps_samplig = atoi( argv[10] );
	
}
