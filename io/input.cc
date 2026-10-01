/**
 * @file input.cc
 * @brief Implementation of input.h (the functions are documented in the header).
 */
#include "input.h"
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
	*total_time = atof( argv[8] );
	*delta_t = atof( argv[9] );
	*number_steps_samplig = atoi( argv[10] );
	
}


//----------------------------------------------------------------------
//input for time evolution of Ising chain in longitudinal (hx) e transversal magnetic field (hz)	

void 
get_data( char* argv[] , int *state , int *N , double *J , double *hx , double *hz, double *ttotal , double *tstep , int *nmeas , int *bonddim, int *localvscluster)
	{
	double hxArray[] { 0. , 0.1 , 0.2 , 0.4 , -0.1 , -0.2 , -0.4};
	double hzArray[] { 0.25 , 0.5 , 0.75 , 1. , 1.25 , 1.5 , 1.75 , 2. , 3.};
	
	int hxChoice = atoi( argv[4] );
	int hzChoice = atoi( argv[5] );	
		
	*state = atoi( argv[1] );	
	*N = atoi( argv[2] );
	*J = atof( argv[3] );
	*hx = hxArray[ hxChoice ];
	*hz = hzArray[ hzChoice ];
	*ttotal =  atof( argv[6] );
	*tstep = atof( argv[7] );
	*nmeas =  atoi( argv[8] );
	*bonddim =  atoi( argv[9] );
	*localvscluster = atoi( argv[10] );
	}


//----------------------------------------------------------------------
//input measuring full counting statistics of a certain oberservable 

void 
get_data_meas( char* argv[] , int *N , int *hxChoice , int *hzChoice , double *tstep , int *nmeas , int *numberPoints , int *maxLength , int *localvscluster )
	{
	*N = atoi( argv[1] );												//number of spin
	*hxChoice = atoi( argv[2] );
	*hzChoice = atoi( argv[3] );
	*tstep = atof( argv[4] );											//time step
	*nmeas =  atoi( argv[5] );											//after how many steps is measured the state
	*numberPoints = atoi( argv[6] );									//number of points to evaluate the genering function in the range [-pi,+pi]
	*maxLength = atoi( argv[7] );										//length of subsystem up to which measure the generating function	
	*localvscluster = atoi( argv[8] );
	}


//----------------------------------------------------------------------
//input for measuring entanglement entropy of a one dimensional system

void 
get_data_entropy( char* argv[] , int *N , double *tstep , int *nmeas , int *localvscluster)
	{
	*N = atoi( argv[1] );
	*tstep = atof( argv[2] );
	*nmeas =  atoi( argv[3] );
	*localvscluster = atoi( argv[4] );
	}
