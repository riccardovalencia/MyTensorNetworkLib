#include "get_data.h"
#include <iostream>

using namespace std;
	
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

