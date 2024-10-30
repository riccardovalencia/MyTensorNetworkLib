#ifndef GET_DATA_H
#define GET_DATA_H
//get_data.h

using namespace std;

//----------------------------------------------------------------------
//input for time evolution of Ising chain in longitudinal (hx) e transversal magnetic field (hz)	

void 
get_data( char ** , int * , int * , double * , double * , double * , double * , double * , int * , int * , int *);

//----------------------------------------------------------------------
//input measuring full counting statistics of a certain oberservable 

void 
get_data_meas( char** , int * , int * , int * , double * , int * , int * , int * , int * );
#endif

//----------------------------------------------------------------------
//input for measuring entanglement entropy of a one dimensional system

void 
get_data_entropy( char** , int * , double * , int * , int *);

