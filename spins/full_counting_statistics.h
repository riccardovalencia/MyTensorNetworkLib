#ifndef FULL_COUNTING_STATISTICS_H
#define FULL_COUNTING_STATISTICS_H

#include <itensor/all.h>
#include <vector>
#include <string>
#include <complex>

#include <itensor/all.h>
#include <iostream>
#include <sstream> // for ostringstream
#include <vector>
#include <string>
#include <iomanip>
#include <complex>

using namespace std;
using namespace itensor;

void 
printing_generating_function( const stringstream * ,  const int , const int , vector<vector<double> > *);

//----------------------------------------------------------------------
//return theta where to evaluate generating function

double 
theta_Step( int , int );

//----------------------------------------------------------------------
//measure of the generating function of a certain probability distribution function at time fixed and size of subsystem
void
generaring_function_sim_size( vector<double> * , vector<double> * , int , int , int , MPS* , const SpinHalf );

//----------------------------------------------------------------------
//measure of generating function and saving of the information	
void	
measure_generating_function( MPS* , const SpinHalf , int , int , int , int , stringstream* , stringstream* );

//----------------------------------------------------------------------
//return the MPO of the total magnetization of a spin-1/2 system of a subsystem of size "size"
MPO
build_totalSx( const SpinHalf , const int ,  const int );

//----------------------------------------------------------------------
//measure the first 4 moments of the full counting statistics of the total magnetization of a spin-1/2 system
void
measuring_moments( MPS * , const SpinHalf , const int , const int , const int , stringstream* , stringstream* );

#endif
