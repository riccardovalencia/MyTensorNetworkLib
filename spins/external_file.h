#ifndef EXTERNAL_FILE_H
#define EXTERNAL_FILE_H

#include <itensor/all.h>
#include <iostream>
#include <string>
#include <ctime>

#include <itensor/all.h>
#include <iostream>
#include <sstream> // for ostringstream
#include <vector>
#include <string>
#include <iomanip>
#include <complex>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
//files to save the information about sites and the state
void 
build_file_TEBD( stringstream * ,  stringstream * , const int  );

//----------------------------------------------------------------------
//files to save information for measuring generating function and moments of a certain full counting statistics
void
build_file_full_counting( stringstream * ,  stringstream * , stringstream * , stringstream * ,  stringstream * ,  stringstream * , const int );

//----------------------------------------------------------------------
//file to save entanglement entropy of a one dimensional system
void
build_file_entanglement_entropy( stringstream * ,  stringstream * , stringstream * , const int );
	
//----------------------------------------------------------------------
//print input for time evolution in Ising model with longitudinal and transversal magnetic fields
void
print_input(int , double , double , double , double , double , double , double );

//----------------------------------------------------------------------
//print information during time evolution: time reached, time needed to do a single step, time needed in total, max bond dimension, entanglement entropy
double
print_info( time_t , time_t  , MPS * , const int , const int , const int , const double );

//----------------------------------------------------------------------
//print at time fixed the generating function of a certain probability distribution


#endif
