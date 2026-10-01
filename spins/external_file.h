#ifndef MYTN_SPINS_EXTERNAL_FILE_H
#define MYTN_SPINS_EXTERNAL_FILE_H

// File names and console/file output for the Ising-chain simulations of the spins module.

#include <itensor/all.h>
#include <sstream>
#include <ctime>

using namespace std;
using namespace itensor;

// sites_file = "sites_N<N>", psi_file = "psi_N<N>_nstep".
void
build_file_TEBD( stringstream *sites_file , stringstream *psi_file , const int N );

// As build_file_TEBD, plus "N<N>_GF_real", "N<N>_GF_imag", "N<N>_Moments_real", "N<N>_Moments_imag".
void
build_file_full_counting( stringstream *sites_file , stringstream *psi_file , stringstream *save_real , stringstream *save_imag , stringstream *saveRealMoments , stringstream *saveImagMoments , const int N );

// As build_file_TEBD, plus save_file = "N<N>_entropy.dat".
void
build_file_entanglement_entropy( stringstream *sites_file , stringstream *psi_file , stringstream *save_file , const int N );

// Write the simulation parameters to "input.txt" (ITensor InputGroup format).
void
print_input( int N , double J , double hx , double hz , double ttotal , double tstep , double nmeas , double bonddim );

// Print time reached, timings, max bond dimension and half-chain entanglement entropy (natural log),
// and return the entropy.
double
print_info( time_t time_elapsed_step , time_t time_elapsed_total , MPS *psi , const int N , const int nmeas , const int n , const double tstep );

#endif
