#ifndef MYTN_SPINS_FULL_COUNTING_STATISTICS_H
#define MYTN_SPINS_FULL_COUNTING_STATISTICS_H

// Full counting statistics of the subsystem magnetization S^x_A = sum_{j in A} S^x_j of a spin-1/2 chain,
// with A a block of size+1 sites centered in the chain.
// Generating function G(theta) = <exp(i theta S^x_A)>, evaluated on numberPoints values of theta in [-pi, pi)
// (denser grid close to theta = 0, see theta_step).

#include <itensor/all.h>
#include <sstream>
#include <vector>
#include <complex>

using namespace std;
using namespace itensor;

// Write G[size][k] to the file named *save_file: one row per theta, columns theta, G(size=0), G(size=1), ...
void
printing_generating_function( const stringstream *save_file , const int numberPoints , const int maxLength , vector<vector<double> > &G );

// Step between consecutive values of theta (index col) on the grid of numberPoints points.
double
theta_step( int col , int numberPoints );

// Pure state psi (MPS) of N sites: append Re/Im of G(theta) for a block of size+1 sites to singleGreal/singleGimag.
void
generaring_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPS* psi , const SpinHalf sites );

// Mixed state rho (MPO) of N sites: as above, with G(theta) = Tr(rho exp(i theta S^x_A)).
void
generaring_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPO* rho , const SpinHalf sites );

// Pure state: compute G for block sizes 0..maxLength-1 and write it to "<save_real>n.dat" / "<save_imag>n.dat"
// (n labels the time step). Skipped if the files already exist.
void
measure_generating_function( MPS* psi , const SpinHalf sites , int N , int n , int maxLength , int numberPoints , stringstream* save_real , stringstream* save_imag );

// Mixed state: as above, writing to "Termal_N<N>_hx_<hx>_hz_<hz>_GF_{real,imag}.dat".
void
measure_generating_function( MPO* rho , const SpinHalf sites , int N , int maxLength , int numberPoints , double hx , double hz );

// MPO of sum_{j=start}^{start+size} S^x_j.
MPO
build_totalSx( const SpinHalf sites , const int start , const int size );

// First four cumulants of S^x_A for block sizes 0..maxLength-1, written (real and imaginary parts) to
// "<saveRealMoments>n.dat" / "<saveImagMoments>n.dat", one row per block size. Skipped if the files exist.
void
measuring_moments( MPS *psi , const SpinHalf sites , const int N , const int n , const int maxLength , stringstream* saveRealMoments , stringstream* saveImagMoments );

#endif
