/**
 * @file full_counting_statistics.h
 * @brief Full counting statistics of the subsystem magnetization S^x_A = sum_{j in A} S^x_j of a
 *        spin-1/2 chain, A being a block of size+1 sites centered in the chain.
 *
 * The generating function G(theta) = <exp(i theta S^x_A)> is evaluated on numberPoints values of
 * theta in [-pi, pi), with a denser grid close to theta = 0 (see theta_step).
 * Results are written with printing_generating_function (io/output.h).
 */
#ifndef MYTN_ANALYSIS_FULL_COUNTING_STATISTICS_H
#define MYTN_ANALYSIS_FULL_COUNTING_STATISTICS_H

#include <itensor/all.h>
#include <sstream>
#include <vector>

using namespace std;
using namespace itensor;

/**
 * @brief Step between theta_col and theta_{col+1} on the grid of numberPoints points starting at -pi.
 */
double
theta_step( int col , int numberPoints );

/**
 * @brief Generating function of a pure state for one block size.
 * @param singleGreal, singleGimag Output: Re and Im of G(theta), numberPoints values appended.
 * @param size        Block A has size+1 sites, centered in the chain.
 * @param N           Number of sites.
 * @param numberPoints Number of values of theta.
 * @param psi         Pure state; its orthogonality center is moved.
 * @param sites       Spin-1/2 site set.
 */
void
generating_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPS* psi , const SpinHalf sites );

/** @brief Same for a mixed state rho given as an MPO, G(theta) = Tr(rho exp(i theta S^x_A)). */
void
generating_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPO* rho , const SpinHalf sites );

/**
 * @brief Generating function of a pure state for block sizes 0..maxLength-1, written to
 *        "<save_real>n.dat" and "<save_imag>n.dat". Skipped if the files already exist.
 * @param n Label of the time step, appended to the file names.
 */
void
measure_generating_function( MPS* psi , const SpinHalf sites , int N , int n , int maxLength , int numberPoints , stringstream* save_real , stringstream* save_imag );

/**
 * @brief Generating function of a mixed state for block sizes 0..maxLength-1, written to
 *        "Termal_N<N>_hx_<hx>_hz_<hz>_GF_{real,imag}.dat". Skipped if the files already exist.
 * @param hx, hz Fields, used only in the file names.
 */
void
measure_generating_function( MPO* rho , const SpinHalf sites , int N , int maxLength , int numberPoints , double hx , double hz );

/** @brief MPO of sum_{j=start}^{start+size} S^x_j. */
MPO
make_block_sx_mpo( const SpinHalf sites , const int start , const int size );

/**
 * @brief First four cumulants of S^x_A for block sizes 0..maxLength-1.
 *
 * Writes one row per block size (C1 C2 C3 C4), real and imaginary parts, to
 * "<saveRealMoments>n.dat" and "<saveImagMoments>n.dat". Skipped if the files already exist.
 */
void
measuring_moments( MPS *psi , const SpinHalf sites , const int N , const int n , const int maxLength , stringstream* saveRealMoments , stringstream* saveImagMoments );

#endif
