/**
 * @file full_counting_statistics.h
 * @brief Full counting statistics of the block magnetization S^x_A = sum_{j in A} S^x_j of a
 *        spin-1/2 chain, A being a block of l consecutive sites centered in the chain.
 *
 * The generating function G(theta) = <exp(i theta S^x_A)> is evaluated on a grid of theta in
 * [-pi, pi) that is denser close to theta = 0 (make_theta_grid).
 * Use write_generating_function (io/output.h) to save it.
 */
#ifndef MYTN_ANALYSIS_FULL_COUNTING_STATISTICS_H
#define MYTN_ANALYSIS_FULL_COUNTING_STATISTICS_H

#include <itensor/all.h>
#include <complex>
#include <vector>

using namespace std;
using namespace itensor;

/**
 * @brief Grid of number_points values of theta starting at -pi: a step (pi-1)/(number_points/4) in
 *        the outer quarters and 4/number_points in the central half.
 * @param number_points Number of values (a multiple of 4 gives a symmetric grid).
 */
vector<double>
make_theta_grid( int number_points );

/**
 * @brief First site (1-indexed) of the block of block_size sites centered in a chain of N sites
 *        (shifted to the left when the block cannot be centered exactly).
 * @param N          Number of sites of the chain.
 * @param block_size Number of sites of the block.
 */
int
block_start( int N, int block_size );

/**
 * @brief Generating function G(theta) = <psi| exp(i theta S^x_A) |psi> of a pure state.
 * @param psi        State; its orthogonality center is moved to the start of the block.
 * @param sites      Spin-1/2 site set.
 * @param block_size Number of sites of the block A (centered, see block_start).
 * @param theta      Values of theta (e.g. from make_theta_grid).
 * @return G(theta) for each value of theta.
 */
vector<complex<double> >
compute_generating_function( MPS* psi, const SpinHalf sites, int block_size, const vector<double>& theta );

/**
 * @brief Same for a mixed state given as an MPO: G(theta) = Tr(rho exp(i theta S^x_A)).
 * @param rho        Density matrix (MPO with indices s, s'; normalized to Tr(rho) = 1).
 * @param sites      Spin-1/2 site set.
 * @param block_size Number of sites of the block A (centered, see block_start).
 * @param theta      Values of theta.
 */
vector<complex<double> >
compute_generating_function( MPO* rho, const SpinHalf sites, int block_size, const vector<double>& theta );

/**
 * @brief MPO of sum_{j=start}^{start+block_size-1} S^x_j.
 * @param sites      Spin-1/2 site set.
 * @param start      First site of the block.
 * @param block_size Number of sites of the block.
 */
MPO
make_block_sx_mpo( const SpinHalf sites, int start, int block_size );

/**
 * @brief First four cumulants of S^x_A for a pure state.
 * @param psi        State; its orthogonality center is moved to the start of the block.
 * @param sites      Spin-1/2 site set.
 * @param block_size Number of sites of the block A (centered, see block_start).
 * @return {C1, C2, C3, C4} (complex in general because of the finite numerical precision).
 */
vector<complex<double> >
compute_cumulants( MPS* psi, const SpinHalf sites, int block_size );

#endif
