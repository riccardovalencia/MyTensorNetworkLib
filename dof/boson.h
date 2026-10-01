/**
 * @file boson.h
 * @brief Bosonic degrees of freedom (ITensor Boson sites, truncated Fock space):
 *        product and Gaussian states, local observables.
 *
 * The set_* functions need an existing MPS *psi with link indices (e.g. MPS psi = randomMPS(sites))
 * and overwrite its tensors (see set_site_tensor in mps/mps_tools.h).
 * Here sigma^x_j = a_j + a^dag_j.
 */
#ifndef MYTN_DOF_BOSON_H
#define MYTN_DOF_BOSON_H

#include <itensor/all.h>
#include <complex>
#include <vector>

using namespace std;
using namespace itensor;

/** @name Product (Fock) states */
///@{

/**
 * @brief Set one site to the Fock state |n>, leaving the other sites unchanged.
 * @param psi   State to modify (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set.
 * @param site  Site to set.
 * @param n     Occupation (n < dim of the site).
 */
void
set_site_occupation( MPS* psi, const SiteSet sites, int site, int n );

/**
 * @brief Vacuum |00...0>.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 */
void
set_vacuum_state( MPS* psi, const SiteSet sites );

/**
 * @brief One boson per site, |11...1>.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 */
void
set_unit_filling_state( MPS* psi, const SiteSet sites );

/**
 * @brief Fock state with n0 bosons on `site` and vacuum elsewhere.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 * @param n0    Number of bosons (n0 < dim of the site).
 * @param site  Occupied site.
 */
void
set_fock_excitation( MPS* psi, const SiteSet sites, int n0, int site );

/**
 * @brief Kink |1...1 0...0> with one boson on each of the first number_ones sites.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 * @param number_ones Number of occupied sites, from site 1.
 */
void
set_kink_state( MPS* psi, const SiteSet sites, int number_ones );
///@}

/** @name Coherent, squeezed and cat states */
///@{

/**
 * @brief Coherent state |alpha> (a|alpha> = alpha|alpha>) on one site, vacuum elsewhere.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 * @param site  Site of the coherent state.
 * @param alpha Amplitude (|alpha|^2 bosons on average; truncated at the Fock-space cutoff).
 */
void
set_coherent_state_on_site( MPS* psi, const SiteSet sites, int site, Cplx alpha );

/**
 * @brief Product of coherent states |alpha> on every site.
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 * @param alpha Amplitude, the same on every site.
 */
void
set_coherent_state_all_sites( MPS* psi, const SiteSet sites, Cplx alpha );

/**
 * @brief Squeezed vacuum with squeezing parameter r on one site, vacuum elsewhere (only even Fock states).
 * @param psi   State to overwrite (with link indices, e.g. MPS(sites) or randomMPS(sites)).
 * @param sites Boson site set of psi.
 * @param site  Site of the squeezed state.
 * @param r     Squeezing parameter.
 */
void
set_squeezed_state_on_site( MPS* psi, const SiteSet sites, int site, double r );

/**
 * @brief Even cat state (|alpha> + |-alpha>), normalized, on one site; vacuum elsewhere.
 * @param psi   State to overwrite (replaced by a new MPS, with bond dimension 2).
 * @param sites Boson site set.
 * @param site  Site of the cat state.
 * @param alpha Amplitude of the two coherent states.
 */
void
set_cat_state_on_site( MPS* psi, const SiteSet sites, int site, Cplx alpha );

/** @param n Non-negative integer. @return n! as a double. */
double
factorial(int n);

/**
 * @param alpha Coherent-state amplitude.
 * @param k     Fock state.
 * @return <k|alpha> = exp(-|alpha|^2/2) alpha^k / sqrt(k!).
 */
complex<double>
coherent_state_amplitude( const complex<double> alpha, const int k);

/**
 * @param r Squeezing parameter.
 * @param k Fock state (even; odd states have zero amplitude).
 * @return <k|r> = (-tanh r)^(k/2) sqrt(k!) / (2^(k/2) (k/2)! sqrt(cosh r)) for the squeezed vacuum.
 */
double
squeezed_state_amplitude( const double r, const int k);
///@}

/** @name Local observables */
///@{

/**
 * @brief Local observables of sigma^x_j = a_j + a^dag_j: this one is <(sigma^x_j)^2>.
 * @param state State (normalized); its orthogonality center is moved to j.
 * @param sites Boson site set.
 * @param j     Site.
 */
double
measure_sigma_x_squared( MPS *state , const SiteSet sites , const int j );

/** @brief <sigma^x_j> (arguments as measure_sigma_x_squared). */
double
measure_sigma_x( MPS *state , const SiteSet sites , const int j );

/** @brief <sigma^x_j n_j> (arguments as measure_sigma_x_squared). */
double
measure_sigma_x_n( MPS *state , const SiteSet sites , const int j );

/** @brief <n_j sigma^x_j> (arguments as measure_sigma_x_squared). */
double
measure_n_sigma_x( MPS *state , const SiteSet sites , const int j );

/**
 * @brief Append <n_j> for j = 1..size to occupation_number.
 * @param ground_state State; its orthogonality center is moved.
 * @param sites        Boson site set.
 * @param size         Number of sites.
 * @param occupation_number Output vector (values are appended).
 */
void
measure_occupation_number( MPS *ground_state , const SiteSet sites , const int size , vector<double> &occupation_number);

/**
 * @brief Append <n_j^2> for j = 1..size to square_occupation_number.
 * @param ground_state State; its orthogonality center is moved.
 * @param sites        Boson site set.
 * @param size         Number of sites.
 * @param square_occupation_number Output vector (values are appended).
 */
void
measure_occupation_number_squared( MPS *ground_state , const SiteSet sites , const int size , vector<double> &square_occupation_number );

/**
 * @brief Imbalance (n_k - n_max)/(n_k + n_max), with n_max the largest occupation of the other sites.
 * @param occupation_number Occupations (0-indexed). Element k is erased from the vector.
 * @param k                 Index of the reference site (0-indexed).
 */
double
compute_imbalance( vector<double> &occupation_number, int k);

/**
 * @brief Probabilities of the Fock states, <|n><n|_j> for n = 0..cut_off_fock_space, on every site.
 * @param ground_state        State; its orthogonality center is moved.
 * @param sites               Boson site set.
 * @param size                Number of sites.
 * @param cut_off_fock_space  Maximum occupation.
 * @param projector_all_sites Output: one row per site (appended).
 */
void
measure_fock_probabilities( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space , vector<vector<double> > &projector_all_sites );

/**
 * @brief Largest probability of the Fock state |cut_off - 1> over all sites, to check the
 *        Fock-space truncation.
 * @param projector_all_sites Output of measure_fock_probabilities.
 * @param size                Number of sites.
 * @param cut_off             Column (1-indexed) of projector_all_sites to check.
 */
double
compute_max_cutoff_probability(vector<vector<double> > &projector_all_sites, const int size, const int cut_off);

/**
 * @brief Connected density-density correlations <n_i n_j>_c, compared with the Gaussian
 *        (Wick) approximation built from <a>, <a^dag a>, <a a>.
 * @param psi   State; normalized in place.
 * @param sites Boson site set.
 * @param covariance_matrix_NN_system   Output: <n_i n_j>_c of the state.
 * @param covariance_matrix_NN_gaussian Output: Gaussian approximation.
 * @param relative_error                Output: relative error of the approximation.
 */
void
measure_number_covariance( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN_system, vector<vector<double> > &covariance_matrix_NN_gaussian, vector<vector<double> > &relative_error);

/**
 * @brief Variance <x^2> - <x>^2 of the quadrature x = a + a^dag of one mode.
 * @param psi  State (normalized).
 * @param A    MPO of the annihilation operator a of the mode (e.g. a dressed operator).
 * @param Adag MPO of a^dag.
 */
double
measure_variance_x(MPS *psi, MPO *A, MPO *Adag );

/** @brief Variance of the quadrature p = -i(a - a^dag) (arguments as measure_variance_x). */
double
measure_variance_p(MPS *psi, MPO *A, MPO *Adag );

/**
 * @brief Minimal quadrature variance on site j, 1 + 2<n_j> - 2|<(a^dag_j)^2>| (1 for the vacuum,
 *        < 1 for squeezed states; exact if <a_j> = 0).
 * @param psi   State (normalized); its orthogonality center is moved to j.
 * @param sites Boson site set.
 * @param j     Site.
 * @return The variance, or -1 if j is out of range.
 */
double
measure_squeezing(MPS *psi, const SiteSet sites, const int j);

/**
 * @brief measure_squeezing for a dressed mode given by MPOs: 1 + 2<N> - 2|<A^2>|.
 * @param psi   State (normalized).
 * @param sites Site set (unused, kept for symmetry with measure_squeezing).
 * @param A     MPO of the annihilation operator of the mode.
 * @param N     MPO of its number operator.
 */
double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N);
///@}

#endif
