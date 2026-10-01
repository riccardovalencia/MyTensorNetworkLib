/**
 * @file boson.h
 * @brief Bosonic degrees of freedom (ITensor Boson sites, truncated Fock space):
 *        product and Gaussian states, local observables.
 *
 * Unless stated otherwise, the functions writing a state need an existing MPS *psi with link
 * indices (e.g. MPS psi = randomMPS(sites)) and overwrite its tensors.
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
 * @brief Fock state with n0 bosons on site excitation_position and vacuum elsewhere.
 * @param psi                 State to overwrite.
 * @param sites               Boson site set.
 * @param size                Number of sites.
 * @param n0                  Occupation of the excited site.
 * @param excitation_position Site with n0 bosons.
 */
void
set_fock_excitation( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position);

/**
 * @brief Set site excitation_position to the Fock state |n0>, leaving the other sites unchanged.
 * @param psi, sites, size, n0, excitation_position As in set_fock_excitation.
 */
void
set_fock_excitation_pinned( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position);

/** @brief Vacuum |00...0>, keeping the link indices of *psi. */
void
set_vacuum_state( MPS* psi, const SiteSet sites, const int size );

/** @brief One boson per site, |11...1>, keeping the link indices of *psi. */
void
set_unit_filling_state( MPS* psi, const SiteSet sites, const int size );

/**
 * @brief Kink |1...1 0...0> with one boson on each of the first number_ones sites.
 * @param psi         State to overwrite.
 * @param sites       Boson site set.
 * @param number_ones Number of occupied sites.
 */
void
set_kink_state(MPS *psi, const SiteSet sites, const int number_ones);

/**
 * @brief Set site position to the Fock state |n>, leaving the other sites unchanged.
 * @param psi      State to modify.
 * @param sites    Boson site set.
 * @param position Site to set.
 * @param n        Occupation (n < dim of the site).
 */
void
set_site_occupation(MPS *psi, const SiteSet sites, const int position, const int n);
///@}

/** @name Coherent, squeezed and cat states */
///@{

/**
 * @brief Coherent state |alpha> (a|alpha> = alpha|alpha>) on one site, vacuum elsewhere.
 * @param psi   State to overwrite.
 * @param sites Boson site set.
 * @param size  Number of sites.
 * @param site  Site hosting the coherent state.
 * @param alpha Amplitude.
 */
void
set_coherent_state_on_site( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha);

/** @brief Product of coherent states |alpha> on every site. */
void
set_coherent_state_all_sites( MPS* psi, const SiteSet sites, const int size , const complex<double> alpha );

/**
 * @brief Squeezed vacuum with squeezing parameter r on one site, vacuum elsewhere.
 * @param r Squeezing parameter (only even Fock states are populated).
 */
void
set_squeezed_state_on_site( MPS* psi, const SiteSet sites, const int size , const int site, const double r );

/** @brief Even cat state (|alpha> + |-alpha>), normalized, on one site; vacuum elsewhere. */
void
set_cat_state_on_site( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha );

/** @return n! as a double. */
double
factorial(int n);

/** @return Fock amplitude <k|alpha> = exp(-|alpha|^2/2) alpha^k / sqrt(k!) of a coherent state. */
complex<double>
coherent_state_amplitude( const complex<double> alpha, const int k);

/** @return Fock amplitude <k|r> of the squeezed vacuum with squeezing parameter r (k even). */
double
squeezed_state_amplitude( const double r, const int k);
///@}

/** @name Local observables */
///@{

/** @return <(sigma^x_j)^2>, sigma^x = a + a^dag. The orthogonality center is moved to j. */
double
measure_sigma_x_squared( MPS *state , const SiteSet sites , const int j );

/** @return <sigma^x_j>. */
double
measure_sigma_x( MPS *state , const SiteSet sites , const int j );

/** @return <sigma^x_j n_j>. */
double
measure_sigma_x_n( MPS *state , const SiteSet sites , const int j );

/** @return <n_j sigma^x_j>. */
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

/** @brief Append <n_j^2> for j = 1..size to square_occupation_number. */
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
 * @param occupation_number   Occupations <n_j>, used only to print a consistency check.
 */
void
measure_fock_probabilities( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space , vector<vector<double> > &projector_all_sites , vector<double> &occupation_number);

/**
 * @return The largest probability of the Fock state |cut_off - 1> over all sites
 *         (to check the Fock-space truncation), from the output of measure_fock_probabilities.
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
 * @return Variance of x = a + a^dag, computed from the MPOs of a and a^dag of one mode.
 */
double
measure_variance_x(MPS *psi, MPO *A, MPO *Adag );

/**
 * @return Variance of p = -i(a - a^dag), computed from the MPOs of a and a^dag of one mode.
 */
double
measure_variance_p(MPS *psi, MPO *A, MPO *Adag );

/**
 * @return Minimal quadrature variance on site j, 1 + 2<n_j> - 2|<(a^dag_j)^2>| (1 for the vacuum,
 *         < 1 for squeezed states), or -1 if j is out of range.
 */
double
measure_squeezing(MPS *psi, const SiteSet sites, const int j);

/**
 * @return Same as measure_squeezing with dressed operators given as MPOs: A (annihilation) and N (number).
 */
double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N);
///@}

#endif
