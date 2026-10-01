/**
 * @file purified_state.h
 * @brief Trace and observables of a density matrix stored as an unfolded (purified) MPS.
 *
 * Layout of the 2N sites: bra on sites 1..N in reversed order, ket on sites N+1..2N, so that
 * physical site q corresponds to bra site N-q+1 and ket site N+q (see dof/spin_boson.h).
 * Tr(rho O) is obtained by contracting each bra site with the corresponding ket site,
 * inserting O where needed.
 */
#ifndef MYTN_DYNAMICS_PURIFIED_STATE_H
#define MYTN_DYNAMICS_PURIFIED_STATE_H

#include <itensor/all.h>
#include <complex>

using namespace std;
using namespace itensor;

/** @return Tr(rho) of the purified state psi. */
double
compute_norm_purified_impurity(MPS* psi);

/** @return Tr(rho) of the purified state psi, for MPS with quantum numbers. */
double
compute_norm_purified_impurity_qn(MPS* psi);

/**
 * @brief Tr(rho sigma^direction_q) (spin sites) or Tr(rho n_q) (boson sites).
 * @param psi                   Purified state.
 * @param direction             "x", "y" or "z".
 * @param compute_normalization Divide by Tr(rho).
 * @param q                     Physical site; -1 (default) for all sites.
 */
vector<complex<double> >
measure_magnetization_impurity_first_site(MPS* psi , string direction, bool compute_normalization = true, int q = -1);

/**
 * @brief Tr(rho O_q) for a local operator O.
 * @param O Operator with the indices of physical site 1 (it is moved to site q).
 * @param q Physical site; -1 (default) for all sites.
 */
vector<complex<double> >
measure_local_obs_impurity_first_site(MPS *psi , const ITensor O, bool compute_normalization, int q = -1);

/**
 * @brief Tr(rho O_q1 O_q2), or its connected part.
 * @param q1, q2    Physical sites.
 * @param connected Subtract Tr(rho O_q1) Tr(rho O_q2).
 */
complex<double>
measure_correlation_impurity_first_site(MPS *psi , const ITensor O, bool compute_normalization, int q1, int q2, bool connected = false);

#endif
