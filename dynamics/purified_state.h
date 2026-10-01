/**
 * @file purified_state.h
 * @brief Trace and observables of a density matrix stored as an unfolded (purified) MPS.
 *
 * Layout of the 2N sites: bra on sites 1..N in reversed order, ket on sites N+1..2N, so that
 * physical site q corresponds to bra site N-q+1 and ket site N+q (see dof/spin_boson.h).
 * Tr(rho O) is obtained by contracting each bra site with the corresponding ket site,
 * inserting O where needed. A local operator O with indices (s, s') enters Tr(rho O) =
 * sum_{ab} rho_ab <b|O|a>: its input index s is contracted with the ket, s' with the bra.
 */
#ifndef MYTN_DYNAMICS_PURIFIED_STATE_H
#define MYTN_DYNAMICS_PURIFIED_STATE_H

#include <itensor/all.h>
#include <complex>

using namespace std;
using namespace itensor;

/** @return Tr(rho) of the purified state psi. */
double
compute_trace_purified(MPS* psi);

/**
 * @brief Tr(rho sigma^direction_q) with Pauli matrices (0 on sites that are not spins).
 * @param psi                   Purified state.
 * @param direction             "x", "y" or "z".
 * @param compute_normalization Divide by Tr(rho).
 * @param q                     Physical site; -1 (default) for all sites.
 * @return One value per physical site (q = -1) or the single value on site q.
 * @throws ITError for a direction other than "x", "y", "z".
 */
vector<complex<double> >
measure_magnetization_purified(MPS* psi , string direction, bool compute_normalization = true, int q = -1);

/**
 * @brief Tr(rho O_q) for a local operator O.
 * @param O Operator with the indices of physical site 1 (it is moved to site q).
 * @param q Physical site; -1 (default) for all sites.
 */
vector<complex<double> >
measure_local_operator_purified(MPS *psi , const ITensor O, bool compute_normalization, int q = -1);

/**
 * @brief Tr(rho O_q1 O_q2), or its connected part.
 * @param q1, q2    Different physical sites.
 * @param connected Subtract Tr(rho O_q1) Tr(rho O_q2).
 * @throws ITError if q1 == q2 or a site is outside 1..N.
 */
complex<double>
measure_correlation_purified(MPS *psi , const ITensor O, bool compute_normalization, int q1, int q2, bool connected = false);

#endif
