/**
 * @file mps_tools.h
 * @brief Model-independent manipulation of MPS: embedding states, swapping sites,
 *        density matrices and two-point functions.
 */
#ifndef MYTN_MPS_MPS_TOOLS_H
#define MYTN_MPS_MPS_TOOLS_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Overwrite a block of sites of *psi with the tensors of psi_seed.
 *
 * The sites [start, start + length(psi_seed) - 1] of *psi are replaced; their site indices
 * are kept and the link indices of the inserted block are retagged, so several copies of the
 * same state can be inserted. Used to build purified states (bra and ket halves).
 *
 * @param psi      State to modify.
 * @param psi_seed State to insert.
 * @param start    First site of *psi where psi_seed is inserted.
 * @param inverted Insert psi_seed in reversed site order.
 * @param dagger   Conjugate the tensors of psi_seed (bra half of a purified state).
 * @note For MPS with conserved quantum numbers use insert_QN_state.
 */
void
insert_state(MPS* psi, MPS psi_seed, const int start, bool inverted, bool dagger);

/**
 * @brief Same as insert_state(MPS*, MPS, int, bool, bool), for MPS with quantum numbers
 *        (link indices with a flux direction).
 */
void
insert_QN_state(MPS* psi, MPS psi_seed, const int start, bool inverted, bool dagger);

/**
 * @brief Insert a state defined on its own site set into a larger MPS.
 *
 * @param psi_t0                Target state (N sites, site set sites); modified in place.
 * @param state_to_insert       State of L sites, with site set sites_state_to_insert.
 * @param sites                 Site set of *psi_t0.
 * @param sites_state_to_insert Site set of state_to_insert.
 * @param start                 First site of *psi_t0 where the state is inserted.
 * @param L                     Number of sites of state_to_insert.
 * @param N                     Number of sites of *psi_t0.
 */
void
insert_state(MPS* psi_t0, MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int start, const int L, const int N);

/**
 * @brief Move the site at position j1 to position j2 with a sequence of SVDs.
 *
 * "Virtual" swap: the physical state is unchanged, only the order of the sites along the MPS.
 * Used to apply long-range gates, see Phys. Rev. Research 2, 043255 (2020).
 *
 * @param psi     State to modify.
 * @param j1,j2   Positions to swap (any order).
 * @param cut_off SVD truncation cutoff.
 * @param maxDim  Maximum bond dimension.
 */
void
swap_gate( MPS *psi, int j1, int j2, double cut_off, int maxDim);

/**
 * @brief Density matrix |psi><psi| of a pure state, as an MPO.
 * @param psi Pure state.
 * @return MPO with unprimed (ket) and primed (bra) site indices; bond dimension chi^2.
 */
MPO
from_MPS_to_MPDO(MPS psi);

/**
 * @brief Same as from_MPS_to_MPDO, fusing the ket and bra link indices into a single index.
 */
MPO
from_MPS_to_MPDO_v2(MPS psi);

/**
 * @brief Reduced density matrix of psi on the sites i..j (inclusive).
 * @param psi State; its orthogonality center is moved to min(i,j).
 * @param i,j First and last site (any order).
 * @return ITensor with unprimed (ket) and primed (bra) site indices of the sites i..j.
 */
ITensor
extract_reduced_density_matrix(MPS * psi, int i, int j);

/**
 * @brief Two-point function <psi| op_i op_j |psi> for operators on sites i != j.
 * @param psi   State; its orthogonality center is moved to min(i,j).
 * @param sites Site set of psi.
 * @param op_i  Operator on site i (indices s_i, s_i').
 * @param op_j  Operator on site j.
 * @param i,j   Sites (any order).
 */
complex<double>
compute_two_point( MPS *psi, const SiteSet sites, ITensor op_i, ITensor op_j, int i, int j);

#endif
