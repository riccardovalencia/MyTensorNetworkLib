/**
 * @file lindblad.h
 * @brief Dissipative gates for Lindblad dynamics, d rho/dt = -i[H, rho] + sum_k gamma_k D[L_k] rho,
 *        with D[L] rho = L rho L^dag - 1/2 {L^dag L, rho}, acting on the vectorized density matrix.
 *
 * The dissipative gates are first order in dt (DissipativeGate stores dt * superoperator) and are
 * applied as rho -> rho + gate * rho, see dynamics/time_evolution.h.
 */
#ifndef MYTN_DYNAMICS_LINDBLAD_H
#define MYTN_DYNAMICS_LINDBLAD_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Local dissipators gamma_k D[L_k] on the listed sites.
 * @param sites    Site set of the vectorized density matrix.
 * @param Lj       Jump operators (Lj[j-1] acts on site j).
 * @param lj_sites Sites where dissipation acts.
 * @param gammaj   Rates.
 * @param dt       Time step.
 */
vector<DissipativeGate>
make_local_dissipative_gates(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt);

/** @brief Local dissipators on every site: Lj[j-1] with rate gammaj[j-1] acts on site j. */
vector<DissipativeGate>
make_local_dissipative_gates(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt);

/**
 * @brief Dissipators of operator pairs (Ti on site i, Tj on site j, rate gamma, see OperatorPair):
 *        gamma [Ti rho Tj^dag + Tj rho Ti^dag - 1/2 ({Tj^dag Ti, rho} + {Ti^dag Tj, rho})], which is
 *        gamma D[T] for i = j. Only i = j and |i - j| = 1 are implemented (ITError otherwise).
 *
 * The gates are built with dt/2 and returned as a symmetric sequence (make_symmetric_sweep).
 */
vector<DissipativeGate>
make_two_site_dissipative_gates(const SiteSet sites , vector<OperatorPair> TTrain, const double dt);

/**
 * @brief Dissipators gamma_k D[L_k] with two-site jump operators L_k (ITError for other sizes).
 *
 * The gates are built with dt/2 and returned as a symmetric sequence (make_symmetric_sweep).
 * @param Lij_list Jump operators.
 * @param Lj_sites Sites of each jump operator (two neighbouring sites).
 * @param gammaj   Rates.
 * @warning Not tested.
 */
vector<DissipativeGate>
make_multisite_dissipative_gates(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt);

/**
 * @brief Dissipator gamma D[L] for each L in Lj, acting on the impurity of the unfolded purified
 *        state (central bond, see models/impurity.h), to first order in dt.
 */
vector<DissipativeGate>
make_impurity_dissipative_gates(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

/** @brief As make_impurity_dissipative_gates, with a higher-order (Pade) approximation of the exponential. */
vector<TebdGate>
make_impurity_dissipative_gates_pade(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

/**
 * @brief Apply a dissipative gate to the vectorized density matrix, rho -> rho + gate * rho, on the
 *        two sites (j, j+1) with j = gate.ket_sites()[0], and split them back with a truncated SVD.
 *
 * The orthogonality center is not moved and the state is not normalized.
 * @code
 * for(DissipativeGate g : gates_D) psi = apply_dissipative_gate(psi, g, {"Cutoff=",1E-14,"MaxDim=",256});
 * @endcode
 * @param psi  Purified state (taken by value).
 * @param gate Dissipative gate (e.g. from make_impurity_dissipative_gates).
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The updated MPS.
 */
MPS
apply_dissipative_gate(MPS psi, DissipativeGate gate, const Args args);

/**
 * @brief Copy gates built on the physical chain onto the unfolded bra-ket chain: each gate acts on
 *        the ket and, complex conjugated, on the mirrored bra.
 *
 * Physical site i (1..N) is mapped to the ket site N + i and to the bra site N + 1 - i, both for the
 * positions of the gates and for their site indices. The ket gates come first, in the order of
 * gates_single, followed by the bra gates. swap_after() is preserved.
 * @param gates_single  Gates on the physical site set.
 * @param sites_single  Physical site set (N sites), e.g. make_spin_boson_sites(N, max_occ).
 * @param sites_doubled Doubled site set (2N sites), e.g. make_purified_spin_boson_sites(N, max_occ).
 */
vector<TebdGate>
make_purified_gates(const vector<TebdGate> gates_single, const SiteSet sites_single, const SiteSet sites_doubled);

#endif
