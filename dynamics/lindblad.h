/**
 * @file lindblad.h
 * @brief Dissipative gates for Lindblad dynamics, d rho/dt = -i[H, rho] + sum_k gamma_k D[L_k] rho,
 *        with D[L] rho = L rho L^dag - 1/2 {L^dag L, rho}, acting on the vectorized density matrix.
 *
 * The dissipative gates are first order in dt (MyBondGateDiss stores dt * superoperator) and are
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
vector<MyBondGateDiss>
gates_local_lindblad(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt);

/** @brief Local dissipators on every site: Lj[j-1] with rate gammaj[j-1] acts on site j. */
vector<MyBondGateDiss>
gates_local_lindblad(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt);

/** @brief Dissipators with two-site jump operators L = Ti Tj (see MyTrainITensor). */
vector<MyBondGateDiss>
gates_nearest_neighbour_local_lindblad(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt);

/**
 * @brief Dissipators with jump operators acting on arbitrary groups of sites.
 * @param Lij_list Jump operators.
 * @param Lj_sites Sites of each jump operator.
 * @warning Not tested.
 */
vector<MyBondGateDiss>
gates_local_nsites_lindblad(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt);

/**
 * @brief Dissipator gamma D[L] for each L in Lj, acting on the impurity of the unfolded purified
 *        state (central bond, see models/impurity.h), to first order in dt.
 */
vector<MyBondGateDiss>
gates_dissipative_impurity(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

/** @brief As gates_dissipative_impurity, with a higher-order (Pade) approximation of the exponential. */
vector<BondGate>
gates_dissipative_impurity_high_pade(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

/**
 * @brief Copy gates built on the physical chain onto the unfolded bra-ket chain: each gate acts on
 *        the ket and, complex conjugated, on the mirrored bra.
 * @param gates_single  Gates on the physical site set.
 * @param sites_single  Physical site set (N sites).
 * @param sites_doubled Doubled site set (2N sites).
 */
vector<MyBondGate>
doubling_space_gates(const vector<BondGate> gates_single, const SiteSet sites_single, const SiteSet sites_doubled);

#endif
