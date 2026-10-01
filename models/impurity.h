/**
 * @file impurity.h
 * @brief Coherent gates of impurity models on the unfolded (purified) density matrix.
 *
 * The 2N sites of the purified state are arranged as (see dof/spin_boson.h):
 * sites 1..N hold the bra in reversed order and evolve with -H, sites N+1..2N hold the ket
 * and evolve with +H. The impurity (first physical site) sits on the central bond (N, N+1),
 * where the dissipative gates of dynamics/lindblad.h act.
 */
#ifndef MYTN_MODELS_IMPURITY_H
#define MYTN_MODELS_IMPURITY_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Coherent gates of a Kondo-like impurity model of spinful fermions.
 * @param sites Doubled Electron site set.
 * @param J     Hoppings of the physical chain (site 1 is the impurity).
 * @param hup   On-site fields for up electrons.
 * @param hdn   On-site fields for down electrons.
 * @param dt    Time step.
 */
vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

/** @brief Same model with the bath in its energy basis (inputs ordered on the physical sites 1..N). */
vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model_energy_basis(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

/**
 * @brief Coherent gates of a spin-1/2 chain (gates_spin_model conventions, models/spin_chain.h)
 *        with a dissipative impurity on the first physical site.
 * @param J     {Jxx, Jyy, Jzz}.
 * @param h     {hx, hy, hz}.
 * @param Lj    Jump operators of the impurity.
 * @param gamma Dissipation rate.
 */
vector<BondGate>
gates_coherent_part_spin_dissipative_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const double gamma, const double dt);

/** @brief As above, with next-nearest-neighbour couplings J_NNN (three-site gates). */
vector<MyBondGate>
gates_coherent_part_spin_dissipative_nnn_interactions_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> J_NNN, const vector<double> h, const double dt);

#endif
