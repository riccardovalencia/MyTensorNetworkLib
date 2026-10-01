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
vector<TebdGate>
make_kondo_impurity_gates(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

/**
 * @brief Coherent gates of a spin-1/2 chain (make_spin_chain_gates conventions, models/spin_chain.h)
 *        with a dissipative impurity on the first physical site.
 *
 * The dissipation of the impurity is added with make_impurity_dissipative_gates (dynamics/lindblad.h).
 * @param J  {Jxx, Jyy, Jzz}.
 * @param h  {hx, hy, hz}.
 * @param dt Time step.
 */
vector<TebdGate>
make_spin_impurity_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt);

/**
 * @brief As make_spin_impurity_gates, with next-nearest-neighbour couplings
 *        J_NNN = {Jxx', Jyy', Jzz'} (three-site gates; at least 3 physical sites).
 */
vector<TebdGate>
make_spin_impurity_nnn_gates(const SiteSet sites , const vector<double> J, const vector<double> J_NNN, const vector<double> h, const double dt);

#endif
