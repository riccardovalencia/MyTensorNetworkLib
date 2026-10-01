/**
 * @file light_matter.h
 * @brief Collective light-matter models: one bosonic mode a (site 1, see make_spin_boson_sites)
 *        coupled to all the spin-1/2 (sites 2..N). S^a = sum_j s^a_j are collective spins.
 *
 * The all-to-all coupling is split into "long-range" gates between the boson and each spin,
 * applied with swap_sites (mps/mps_tools.h); see the collective_light_matter examples.
 */
#ifndef MYTN_MODELS_LIGHT_MATTER_H
#define MYTN_MODELS_LIGHT_MATTER_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Gates of the Tavis-Cummings model H = omega0 n_a + h S^z + g (S^+ a + S^- a^dag).
 * @param sites            Site set from make_spin_boson_sites.
 * @param omega0           Frequency of the boson.
 * @param h                Field on the spins.
 * @param g                Light-matter coupling.
 * @param dt               Time step.
 * @param photon_or_matter "short-range" for the local terms; otherwise the photon-matter gates.
 */
vector<BondGate>
make_tavis_cummings_gates(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string photon_or_matter);

/**
 * @brief Gates of a general light-matter model,
 *        H = omega0 n_a + h S^z + (coupling) + (nearest-neighbour spin interaction).
 *
 * @param sites            Site set from make_spin_boson_sites.
 * @param omega0           Frequency of the boson.
 * @param h                Field on the spins.
 * @param g                Light-matter coupling (include the 1/sqrt(N) scaling if wanted).
 * @param dt               Time step.
 * @param matter_or_photon "short-range" for the local terms, "long-range" for the photon-matter gates.
 * @param type_of_coupling "dicke": g S^x (a + a^dag); "tavis" (default): g (S^+ a + S^- a^dag).
 * @param V                Strength of the spin-spin interaction (default 0).
 * @param interaction_axis "z" (default): V sum_j n_j n_{j+1}, n = (1-Z)/2 (Rydberg interaction);
 *                         "x": V sum_j X_j X_{j+1}.
 */
vector<BondGate>
make_light_matter_gates(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling = "tavis", const double V = 0, string interaction_axis = "z");

/** @name Purified (bra-ket) states
 *  Gates mapped onto the doubled chain with make_purified_gates (dynamics/lindblad.h): the two
 *  bosons sit on the central bond, the spins of the ket to their right and those of the bra
 *  (mirrored) to their left.
 */
///@{

/**
 * @brief Apply a gate acting on two neighbouring sites of the purified chain, given in any order.
 *
 * Moves the orthogonality center to the left site, applies the gate and splits back with a
 * truncated SVD.
 * @param psi  Purified state (taken by value).
 * @param gate Gate on two neighbouring sites.
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 */
MPS
apply_local_gate_purified(MPS psi, TebdGate gate, const Args args);

/**
 * @brief Apply a photon-matter gate and move the boson one site outward with a swap.
 *
 * The gate couples the boson, adjacent to the spin it acts on, and that spin. Applying the
 * "long-range" gates of make_light_matter_gates in order lets the boson travel through the ket (or the
 * bra) and interact with every spin.
 * @param psi  Purified state (taken by value).
 * @param gate Photon-matter gate mapped onto the doubled chain; jn = {boson site, spin site}.
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 */
MPS
apply_photon_matter_gate_purified(MPS psi, TebdGate gate, const Args args);
///@}

#endif
