#ifndef MYTN_SPIN_BOSON_TEBD_EDGE_DISSIPATION_H
#define MYTN_SPIN_BOSON_TEBD_EDGE_DISSIPATION_H

// Impurity problems with dissipation acting on the first physical site, simulated on the
// unfolded (purified) density matrix of 2N sites (see custom_spin_boson_doubling):
//   sites 1..N   : bra, in reversed order, evolved with -H
//   sites N+1..2N: ket, evolved with +H
// The jump / non-unitary part acts on the central bond (N, N+1).

#include <itensor/all.h>
#include "../core/MyClasses.h"

using namespace std;
using namespace itensor;

// Coherent part of the unfolded Kondo-like impurity model (spinful fermions): hoppings J,
// fields hup, hdn, given for the physical sites 1..N.
vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

// Same model in the energy basis of the bath.
vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model_energy_basis(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

// Dissipator gamma D[L] for each L in Lj, acting on the central bond, to first order in dt.
vector<MyBondGateDiss>
gates_dissipative_impurity(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

// As gates_dissipative_impurity, with a higher-order (Pade) approximation of the exponential.
vector<BondGate>
gates_dissipative_impurity_high_pade(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt);

// Coherent part of a spin-1/2 chain (spin_model conventions, see TEBD.h) with a dissipative impurity on site 1.
vector<BondGate>
gates_coherent_part_spin_dissipative_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const double gamma, const double dt);

// As above, including next-nearest-neighbour couplings J_NNN (three-site gates).
vector<MyBondGate>
gates_coherent_part_spin_dissipative_NNN_interactions_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> J_NNN, const vector<double> h, const double dt);

// Map gates built on the physical sites (sites_single) to the doubled bra-ket chain (sites_doubled):
// each gate is copied on the ket and, conjugated, on the bra.
vector<MyBondGate>
doubling_space_gates(const vector<BondGate> gates_single, const SiteSet sites_single,  const SiteSet sites_doubled);

#endif
