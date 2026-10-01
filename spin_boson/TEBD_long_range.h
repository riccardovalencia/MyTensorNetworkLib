#ifndef MYTN_SPIN_BOSON_TEBD_LONG_RANGE_H
#define MYTN_SPIN_BOSON_TEBD_LONG_RANGE_H

// Collective light-matter models: one bosonic mode (site 1, see custom_spin_boson) coupled
// to all the spin-1/2 (sites 2..N). The all-to-all coupling is applied with swap gates.

#include <itensor/all.h>
#include "../core/MyClasses.h"

using namespace std;
using namespace itensor;

// Tavis-Cummings model H = omega0 n_a + h S^z + g (S^+ a + S^- a^dag), S^a = sum_j s^a_j.
// photon_or_matter = "short-range": local terms only; otherwise the long-range photon-matter gates.
vector<BondGate>
gates_tavis_cummings(const SiteSet sites , const double omega0 , const double h , const double g,const double dt,string photon_or_matter);

// General light-matter model, H = omega0 n_a + h S^z + coupling + nearest-neighbour spin term, with
//   type_of_coupling = "dicke": g S^x (a + a^dag);  "tavis": g (S^+ a + S^- a^dag);
//   interaction_axis = "z" (default): V sum_j n_j n_{j+1}, n = (1-Z)/2 (Rydberg interaction);
//                      "x": V sum_j X_j X_{j+1}.
// matter_or_photon = "short-range" returns the local gates, "long-range" the photon-matter gates
// (to be applied together with swap_gate, see examples/collective_light_matter_*).
vector<BondGate>
gates_photon_matter(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling = "tavis", const double V = 0, string interaction_axis="z");

// Lindblad evolution of the purified state with a coherent MPO H applied to first order in dt
// and local dissipation on the central (photon) bond.
MPS
MPO_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

// Lindblad evolution of the purified state with long-range coherent gates (applied with swap gates)
// and local dissipation on the central (photon) bond.
MPS
TEBD_long_range_int_lindbland_time_evolve(MPS psi_t, vector<BondGate> gates_H, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

#endif
