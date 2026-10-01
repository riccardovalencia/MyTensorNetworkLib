#ifndef MYTN_SPIN_BOSON_TEBD_H
#define MYTN_SPIN_BOSON_TEBD_H

// TEBD gates for short-range spin-1/2 models (spin chains, Rydberg arrays, PXP),
// free spinful fermions, and local Lindblad dissipators on the purified (bra-ket) state.
//
// Conventions (unless stated otherwise):
//  - sites are 1-indexed; X,Y,Z are Pauli matrices, S^a = sigma^a/2 (ITensor "Sx","Sy","Sz");
//  - n_j = (1 - Z_j)/2 = |down_z><down_z| is the Rydberg (excited-state) projector;
//  - the returned gate list implements one second-order Trotter step of length dt:
//    a forward sweep with dt/2 followed by the reversed sweep with dt/2.
//    Apply it with apply_gate (core/apply_gate.h) or ITensor's gateTEvol (for BondGate).

#include <itensor/all.h>
#include "../core/MyClasses.h"

using namespace std;
using namespace itensor;

// Nearest-neighbour spin-1/2 chain:
//   H = sum_j (hx X_j + hy Y_j + hz Z_j) + sum_j (Jxx X_j X_{j+1} + Jyy Y_j Y_{j+1} + Jzz Z_j Z_{j+1})
// J = {Jxx, Jyy, Jzz}, h = {hx, hy, hz}.
vector<MyBondGate>
gates_spin_model(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt);

// Same as gates_spin_model, returned as ITensor BondGate (for gateTEvol).
vector<BondGate>
gates_spin_model_bondgate(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt);

// gates_spin_model plus the anti-hermitian part -i/2 sum_k gamma_k L_k^dag L_k of local jump
// operators Lj acting on sites Lj_sites (effective non-hermitian Hamiltonian of quantum trajectories).
vector<BondGate>
gates_spin_eff_model_bondgate(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const vector<int> Lj_sites, const vector<double> gamma, const double dt);

// Free spinful fermions (Electron sites): hopping J_j between j and j+1 for both spin species,
// on-site fields hup_j n_up,j + hdn_j n_dn,j.
vector<BondGate>
gates_free_spinful_fermions(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

// PXP model, H = omega sum_j P_{j-1} X_j P_{j+1} with P = (1+Z)/2 (three-site gates).
vector<MyBondGate>
gates_pxp(const SiteSet sites , const double omega, const double dt);

// Rydberg chain with nearest-neighbour interactions:
//   H = sum_j Omegaj[j] S^x_j + sum_j Deltaj[j] n_j + sum_j Vj[j] n_j n_{j+1}
// Deltaj, Omegaj have N entries, Vj has N-1 (Vj[j-1] couples sites j and j+1).
vector<MyBondGate>
gates_rydberg_up_to_VNN(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

// As gates_rydberg_up_to_VNN, plus next-nearest-neighbour interactions (three-site gates).
// The NNN coupling is derived from the NN ones assuming V(r) = 1/r^6 on a line:
//   V_{j,j+2} = 1/(r_j + r_{j+1})^6 with r_j = V_j^(-1/6).
// Benchmarked against exact diagonalization (arXiv:2309.12392).
vector<MyBondGate>
gates_rydberg_up_to_VNNN(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

// Older implementation of gates_rydberg_up_to_VNNN (gates exponentiated via SVD). Prefer the one above.
vector<MyBondGate>
gates_rydberg_up_to_VNNN_deprecated(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

// Single-site field H = sum_j (w_x X_j + w_y Y_j + w_z Z_j), omegaj = {w_x, w_y, w_z}.
vector<MyBondGate>
gates_spin_local_field(const SiteSet sites , vector<double> omegaj, const double dt);

// Local Lindblad dissipators gammaj[k] D[Lj[k]] on the vectorized density matrix.
// With lj_sites the jump operators act only on the listed sites; without it, Lj[j-1] acts on site j.
vector<MyBondGateDiss>
gates_local_lindbland(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt);

vector<MyBondGateDiss>
gates_local_lindbland(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt);

// Two-site jump operators L = Ti Tj (see MyTrainITensor) with rate gamma.
vector<MyBondGateDiss>
gates_nearest_neighbour_local_lindbland(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt);

// Jump operators acting on arbitrary groups of sites Lj_sites[k]. Not tested.
vector<MyBondGateDiss>
gates_local_nsites_lindbland(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt);

// Time evolution of the purified density matrix (impurity geometry, dissipation on the central bond):
// unitary gates via gateTEvol, dissipative gates_D to first order in dt.
// Evolves for a time T from t_start; every steps_save_state steps writes the state to
// "<file_root>_psi_t<t>" and stops early if the bond dimension exceeds "MaxDim" in TEBD_args.
// normalize: rescale to Tr(rho) = 1 after each step.
MPS
TEBD_lindbland_time_evolve(MPS psi_t, vector<BondGate> gates , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

#endif
