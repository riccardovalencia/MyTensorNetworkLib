/**
 * @file spin_chain.h
 * @brief Nearest-neighbour spin-1/2 chains: TEBD gates.
 *
 * X, Y, Z are Pauli matrices. Unless stated otherwise the gate lists implement one
 * second-order Trotter step of length dt: a forward sweep with dt/2 followed by the reversed
 * sweep. Apply them with tebd_step (dynamics/time_evolution.h) or apply_gates (mps/gates.h).
 */
#ifndef MYTN_MODELS_SPIN_CHAIN_H
#define MYTN_MODELS_SPIN_CHAIN_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Two-site term on the bond (j, j+1) of the spin chain of make_spin_chain_gates.
 *
 * The bond couplings enter in full; the fields of sites j and j+1 are divided by the number of
 * gates sharing them (count_left, count_right; see count_gates_containing in mps/gates.h), so
 * that the sum of the bond terms is H.
 * @param sites       Spin-1/2 site set.
 * @param J           Couplings {Jxx, Jyy, Jzz}.
 * @param h           Fields {hx, hy, hz}.
 * @param j           Left site of the bond.
 * @param count_left  Number of bonds sharing the field of site j.
 * @param count_right Number of bonds sharing the field of site j+1.
 */
ITensor
make_spin_chain_bond_hamiltonian(const SiteSet sites, const vector<double> J, const vector<double> h, const int j, const int count_left, const int count_right);

/**
 * @brief Gates of the nearest-neighbour spin chain
 *        H = sum_j (hx X_j + hy Y_j + hz Z_j) + sum_j (Jxx X_j X_{j+1} + Jyy Y_j Y_{j+1} + Jzz Z_j Z_{j+1}).
 * @param sites Spin-1/2 site set.
 * @param J     Couplings {Jxx, Jyy, Jzz}.
 * @param h     Fields {hx, hy, hz}.
 * @param dt    Time step.
 */
vector<TebdGate>
make_spin_chain_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt);

/**
 * @brief Gates of make_spin_chain_gates plus the anti-hermitian term -i/2 sum_k gamma_k L_k^dag L_k
 *        (effective non-hermitian Hamiltonian of quantum trajectories).
 * @param Lj       Local jump operators.
 * @param Lj_sites Site of each jump operator.
 * @param gamma    Rate of each jump operator.
 */
vector<TebdGate>
make_spin_chain_effective_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const vector<int> Lj_sites, const vector<double> gamma, const double dt);

/**
 * @brief Single-site gates of H = sum_j (w_x X_j + w_y Y_j + w_z Z_j).
 * @param omegaj Field {w_x, w_y, w_z}, the same on every site.
 */
vector<TebdGate>
make_local_field_gates(const SiteSet sites , vector<double> omegaj, const double dt);

/**
 * @brief Gates of the Ising chain in longitudinal (hx) and transverse (hz) fields,
 *        H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ]
 *        (make_spin_chain_gates with J = {-J, 0, 0} and h = {-J hx, 0, -J hz}).
 * @param sites Spin-1/2 site set.
 * @param J     Overall energy scale.
 * @param hx    Longitudinal field (along the Ising axis x).
 * @param hz    Transverse field.
 * @param dt    Time step.
 */
vector<TebdGate>
make_ising_gates( const SiteSet sites , const double J , const double hx , const double hz , const double dt );

#endif
