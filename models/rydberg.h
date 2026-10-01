/**
 * @file rydberg.h
 * @brief Rydberg atom arrays and the PXP model: TEBD gates, MPO, interactions from positions.
 *
 * Each atom is a spin-1/2 with |up_z> the ground state and |down_z> the Rydberg state;
 * n_j = (1 - Z_j)/2 and S^x = X/2. Gate lists implement one second-order Trotter step of
 * length dt (forward sweep with dt/2, then the reversed sweep); apply them with apply_gate.
 */
#ifndef MYTN_MODELS_RYDBERG_H
#define MYTN_MODELS_RYDBERG_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Gates of the Rydberg chain with nearest-neighbour interactions,
 *        H = sum_j Omega_j S^x_j + sum_j Delta_j n_j + sum_j V_j n_j n_{j+1}.
 * @param sites  Spin-1/2 site set of N sites.
 * @param Deltaj Detunings (N values).
 * @param Omegaj Rabi frequencies (N values).
 * @param Vj     Interactions (N-1 values); Vj[j-1] couples sites j and j+1.
 * @param dt     Time step.
 */
vector<TebdGate>
make_rydberg_gates_nn(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

/**
 * @brief As make_rydberg_gates_nn, plus next-nearest-neighbour interactions (three-site gates).
 *
 * The next-nearest-neighbour coupling is derived from the nearest-neighbour ones assuming
 * V(r) = 1/r^6 on a line: V_{j,j+2} = 1/(r_j + r_{j+1})^6 with r_j = V_j^(-1/6).
 * Benchmarked against exact diagonalization (arXiv:2309.12392).
 */
vector<TebdGate>
make_rydberg_gates_nnn(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

/** @deprecated Older implementation of make_rydberg_gates_nnn; use make_rydberg_gates_nnn. */
vector<TebdGate>
make_rydberg_gates_nnn_deprecated(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt);

/**
 * @brief Three-site gates of the PXP model, H = omega sum_j P_{j-1} X_j P_{j+1}, P = (1+Z)/2.
 * @param sites Spin-1/2 site set.
 * @param omega Rabi frequency.
 * @param dt    Time step.
 */
vector<TebdGate>
make_pxp_gates(const SiteSet sites , const double omega, const double dt);

/** @brief MPO of the PXP Hamiltonian H = omega sum_j P_{j-1} X_j P_{j+1}. */
MPO
make_pxp_mpo(const SiteSet s, const double omega);

/**
 * @brief Nearest-neighbour couplings V_j = 1/|r_j - r_{j+1}|^alpha from atomic positions.
 * @param rj    Positions {x, y, z} of the N atoms, in order along the chain.
 * @param alpha Power-law exponent (6 for van der Waals interactions).
 * @return N-1 couplings.
 */
vector<double>
compute_power_law_couplings(const vector< vector<double> > rj , const double alpha);

#endif
