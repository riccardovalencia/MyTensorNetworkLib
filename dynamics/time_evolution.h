/**
 * @file time_evolution.h
 * @brief Time-evolution drivers for the purified density matrix (impurity geometry, see
 *        models/impurity.h): coherent gates plus dissipation on the central bond.
 *
 * Common parameters:
 * - psi_t            : initial purified state (returned evolved);
 * - gates_D          : dissipative gates, applied to first order in dt;
 * - TEBD_args        : SVD parameters, "Cutoff" and "MaxDim" required;
 * - dissipative      : apply gates_D (otherwise unitary evolution);
 * - dt, T, t_start   : time step, total time, initial time;
 * - steps_save_state : every steps_save_state steps the state is written to "<file_root>_psi_t<t>";
 *                      the evolution stops early if the bond dimension exceeds "MaxDim";
 * - normalize        : rescale to Tr(rho) = 1 after each step.
 */
#ifndef MYTN_DYNAMICS_TIME_EVOLUTION_H
#define MYTN_DYNAMICS_TIME_EVOLUTION_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Lindblad evolution with short-range coherent gates and dissipation: each step applies
 *        the coherent gates (gateTEvol), the dissipative gates and, if dissipative, the coherent gates again.
 * @param gates Coherent gates on the purified state.
 */
MPS
TEBD_lindblad_time_evolve(MPS psi_t, vector<BondGate> gates , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

/**
 * @brief Lindblad evolution with long-range coherent gates (boson coupled to every spin, applied
 *        with swap_gate) and dissipation on the boson.
 * @param gates_H Coherent gates on the purified state.
 */
MPS
TEBD_long_range_int_lindblad_time_evolve(MPS psi_t, vector<BondGate> gates_H, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

/**
 * @brief Lindblad evolution with the coherent part given as an MPO H, applied to first order in dt.
 * @param H Hamiltonian MPO acting on the purified state.
 */
MPS
MPO_lindblad_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start = 0.);

#endif
