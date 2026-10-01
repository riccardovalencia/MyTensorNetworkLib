/**
 * @file time_evolution.h
 * @brief One TEBD time step, for pure states and for purified density matrices.
 */
#ifndef MYTN_DYNAMICS_TIME_EVOLUTION_H
#define MYTN_DYNAMICS_TIME_EVOLUTION_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief One time step: the coherent gates and, if dissipative_gates is not empty, the dissipative
 *        gates followed by the coherent gates again (symmetric splitting; build the coherent gates
 *        with dt/2 in that case).
 *
 * @code
 * for(int k = 0 ; k < steps ; k++)
 * {
 *     psi = tebd_step(psi, gates, gates_D, {"Cutoff=", 1E-12, "MaxDim=", 64});
 *     // measure ...
 * }
 * @endcode
 *
 * @param psi               State (pure state or purified density matrix).
 * @param gates             Coherent gates, applied in order (apply_gate).
 * @param dissipative_gates Dissipative gates of the purified density matrix (apply_dissipative_gate).
 * @param args              SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved state (not normalized).
 */
MPS
tebd_step(MPS psi, const vector<TebdGate>& gates, const vector<DissipativeGate>& dissipative_gates, const Args args);

/** @brief One time step of a closed system: the coherent gates in order. */
MPS
tebd_step(MPS psi, const vector<TebdGate>& gates, const Args args);

/**
 * @brief Number of time steps dt between two measurements separated by t_measure.
 * @throws ITError if t_measure is not a positive multiple of dt (the measurement times would
 *         not be multiples of t_measure).
 */
int
compute_steps_per_measure(const double t_measure, const double dt);

#endif
