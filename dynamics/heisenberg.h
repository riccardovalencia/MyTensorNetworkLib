/**
 * @file heisenberg.h
 * @brief Time evolution of operators (Heisenberg picture), O -> U^dag O U, for operators stored as MPOs.
 *
 * MPO convention (as ITensor op): an operator tensor has input indices s (prime level 0) and output
 * indices s' (prime level 1). One step uses either the TEBD gates of a model (second order, the same
 * lists as in the Schroedinger picture, see dynamics/time_evolution.h) or a propagator given as an MPO
 * (e.g. ITensor toExpH). Then <psi| O(t) |psi> = <psi(t)| O |psi(t)>.
 *
 * For a time-dependent Hamiltonian, U(t) = U_n ... U_1 and O(t) = U_1^dag ... U_n^dag O U_n ... U_1:
 * apply the steps from the last to the first. Applying them from the first to the last with the
 * inverse steps gives instead the dressed operator U O U^dag (see make_dressed_operator).
 */
#ifndef MYTN_DYNAMICS_HEISENBERG_H
#define MYTN_DYNAMICS_HEISENBERG_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief Adjoint U^dag of an operator MPO (complex conjugate, input and output site indices exchanged).
 * @param U Operator MPO.
 */
MPO
make_adjoint_mpo(const MPO& U);

/**
 * @brief One gate in the Heisenberg picture, O -> g^dag O g, on the sites of the gate.
 *
 * The tensors of the gate sites are contracted with g and g^dag and split back with truncated SVDs
 * (orthogonality center first moved to the first site). swap_after() gates also exchange the two
 * sites of the operator, as apply_gate does for states.
 * @param O    Operator.
 * @param gate One-, two- or three-site gate.
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved operator.
 */
MPO
apply_gate_to_operator(MPO O, TebdGate gate, const Args args);

/**
 * @brief One time step in the Heisenberg picture, O -> U^dag O U, with U the product of the gates
 *        as applied to a state (first gate first): the gates are applied to O in reverse order.
 * @code
 * vector<TebdGate> gates = make_spin_chain_gates(sites, J, h, dt);
 * for(int k = 0 ; k < steps ; k++) O = heisenberg_step(O, gates, {"Cutoff=", 1E-12, "MaxDim=", 256});
 * Cplx z = innerC(psi0, O, psi0);   // = <psi(t)| O |psi(t)>
 * @endcode
 * @param O     Operator.
 * @param gates Gates of one time step (e.g. a symmetric sweep).
 * @param args  SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved operator.
 */
MPO
heisenberg_step(MPO O, const vector<TebdGate>& gates, const Args args);

/**
 * @brief One time step in the Heisenberg picture with a propagator MPO U: O -> U^dag O U.
 * @param O    Operator.
 * @param U    Propagator of one step (e.g. toExpH(ampo, Cplx_i * dt) = exp(-i dt H) to first order).
 * @param args Truncation of the MPO products (nmultMPO), e.g. "Cutoff", "MaxDim".
 * @return The evolved operator.
 */
MPO
heisenberg_step(MPO O, const MPO& U, const Args args);

#endif
