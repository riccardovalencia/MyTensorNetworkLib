/**
 * @file gates.h
 * @brief Gate containers for TEBD and their application to an MPS.
 *
 * Sites are 1-indexed, as in ITensor. A gate acting on n sites carries the unprimed (input)
 * and primed (output) site indices of those sites.
 */
#ifndef MYTN_MPS_GATES_H
#define MYTN_MPS_GATES_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Gate acting on one, two or three consecutive sites of an MPS.
 *
 * Built either from a local Hamiltonian term h, as exp(-i dt h) (h hermitian), or from an already
 * exponentiated gate. If swap_after() is true, apply_gate swaps the first two sites after applying
 * the gate (used to move a site, e.g. a cavity mode, along the chain, see models/light_matter.h).
 */
class TebdGate
{
private:
    ITensor gate_;
    vector<int> sites_;
    bool swap_after_;
public:
    /**
     * @param sites Sites the gate acts on, e.g. {j, j+1} or {j, j+1, j+2}.
     * @param dt    Time step.
     * @param h     Local hermitian Hamiltonian term, with indices (s_j, s_j') for every site in sites.
     */
    TebdGate(vector<int> sites, const double dt, const ITensor h);
    /**
     * @param sites      Sites the gate acts on.
     * @param gate       Gate with unprimed (input) and primed (output) site indices.
     * @param swap_after Swap the first two sites after applying the gate.
     */
    TebdGate(vector<int> sites, const ITensor gate, bool swap_after = false);
    /** @return The gate. */
    ITensor gate();
    /** @return The sites the gate acts on. */
    vector<int> sites();
    /** @return Whether the first two sites are swapped after the gate. */
    bool swap_after();
    /** @brief Replace the stored gate, e.g. after mapping it onto different site indices. */
    void set_gate(ITensor new_gate);
};

/**
 * @brief Dissipative gate for the vectorized (purified) density matrix.
 *
 * Stores the first-order term dt*h of exp(dt*h), acting on the sites jket of the ket and jbra
 * of the bra. Apply it as psi -> psi + gate*psi (see dynamics/time_evolution.h).
 */
class DissipativeGate
{
private:
    ITensor gate_;
    vector<int> ket_sites_;
    vector<int> bra_sites_;
public:
    /**
     * @param jket  Sites acted on in the ket.
     * @param jbra  Sites acted on in the bra.
     * @param dt    Time step.
     * @param h     Lindblad superoperator term acting on jket and jbra.
     */
    DissipativeGate(vector<int> jket, vector<int> jbra, const double dt, const ITensor h);
    /** @return dt*h (linear approximation of the exponential). */
    ITensor gate();
    /** @return The sites acted on in the ket. */
    vector<int> ket_sites();
    /** @return The sites acted on in the bra. */
    vector<int> bra_sites();
};

/**
 * @brief Pair of local operators (Ti on site i, Tj on site j) with a rate gamma.
 *
 * Describes a two-site jump operator L = Ti Tj; see make_two_site_dissipative_gates.
 */
class OperatorPair
{
private:
    ITensor Ti_;
    ITensor Tj_;
    int i_;
    int j_;
    double gamma_;

public:
    /**
     * @param T     Operators {Ti, Tj}.
     * @param j     Sites {i, j}.
     * @param gamma Rate of the jump operator.
     */
	OperatorPair(vector<ITensor> T, vector<int> j, double gamma);
    ITensor op_i();     ///< Operator on site i.
    ITensor op_j();     ///< Operator on site j.
    int site_i();       ///< Site i.
    int site_j();       ///< Site j.
    double rate();      ///< Rate.
};

/**
 * @brief Symmetric (second-order Trotter) sequence: the gates followed by the same gates in reverse
 *        order. Build the gates with half the time step.
 * @param gates Gates of one sweep (TebdGate or DissipativeGate), each with dt/2.
 * @return The 2 * gates.size() gates of one step dt.
 * @code
 * vector<TebdGate> half_step;
 * for(int j = 1 ; j < N ; j++) half_step.push_back(TebdGate({j,j+1}, dt/2., h[j]));
 * vector<TebdGate> gates = make_symmetric_sweep(half_step);   // one step dt
 * @endcode
 */
template <class Gate>
vector<Gate>
make_symmetric_sweep(vector<Gate> gates)
{
    vector<Gate> reversed(gates.rbegin(), gates.rend());
    gates.insert(gates.end(), reversed.begin(), reversed.end());
    return gates;
}

/**
 * @brief Number of gates on gate_size consecutive sites (starting at 1, 2, ..., N - gate_size + 1)
 *        containing the term on the term_size sites first, ..., first + term_size - 1.
 *
 * When the terms of a Hamiltonian are distributed over overlapping gates, each gate takes the
 * fraction 1/count of the term. E.g. for two-site gates, an on-site term (term_size = 1) is shared
 * by 2 gates in the bulk and belongs to a single gate on the edges.
 * @param first     First site of the term.
 * @param term_size Number of sites of the term.
 * @param gate_size Number of sites of each gate.
 * @param N         Number of sites of the chain.
 * @return The number of gates containing all the sites of the term.
 */
int
count_gates_containing(int first, int term_size, int gate_size, int N);

/**
 * @brief Apply a one-, two- or three-site gate on consecutive sites of an MPS.
 *
 * Moves the orthogonality center to the first site, applies the gate and splits the result back
 * with truncated SVDs. The state is not normalized.
 *
 * @param psi   State to evolve (taken by value).
 * @param gate  Gate with unprimed/primed site indices of the sites it acts on.
 * @param sites Consecutive sites the gate acts on (1 to 3 sites, in any order).
 * @param args  SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved MPS.
 */
MPS
apply_gate(MPS psi, const ITensor gate, vector<int> sites, const Args args);

/**
 * @brief Apply a TebdGate (followed by a swap of its first two sites if gate.swap_after()).
 * @param psi  State to evolve (taken by value).
 * @param gate Gate and the sites it acts on.
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved MPS (not normalized).
 */
MPS
apply_gate(MPS psi, TebdGate gate, const Args args);

/**
 * @brief Apply a list of gates in order.
 * @param psi   State to evolve (taken by value).
 * @param gates Gates, applied from first to last.
 * @param args  SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved MPS (not normalized).
 * @code
 * psi = apply_gates(psi, make_rydberg_gates_nnn(sites, Delta, Omega, V, dt), {"Cutoff=",1E-12,"MaxDim=",64});
 * @endcode
 */
MPS
apply_gates(MPS psi, vector<TebdGate> gates, const Args args);

#endif
