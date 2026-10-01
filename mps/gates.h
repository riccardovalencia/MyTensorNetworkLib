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
 * @brief Unitary gate exp(-i dt h) acting on one, two or three consecutive sites.
 *
 * The exponential is computed once, at construction.
 */
class TebdGate
{
private:
    ITensor gate_;
    vector<int> jn_;
    SiteSet sites_;
public:
    /**
     * @param sites Site set of the MPS the gate acts on.
     * @param j     Sites the gate acts on, e.g. {j, j+1} or {j, j+1, j+2}.
     * @param dt    Time step.
     * @param h     Local Hamiltonian term, with indices (s_j, s_j') for every site in j.
     */
	TebdGate(const SiteSet sites, vector<int> j, const double dt, const ITensor h);
    /** @return The gate exp(-i dt h). */
    ITensor gate();
    /** @return The sites the gate acts on. */
    vector<int> sites();
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
    vector<int> jnket_;
    vector<int> jnbra_;
    SiteSet sites_;
public:
    /**
     * @param sites Site set of the purified state.
     * @param jket  Sites acted on in the ket.
     * @param jbra  Sites acted on in the bra.
     * @param dt    Time step.
     * @param h     Lindblad superoperator term acting on jket and jbra.
     */
	DissipativeGate(const SiteSet sites, vector<int> jket, vector<int> jbra, const double dt, const ITensor h);
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
 * @brief Apply a one-, two- or three-site gate on consecutive sites of an MPS.
 *
 * Moves the orthogonality center to jn[0], applies the gate and splits the result back with
 * truncated SVDs. The state is not normalized.
 *
 * @code
 * for(TebdGate g : gates) psi = apply_gate(psi, g.gate(), g.sites(), {"Cutoff=",1E-12,"MaxDim=",64});
 * @endcode
 *
 * @param psi  State to evolve (taken by value).
 * @param gate Gate with unprimed/primed site indices of the sites jn.
 * @param jn   Consecutive sites the gate acts on (1 to 3 sites).
 * @param args SVD parameters; "Cutoff" and "MaxDim" are required.
 * @return The evolved MPS.
 */
MPS
apply_gate(MPS psi, const ITensor gate, const vector<int> jn, const Args args);

#endif
