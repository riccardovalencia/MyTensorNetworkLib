/**
 * @file gates.cc
 * @brief Implementation of gates.h (interfaces documented in the header, logic commented here).
 */
#include "gates.h"
#include "mps_tools.h"
#include <algorithm>
#include <itensor/all.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// gate containers: the constructors store the data, the accessors return it

// gate = exp(-i dt h), exponentiating the hermitian h on its pairs of indices (s, s')
TebdGate::TebdGate(vector<int> sites, const double dt, const ITensor h)
{
    sites_      = sites;
    gate_       = expHermitian(h, -1_i * dt);
    swap_after_ = false;
}

// gate stored as given (e.g. a swap gate, or a gate exponentiated elsewhere)
TebdGate::TebdGate(vector<int> sites, const ITensor gate, bool swap_after)
{
    sites_      = sites;
    gate_       = gate;
    swap_after_ = swap_after;
}

ITensor TebdGate::gate()
{
    return gate_;
}

vector<int> TebdGate::sites()
{
    return sites_;
}

bool TebdGate::swap_after()
{
    return swap_after_;
}

void TebdGate::set_gate(ITensor new_gate)
{
    gate_ = new_gate;
}


// first-order term dt h of exp(dt h): tebd_step applies psi + dt h psi
DissipativeGate::DissipativeGate(vector<int> jket, vector<int> jbra, double dt, ITensor h)
{
	ket_sites_ = jket;
	bra_sites_ = jbra;

	gate_ = h * dt;
}


ITensor DissipativeGate::gate()
{
	return gate_;
}


vector<int> DissipativeGate::ket_sites()
{
	return ket_sites_;
}


vector<int> DissipativeGate::bra_sites()
{
	return bra_sites_;
}


OperatorPair::OperatorPair(vector<ITensor> T, vector<int> j, double gamma)
{
	Ti_   = T[0];
    Tj_   = T[1];
    i_    = j[0];
    j_    = j[1];
	gamma_ = gamma; 
}


ITensor OperatorPair::op_i()
{
	return Ti_;
}


ITensor OperatorPair::op_j()
{
	return Tj_;
}


int OperatorPair::site_i()
{
	return i_;
}


int OperatorPair::site_j()
{
	return j_;
}


double OperatorPair::rate()
{
	return gamma_;
}


// ----------------------------------------------------------
// gate application

// Contract the gate with the tensors of its (sorted, consecutive) sites, with the orthogonality
// center on the first one, then split the result back into one tensor per site: one truncated SVD
// for two sites, two successive SVDs (left site first) for three. The singular values are
// absorbed to the right, so the orthogonality center ends on the last site.
MPS
apply_gate(MPS psi, const ITensor gate, vector<int> sites, const Args args)
{
    double cut_off = args.getReal("Cutoff");
    int max_dim    = args.getInt("MaxDim");
    sort(sites.begin(), sites.end());
    int j = sites[0];
    if(j < 1 || sites.back() > length(psi) || sites.back() - j != int(sites.size()) - 1)
        throw ITError("apply_gate: the sites must be consecutive sites of the MPS");

    psi.position(j);
    ITensor AA = gate;
    for(int q : sites) AA *= psi(q);
    AA.mapPrime(1,0);

    if(sites.size() == 1)
    {
        psi.set(j, AA);
    }
    else if(sites.size() == 2)
    {
        auto [U,S,V] = svd(AA, inds(psi(j)), {"Cutoff=", cut_off, "MaxDim=", max_dim});
        psi.set(j, U);
        psi.set(j+1, S*V);
    }
    else if(sites.size() == 3)
    {
        auto [U,S,V] = svd(AA, inds(psi(j)), {"Cutoff=", cut_off, "MaxDim=", max_dim});
        Index l = commonIndex(U, S);
        psi.set(j, U);
        auto [U1,S1,V1] = svd(S*V, {inds(psi(j+1)), l}, {"Cutoff=", cut_off, "MaxDim=", max_dim});
        psi.set(j+1, U1);
        psi.set(j+2, S1*V1);
    }
    else throw ITError("apply_gate: gates on more than 3 sites are not implemented");

    return psi;
}


// apply the stored gate, then exchange its first two sites if the gate asks for it
MPS
apply_gate(MPS psi, TebdGate gate, const Args args)
{
    psi = apply_gate(psi, gate.gate(), gate.sites(), args);
    if(gate.swap_after())
    {
        vector<int> sites = gate.sites();
        int j = *min_element(sites.begin(), sites.end());
        swap_sites(&psi, j, j+1, args.getReal("Cutoff"), args.getInt("MaxDim"));
    }
    return psi;
}


// gates applied one after the other, in the order of the list
MPS
apply_gates(MPS psi, vector<TebdGate> gates, const Args args)
{
    for(TebdGate gate : gates) psi = apply_gate(psi, gate, args);
    return psi;
}


// The gate starting at site g (1 <= g <= N - gate_size + 1) contains the term if
// g <= first and g + gate_size - 1 >= first + term_size - 1: count the g in both ranges.
int
count_gates_containing(int first, int term_size, int gate_size, int N)
{
    int first_gate = max(1, first + term_size - gate_size);
    int last_gate  = min(first, N - gate_size + 1);
    return last_gate - first_gate + 1;
}
