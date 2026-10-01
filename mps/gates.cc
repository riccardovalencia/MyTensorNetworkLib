/**
 * @file gates.cc
 * @brief Implementation of gates.h (the functions are documented in the header).
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


TebdGate::TebdGate(vector<int> sites, const double dt, const ITensor h)
{
    sites_      = sites;
    gate_       = expHermitian(h, -1_i * dt);
    swap_after_ = false;
}

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


DissipativeGate::DissipativeGate(vector<int> jket, vector<int> jbra, double dt, ITensor h)
{
	ket_sites_ = jket;
	bra_sites_ = jbra;

	// linear approximation - tested and works well for our purposes
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
// Apply a gate on the jn sites of an MPS psi 
// It is possible to apply up to 3-sites gates.

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


MPS
apply_gates(MPS psi, vector<TebdGate> gates, const Args args)
{
    for(TebdGate gate : gates) psi = apply_gate(psi, gate, args);
    return psi;
}


int
count_gates_containing(int first, int term_size, int gate_size, int N)
{
    int first_gate = max(1, first + term_size - gate_size);
    int last_gate  = min(first, N - gate_size + 1);
    return last_gate - first_gate + 1;
}
