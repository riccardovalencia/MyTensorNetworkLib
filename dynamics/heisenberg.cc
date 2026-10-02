/**
 * @file heisenberg.cc
 * @brief Implementation of heisenberg.h (interfaces documented in the header, logic commented here).
 */
#include "heisenberg.h"
#include <itensor/all.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace itensor;


// conjugate and exchange s <-> s' on every tensor; dag also flips the arrows of quantum numbers
MPO
make_adjoint_mpo(const MPO& U)
{
    MPO Ud = U;
    for(int j = 1 ; j <= length(U) ; j++) Ud.ref(j) = swapPrime(dag(U(j)), 0, 1, "Site");
    return Ud;
}


// ----------------------------------------------------------
// one gate on an MPO

// g^dag T g for the tensor T of the gate sites (site indices s, s'; links untouched):
// O g contracts the input of O with the output of g (both moved to prime level 2), then g^dag
// contracts its input with the output of O g, g^dag(a -> s') = conj(g(s' -> a)).
static ITensor
conjugate_with_gate(const ITensor& T, const ITensor& g)
{
    ITensor Og = mapPrime(T, 0, 2, "Site") * mapPrime(g, 1, 2);
    ITensor gd = conj(g);
    gd.mapPrime(1, 2);
    gd.mapPrime(0, 1);
    return mapPrime(Og, 1, 2, "Site") * gd;
}


// unprimed site index of the MPO tensor at position j
static Index
find_site_index(const MPO& O, const int j)
{
    for(const Index& s : siteInds(O, j)) if(primeLevel(s) == 0) return s;
    throw ITError("apply_gate_to_operator: MPO tensor without site index");
}


// Split T over the positions first, first+1, ...: position first+k receives the site indices
// (s, s') of site_order[k], with one truncated SVD per bond, from left to right.
static void
split_over_sites(MPO& O, ITensor T, const int first, const vector<Index>& site_order, const Args& args)
{
    Index left = (first > 1) ? leftLinkIndex(O, first) : Index();
    for(size_t k = 0 ; k + 1 < site_order.size() ; k++)
    {
        IndexSet left_inds = left ? IndexSet(left, site_order[k], prime(site_order[k]))
                                  : IndexSet(site_order[k], prime(site_order[k]));
        ITensor U, S, V;
        tie(U, S, V) = svd(T, left_inds, args);
        O.set(first + k, U);
        left = commonIndex(U, S);
        T = S * V;
    }
    O.set(first + site_order.size() - 1, T);
}


// contract the tensors of the gate sites, conjugate with the gate, split back (with the first two
// site indices exchanged for a swap_after gate)
MPO
apply_gate_to_operator(MPO O, TebdGate gate, const Args args)
{
    vector<int> sites = gate.sites();
    sort(sites.begin(), sites.end());
    int first = sites[0];
    if(first < 1 || sites.back() > length(O) || sites.back() - first != int(sites.size()) - 1 || sites.size() > 3)
        throw ITError("apply_gate_to_operator: the gate must act on 1 to 3 consecutive sites of the MPO");

    O.position(first);
    ITensor T = O(first);
    vector<Index> site_order = {find_site_index(O, first)};
    for(size_t k = 1 ; k < sites.size() ; k++)
    {
        T *= O(first + k);
        site_order.push_back(find_site_index(O, first + k));
    }
    if(gate.swap_after()) swap(site_order[0], site_order[1]);

    Args svd_args = {"Cutoff=", args.getReal("Cutoff"), "MaxDim=", args.getInt("MaxDim")};
    split_over_sites(O, conjugate_with_gate(T, gate.gate()), first, site_order, svd_args);
    return O;
}


// ----------------------------------------------------------
// one time step

// U = g_n ... g_1, so U^dag O U = g_1^dag ( ... (g_n^dag O g_n) ... ) g_1: the last gate goes first
MPO
heisenberg_step(MPO O, const vector<TebdGate>& gates, const Args args)
{
    for(int k = int(gates.size()) - 1 ; k >= 0 ; k--) O = apply_gate_to_operator(O, gates[k], args);
    return O;
}


// nmultMPO(A, prime(B)) contracts the output of A with the input of B: it is the product B A,
// so O U = nmultMPO(U, prime(O)) and U^dag (O U) = nmultMPO(O U, prime(U^dag))
MPO
heisenberg_step(MPO O, const MPO& U, const Args args)
{
    O = nmultMPO(U, prime(O), args);
    O.mapPrime(2, 1);
    O = nmultMPO(O, prime(make_adjoint_mpo(U)), args);
    O.mapPrime(2, 1);
    return O;
}
