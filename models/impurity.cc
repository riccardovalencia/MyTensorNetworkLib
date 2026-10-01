/**
 * @file impurity.cc
 * @brief Implementation of impurity.h (the functions are documented in the header).
 */
#include "impurity.h"
#include "spin_chain.h"
#include "tight_binding.h"
#include "../mps/gates.h"
#include <itensor/all.h>
#include <algorithm>
#include <vector>

using namespace std;
using namespace itensor;


// position (1..N) of site j of the doubled chain inside its half: bra 1..N, ket N+1..2N

static int
position_in_half(const int j, const int N)
{
    return (j > N) ? j - N : j;
}


// number of two-site gates of the half sharing the on-site terms of site j of the doubled chain

static int
count_bonds_in_half(const int j, const int N)
{
    return count_gates_containing(position_in_half(j, N), 1, 2, N);
}


// the bra evolves with exp(+i H^T t), i.e. with the generator -H^T

static ITensor
transpose_for_bra(const ITensor& H)
{
    return swapPrime(-H, 0, 1);
}


vector<TebdGate>
make_kondo_impurity_gates(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{
    int N = length(sites)/2;

    // couplings on the doubled chain: bra (reversed), central bra-ket bond (no hopping), ket
    vector<double> J2(J.rbegin(), J.rend());
    J2.push_back(0.);
    J2.insert(J2.end(), J.begin(), J.end());

    vector<double> hup2(hup.rbegin(), hup.rend());
    hup2.insert(hup2.end(), hup.begin(), hup.end());
    vector<double> hdn2(hdn.rbegin(), hdn.rend());
    hdn2.insert(hdn2.end(), hdn.begin(), hdn.end());

    vector<TebdGate> gates;
    for(int j = 1 ; j <= 2*N-1 ; j++)
    {
        if(j == N) continue;
        bool bra = (j < N);
        ITensor H = make_free_fermion_bond_hamiltonian(sites, j, J2[j-1], {hup2[j-1], hup2[j]}, {hdn2[j-1], hdn2[j]},
                                                       count_bonds_in_half(j, N), count_bonds_in_half(j+1, N), bra);
        if(bra) H = transpose_for_bra(H);
        gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()));
    }
    return make_symmetric_sweep(gates);
}


vector<TebdGate>
make_spin_impurity_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt)
{
    int N = length(sites)/2;

    vector<TebdGate> gates;
    for(int j = 1 ; j <= 2*N-1 ; j++)
    {
        if(j == N) continue;
        ITensor H = make_spin_chain_bond_hamiltonian(sites, J, h, j, count_bonds_in_half(j, N), count_bonds_in_half(j+1, N));
        if(j < N) H = transpose_for_bra(H);
        gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()));
    }
    return make_symmetric_sweep(gates);
}


// Three-site term on (j, j+1, j+2), p = position of j in its half: the next-nearest-neighbour
// couplings in full, the fields and the nearest-neighbour couplings divided by the number of
// three-site gates sharing them

static ITensor
make_three_site_spin_hamiltonian(const SiteSet& sites, const int j, const int p, const int N, const vector<double>& J, const vector<double>& J_NNN, const vector<double>& h)
{
    vector<ITensor> I;
    vector<vector<ITensor> > pauli(3);   // pauli[c][a]: X, Y or Z (c = 0, 1, 2) on site j + a
    for(int q = j ; q <= j+2 ; q++)
    {
        I.push_back(op(sites,"Id",q));
        pauli[0].push_back(2 * op(sites,"Sx",q));
        pauli[1].push_back(2 * op(sites,"Sy",q));
        pauli[2].push_back(2 * op(sites,"Sz",q));
    }
    // product of the operators O[a] on the sites a in `on` and of identities elsewhere
    auto embed = [&](const vector<ITensor>& O, const vector<int>& on)
    {
        auto factor = [&](int a) { return (find(on.begin(), on.end(), a) != on.end()) ? O[a] : I[a]; };
        return factor(0) * factor(1) * factor(2);
    };

    ITensor H;
    for(int c = 0 ; c < 3 ; c++)
        for(int a = 0 ; a < 3 ; a++)
            H += h[c] / double(count_gates_containing(p+a, 1, 3, N)) * embed(pauli[c], {a});
    for(int c = 0 ; c < 3 ; c++)
    {
        H += J[c] / double(count_gates_containing(p,   2, 3, N)) * embed(pauli[c], {0, 1});
        H += J[c] / double(count_gates_containing(p+1, 2, 3, N)) * embed(pauli[c], {1, 2});
    }
    for(int c = 0 ; c < 3 ; c++) H += J_NNN[c] * embed(pauli[c], {0, 2});
    return H;
}


vector<TebdGate>
make_spin_impurity_nnn_gates(const SiteSet sites , const vector<double> J, const vector<double> J_NNN, const vector<double> h, const double dt)
{
    int N = length(sites)/2;
    if(N < 3) throw ITError("make_spin_impurity_nnn_gates: three-site gates need at least 3 physical sites");

    // three layers of non-overlapping gates in each half: [1,2,3] [4,5,6] ..., [2,3,4] ..., [3,4,5] ...
    vector<TebdGate> gates_bra, gates_ket;
    for(int layer = 1 ; layer <= 3 ; layer++)
    {
        for(int j = layer ; j <= N-2 ; j += 3)
            gates_bra.push_back(TebdGate({j,j+1,j+2}, dt/2., transpose_for_bra(make_three_site_spin_hamiltonian(sites, j, j, N, J, J_NNN, h))));
        for(int j = N + layer ; j <= 2*N-2 ; j += 3)
            gates_ket.push_back(TebdGate({j,j+1,j+2}, dt/2., make_three_site_spin_hamiltonian(sites, j, j-N, N, J, J_NNN, h)));
    }

    vector<TebdGate> gates = gates_bra;
    gates.insert(gates.end(), gates_ket.begin(), gates_ket.end());
    return make_symmetric_sweep(gates);
}
