/**
 * @file lindblad.cc
 * @brief Implementation of lindblad.h (the functions are documented in the header).
 */
#include "lindblad.h"
#include "../mps/gates.h"
#include "../mps/mps_tools.h"
#include <itensor/all.h>
#include <cmath>
#include <vector>

using namespace std;
using namespace itensor;

// ----------------------------------------------------------
// Conventions. On each site the ket has index s and the bra s'; a superoperator maps
// (s, s') (input) to (s'', s''') (output):
//              _
//      s'' - |   | - s
//            |   |
//     s''' - |_ _| - s'
//
// so that D[L] = L (x) L^* - 1/2 (L^dag L (x) I) - 1/2 (I (x) (L^dag L)^T).

// gamma D[L] for a jump operator L acting on the sites with indices site_inds (one or more sites)
static ITensor
make_dissipator(const ITensor& L, const vector<Index>& site_inds, const double gamma)
{
    // identities on the ket (s -> s'') and on the bra (s' -> s''')
    ITensor Idket, Idbra;
    for(const Index& s : site_inds)
    {
        Idket = Idket ? Idket * make_identity_operator(s, prime(s,2))        : make_identity_operator(s, prime(s,2));
        Idbra = Idbra ? Idbra * make_identity_operator(prime(s), prime(s,3)) : make_identity_operator(prime(s), prime(s,3));
    }

    ITensor Ld = conj(L);

    // non-hermitian part on the ket: L^dag L = (L^*)^T L (the 'row' index of L^* becomes the 'column' one)
    ITensor LdL_I = mapPrime(Ld, 0, 2) * L * Idbra;

    // non-hermitian part on the bra: (L^dag L)^T, (0,1) -> (2,3)
    ITensor I_LdL = mapPrime(Ld, 0, 3) * L;    // indices (3, 0)
    I_LdL.mapPrime(3,1);                        // (3,0) -> (1,0)
    I_LdL.mapPrime(0,3);                        // (1,0) -> (1,3): transposition
    I_LdL *= Idket;

    // jumps L rho L^dag
    ITensor L_Ld = mapPrime(L, 1, 2) * mapPrime(mapPrime(Ld, 1, 3), 0, 1);

    return gamma * (L_Ld - 0.5 * LdL_I - 0.5 * I_LdL);
}


vector<DissipativeGate>
make_local_dissipative_gates(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt)
{
    vector<DissipativeGate> gates;
    for(int j : lj_sites)
        gates.push_back(DissipativeGate({j}, {j}, dt, make_dissipator(Lj[j-1], {sites(j)}, gammaj[j-1])));
    return gates;
}


vector<DissipativeGate>
make_local_dissipative_gates(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
{
    vector<int> all_sites;
    for(int j = 1 ; j <= length(sites) ; j++) all_sites.push_back(j);
    return make_local_dissipative_gates(sites, Lj, all_sites, gammaj, dt);
}


// Jumps of an OperatorPair (Li on site i, Lj on site j != i, |i - j| = 1):
//   Li rho Lj^dag + Lj rho Li^dag - 1/2 ( {Lj^dag Li, rho} + {Li^dag Lj, rho} )

static ITensor
make_cross_dissipator(const SiteSet& sites, OperatorPair T)
{
    ITensor li  = T.op_i();
    ITensor lj  = T.op_j();
    ITensor lid = dag(li);
    ITensor ljd = dag(lj);

    Index si = sites(T.site_i());
    Index sj = sites(T.site_j());

    // identities on (ket i, ket j), (bra i, bra j), (ket i, bra j), (bra i, ket j)
    ITensor Idket        = make_identity_operator(si, prime(si,2))        * make_identity_operator(sj, prime(sj,2));
    ITensor Idbra        = make_identity_operator(prime(si), prime(si,3)) * make_identity_operator(prime(sj), prime(sj,3));
    ITensor Id_iket_jbra = make_identity_operator(si, prime(si,2))        * make_identity_operator(prime(sj), prime(sj,3));
    ITensor Id_ibra_jket = make_identity_operator(sj, prime(sj,2))        * make_identity_operator(prime(si), prime(si,3));

    // non-hermitian part on the ket
    ITensor LdL_I = (mapPrime(mapPrime(ljd,0,2),1,0) * mapPrime(li,1,2) + mapPrime(mapPrime(lid,0,2),1,0) * mapPrime(lj,1,2)) * Idbra;

    // non-hermitian part on the bra (transposed)
    ITensor I_LdL = mapPrime(ljd,0,3) * mapPrime(mapPrime(li,1,3),0,1) + mapPrime(lid,0,3) * mapPrime(mapPrime(lj,1,3),0,1);
    I_LdL = swapPrime(I_LdL,1,3) * Idket;

    // jumps
    ITensor L_Ld = mapPrime(li,1,2) * mapPrime(mapPrime(ljd,1,3),0,1) * Id_ibra_jket
                 + mapPrime(lj,1,2) * mapPrime(mapPrime(lid,1,3),0,1) * Id_iket_jbra;

    return T.rate() * (L_Ld - 0.5 * LdL_I - 0.5 * I_LdL);
}


vector<DissipativeGate>
make_two_site_dissipative_gates(const SiteSet sites , vector<OperatorPair> TTrain, const double dt)
{
    vector<DissipativeGate> gates;
    for(OperatorPair T : TTrain)
    {
        int i = T.site_i();
        int j = T.site_j();
        if(abs(i-j) > 1) throw ITError("make_two_site_dissipative_gates: only on-site and nearest-neighbour jump operators are implemented");

        if(i == j) gates.push_back(DissipativeGate({j}, {j}, dt/2., make_dissipator(T.op_j(), {sites(j)}, T.rate())));
        else       gates.push_back(DissipativeGate({i,j}, {i,j}, dt/2., make_cross_dissipator(sites, T)));
    }
    return make_symmetric_sweep(gates);
}


vector<DissipativeGate>
make_multisite_dissipative_gates(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt)
{
    if(Lj_sites.size() != Lij_list.size())
        throw ITError("make_multisite_dissipative_gates: Lij_list and Lj_sites have different lengths");

    vector<DissipativeGate> gates;
    for(size_t k = 0 ; k < Lj_sites.size() ; k++)
    {
        vector<int> jn = Lj_sites[k];
        if(jn.size() != 2) throw ITError("make_multisite_dissipative_gates: only two-site jump operators are implemented");
        ITensor D = make_dissipator(Lij_list[k], {sites(jn[0]), sites(jn[1])}, gammaj[k]);
        gates.push_back(DissipativeGate(jn, jn, dt/2., D));
    }
    return make_symmetric_sweep(gates);
}


// ----------------------------------------------------------
// Dissipator on the central bond (N, N+1) of the unfolded density matrix of an impurity problem
// (sites 1..N: bra, evolving with -H; sites N+1..2N: ket, evolving with +H). Lj[0] acts on the bra
// (site N), Lj[1] on the ket (site N+1).

static ITensor
make_impurity_dissipator(const SiteSet& sites, const vector<ITensor>& Lj, const double gamma)
{
    int N = length(sites)/2;
    ITensor lj1 = Lj[0];
    ITensor lj2 = Lj[1];

    // (L^dag L)^T = L^T L^* on the bra, L^dag L on the ket
    ITensor ljdlj1 = conj(lj1) * mapPrime(lj1, 0, 2);
    ITensor ljdlj2 = mapPrime(conj(lj2), 0, 2) * lj2;
    ljdlj1.mapPrime(2,1);
    ljdlj2.mapPrime(2,1);

    // L^dag L (x) I: the first operator acts on the ket (second half of the chain)
    ITensor LdL_I = op(sites,"Id",N) * ljdlj2;
    ITensor I_LdL = ljdlj1 * op(sites,"Id",N+1);

    return gamma * (conj(lj1) * lj2 - 0.5 * LdL_I - 0.5 * I_LdL);
}


vector<DissipativeGate>
make_impurity_dissipative_gates(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt)
{
    int N = length(sites)/2;
    return {DissipativeGate({N,N+1}, {N,N+1}, dt, make_impurity_dissipator(sites, Lj, gamma))};
}


// BondGate exponentiates also non-hermitian generators (Taylor series), beyond the first order in dt
// of make_impurity_dissipative_gates

vector<TebdGate>
make_impurity_dissipative_gates_pade(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt)
{
    int N = length(sites)/2;
    ITensor D = make_impurity_dissipator(sites, Lj, gamma);
    return {TebdGate({N,N+1}, BondGate(sites,N,N+1,BondGate::tImag,-1*dt,D).gate())};
}


// ----------------------------------------------------------

vector<TebdGate>
make_purified_gates(const vector<TebdGate> gates_single, const SiteSet sites_single,  const SiteSet sites_doubled)
{
    // physical site i (1..N) -> ket site N + i, mirrored bra site N + 1 - i
    int N = length(sites_single);
    auto physical_site = [&](const Index& s)
    {
        for(int i = 1 ; i <= N ; i++) if(sites_single(i) == s) return i;
        throw ITError("make_purified_gates: gate index not in sites_single");
    };

    vector<TebdGate> gates_doubled;
    for(bool ket : {true, false})
    {
        for(TebdGate g : gates_single)
        {
            auto map_site = [&](int i) { return ket ? N + i : N + 1 - i; };

            // move every site index of the gate onto the doubled chain
            ITensor original = g.gate();
            ITensor gate     = original;
            for(Index s : inds(original))
            {
                if(primeLevel(s) != 0) continue;
                Index s_new = sites_doubled(map_site(physical_site(s)));
                gate *= delta(s, s_new);
                gate *= delta(prime(s), prime(s_new));
            }

            vector<int> positions;
            for(int i : g.sites()) positions.push_back(map_site(i));

            // the bra evolves with the conjugate gate
            gates_doubled.push_back(TebdGate(positions, ket ? gate : dag(gate), g.swap_after()));
        }
    }
    return gates_doubled;
}


// ----------------------------------------------------------
// rho -> rho + gate * rho on the sites (j, j+1), j = first ket site of the gate

MPS
apply_dissipative_gate(MPS psi, DissipativeGate gate, const Args args)
{
    double cut_off = args.getReal("Cutoff");
    int maxDim     = args.getInt("MaxDim");
    int j          = gate.ket_sites()[0];

    ITensor AA   = psi(j) * psi(j+1);
    ITensor dpsi = gate.gate() * AA;
    dpsi.mapPrime(1,0);
    AA = AA + dpsi;

    auto [U,S,V] = svd(AA,inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
    psi.set(j,U);
    psi.set(j+1,S*V);
    return psi;
}
