/**
 * @file light_matter.cc
 * @brief Implementation of light_matter.h (the functions are documented in the header).
 */
#include "light_matter.h"
#include "../dof/spin_half.h"
#include "../mps/mps_tools.h"
#include <itensor/all.h>
#include <cmath>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// local operators on a site of the spin-boson chain (boson or spin-1/2), indices (s, s')

// spin lowering (lower = true) or raising operator, zero on a boson
static ITensor
make_spin_ladder(const Index& s, bool lower)
{
    ITensor S = ITensor(s, prime(s));
    if(hasTags(s,"Site,S=1/2"))
    {
        if(lower) S.set(s(1), prime(s)(2), 1.);
        else      S.set(s(2), prime(s)(1), 1.);
    }
    return S;
}

// bosonic annihilation (annihilate = true) or creation operator
static ITensor
make_boson_ladder(const Index& s, bool annihilate)
{
    ITensor A = ITensor(s, prime(s));
    for(int d = 1 ; d < dim(s) ; d++)
    {
        if(annihilate) A.set(s(d+1), prime(s)(d), sqrt(d));
        else           A.set(s(d), prime(s)(d+1), sqrt(d));
    }
    return A;
}


// ----------------------------------------------------------
// local terms omega0 n_a + h Z_j on the bond (j, j+1), shared with the neighbouring bonds,
// plus the spin-spin interaction if both sites are spins

static ITensor
make_local_bond_hamiltonian(const SiteSet& sites, const int j, const double omega0, const double h, const double V, const string& interaction_axis)
{
    int N = length(sites);
    Index s1 = sites(j);
    Index s2 = sites(j+1);
    auto field = [&](const Index& s) { return hasTags(s,"Site,Boson") ? omega0 : h; };

    ITensor I1 = make_identity_operator(s1, prime(s1)), I2 = make_identity_operator(s2, prime(s2));
    ITensor D1 = make_magnetization_operator(s1, "z"), D2 = make_magnetization_operator(s2, "z");   // n_a or Z_j

    ITensor H = field(s1) / count_gates_containing(j,   1, 2, N) * D1 * I2
              + field(s2) / count_gates_containing(j+1, 1, 2, N) * I1 * D2;

    if(hasTags(s1,"Site,S=1/2") && hasTags(s2,"Site,S=1/2"))
    {
        if(interaction_axis == "z")      H += V * ((I1 - D1) / 2.) * ((I2 - D2) / 2.);   // V n_j n_{j+1}
        else if(interaction_axis == "x") H += V * make_pauli_operator(s1, prime(s1), "x") * make_pauli_operator(s2, prime(s2), "x");
        else throw ITError("make_light_matter_gates: interaction_axis must be \"z\" or \"x\"");
    }
    return H;
}


// coupling between the boson (site 1) and the spin on site j+1

static ITensor
make_photon_matter_hamiltonian(const SiteSet& sites, const int j, const double g, const string& type_of_coupling)
{
    ITensor A  = make_boson_ladder(sites(1), true);
    ITensor Ad = make_boson_ladder(sites(1), false);
    ITensor Sm = make_spin_ladder(sites(j+1), true);
    ITensor Sp = make_spin_ladder(sites(j+1), false);

    if(type_of_coupling == "dicke") return g * (A + Ad) * (Sp + Sm);
    if(type_of_coupling == "tavis") return g * (A * Sp + Ad * Sm);
    throw ITError("make_light_matter_gates: type_of_coupling must be \"dicke\" or \"tavis\"");
}


// ----------------------------------------------------------
// short-range gates: local terms on (j, j+1); long-range gates: boson-spin coupling on (j, j+1),
// with the boson on position j, followed by a swap that moves the boson forward
// (swap gates with ITensor BondGate do not work with different local dimensions)

vector<TebdGate>
make_light_matter_gates(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling, const double V, string interaction_axis)
{
    int N = length(sites);
    vector<TebdGate> gates;
    for(int j = 1 ; j < N ; j++)
    {
        if(matter_or_photon == "short-range")
        {
            ITensor hj = make_local_bond_hamiltonian(sites, j, omega0, h, V, interaction_axis);
            gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj).gate()));
        }
        else if(matter_or_photon == "long-range")
        {
            ITensor hj = make_photon_matter_hamiltonian(sites, j, g, type_of_coupling);
            gates.push_back(TebdGate({j,j+1}, BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj).gate(), true));
        }
        else throw ITError("make_light_matter_gates: matter_or_photon must be \"short-range\" or \"long-range\"");
    }
    return make_symmetric_sweep(gates);
}
