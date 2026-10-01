/**
 * @file spin_chain.cc
 * @brief Implementation of spin_chain.h (interfaces documented in the header, logic commented here).
 */
#include "spin_chain.h"
#include "../mps/gates.h"
#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// MPO

// Jxx X_j X_k + Jyy Y_j Y_k + Jzz Z_j Z_k, with X = S+ + S-, Y = -i (S+ - S-), Z = 2 Sz;
// terms with zero coefficient are left out (they would break the S^z conservation)

static void
add_spin_coupling(AutoMPO& ampo, const vector<double> J, const int j, const int k)
{
    double same = J[0] - J[1], opposite = J[0] + J[1];
    if(same != 0.)
    {
        ampo += same, "S+", j, "S+", k;
        ampo += same, "S-", j, "S-", k;
    }
    if(opposite != 0.)
    {
        ampo += opposite, "S+", j, "S-", k;
        ampo += opposite, "S-", j, "S+", k;
    }
    if(J[2] != 0.) ampo += 4 * J[2], "Sz", j, "Sz", k;
}


// AutoMPO with fields on every site, J on the bonds (j, j+1) and J2 on (j, j+2)
MPO
make_spin_chain_mpo(const SiteSet sites, const vector<double> J, const vector<double> J2, const vector<double> h)
{
    int N = length(sites);
    auto ampo = AutoMPO(sites);
    for(int j = 1 ; j <= N ; j++)
    {
        if(h[0] != 0.) { ampo += h[0], "S+", j;  ampo += h[0], "S-", j; }
        if(h[1] != 0.) { ampo += -Cplx_i * h[1], "S+", j;  ampo += Cplx_i * h[1], "S-", j; }
        if(h[2] != 0.) ampo += 2 * h[2], "Sz", j;
    }
    for(int j = 1 ; j < N   ; j++) add_spin_coupling(ampo, J,  j, j+1);
    for(int j = 1 ; j < N-1 ; j++) add_spin_coupling(ampo, J2, j, j+2);
    return toMPO(ampo);
}


// ----------------------------------------------------------
// fields of site j divided among the count bonds sharing it, plus the couplings of the bond (j, j+1)

ITensor
make_spin_chain_bond_hamiltonian(const SiteSet sites, const vector<double> J, const vector<double> h, const int j, const int count_left, const int count_right)
{
    vector<ITensor> X, Y, Z, Id;
    for(int q = j ; q <= j+1 ; q++)
    {
        Id.push_back(    op(sites,"Id",q) );
        X.push_back( 2 * op(sites,"Sx",q) );
        Y.push_back( 2 * op(sites,"Sy",q) );
        Z.push_back( 2 * op(sites,"Sz",q) );
    }

    ITensor H_S  = (h[0] * X[0] + h[1] * Y[0] + h[2] * Z[0]) * Id[1] / double(count_left);
    H_S         += Id[0] * (h[0] * X[1] + h[1] * Y[1] + h[2] * Z[1]) / double(count_right);
    ITensor H_SS = J[0] * X[0] * X[1] + J[1] * Y[0] * Y[1] + J[2] * Z[0] * Z[1];
    return H_S + H_SS;
}


// L^dag L of a local jump operator L (indices s, s')

static ITensor
make_jump_norm_operator(const ITensor& L)
{
    ITensor Ld = conj(L);
    Ld.mapPrime(0,2);           // (L^*)^T L = L^dag L: the 'row' index of L^* becomes the 'column' one
    ITensor LdL = Ld * L;
    LdL.mapPrime(2,1);
    return LdL;
}


// bond term of an open chain: the fields of the edge sites are not shared

static ITensor
make_open_chain_bond_hamiltonian(const SiteSet& sites, const vector<double>& J, const vector<double>& h, int j)
{
    int N = length(sites);
    return make_spin_chain_bond_hamiltonian(sites, J, h, j, count_gates_containing(j, 1, 2, N), count_gates_containing(j+1, 1, 2, N));
}


// one gate exp(-i dt/2 h_j) per bond, then the reversed sweep
vector<TebdGate>
make_spin_chain_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt)
{
    vector<TebdGate> gates;
    for(int j = 1 ; j <= length(sites)-1 ; j++)
        gates.push_back(TebdGate({j,j+1}, dt/2., make_open_chain_bond_hamiltonian(sites, J, h, j)));
    return make_symmetric_sweep(gates);
}


// the non-hermitian term -i/2 gamma L^dag L of a jump operator on site j enters the bond (j, j+1),
// and the bond (N-1, N) for j = N

vector<TebdGate>
make_spin_chain_effective_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const vector<int> Lj_sites, const vector<double> gamma, const double dt)
{
    int N = length(sites);
    vector<TebdGate> gates;
    for(int j = 1 ; j <= N-1 ; j++)
    {
        ITensor H = make_open_chain_bond_hamiltonian(sites, J, h, j);
        for(int k = 0 ; k < (int)Lj_sites.size() ; k++)
        {
            ITensor LdL = -0.5 * gamma[k] * Cplx_i * make_jump_norm_operator(Lj[k]);
            if(Lj_sites[k] == j)                H += LdL * op(sites,"Id",j+1);
            if(Lj_sites[k] == N && j == N-1)    H += op(sites,"Id",j) * LdL;
        }
        gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()));
    }
    return make_symmetric_sweep(gates);
}


// one single-site gate exp(-i dt/2 w.sigma) per site, then the reversed sweep
vector<TebdGate>
make_local_field_gates(const SiteSet sites , vector<double> omegaj, const double dt)
{
    vector<TebdGate> gates;
    for(int j = 1 ; j <= length(sites) ; j++)
    {
        ITensor hj = omegaj[0] * 2 * op(sites,"Sx",j) + omegaj[1] * 2 * op(sites,"Sy",j) + omegaj[2] * 2 * op(sites,"Sz",j);
        gates.push_back(TebdGate({j}, dt/2., hj));
    }
    return make_symmetric_sweep(gates);
}


// H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ] is the spin chain with J = {-J, 0, 0}, h = {-J hx, 0, -J hz}

vector<TebdGate>
make_ising_gates( const SiteSet sites , const double J , const double hx , const double hz , const double dt )
{
    return make_spin_chain_gates(sites, {-J, 0., 0.}, {-J * hx, 0., -J * hz}, dt);
}
