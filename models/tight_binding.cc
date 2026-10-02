/**
 * @file tight_binding.cc
 * @brief Implementation of tight_binding.h (interfaces documented in the header, logic commented here).
 */
#include "tight_binding.h"
#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;


// AutoMPO: hopping of both spins on every bond and the fields on every site, with the coefficient
// vectors read cyclically (index modulo their length)
MPO
make_tight_binding_mpo(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup, const vector<double> hdn)
{
    int size_J = J.size();
    int size_hup = hup.size();
    int size_hdn = hdn.size();


    AutoMPO ampo(sites);

    for(int j : range1(N-1))
    {
        
        ampo += J[(j-1)%size_J] , "Cdagup" , j   , "Cup" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagup" , j+1 , "Cup" , j  ;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j   , "Cdn" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j+1 , "Cdn" , j  ;
        
    }

    // on-site fields on every site
    for(int j : range1(N))
    {
        ampo += hup[(j-1)%size_hup] , "Nup" , j;
        ampo += hdn[(j-1)%size_hdn] , "Ndn" , j;
    }


    return toMPO(ampo);
}


// tridiagonal matrix: h on the diagonal, J_first between sites 1 and 2, J_bulk elsewhere

static ITensor
make_hopping_matrix(const int N, const double h, const double J_first, const double J_bulk, const bool spinful)
{
    int L = spinful ? 2 * N : N;
    Index r = Index(L);
    Index l = prime(r);

    ITensor H = ITensor(r,l);
    for(int i = 1 ; i <= L ; i++)
    {
        H.set(l(i),r(i),h);
        if(i < L)
        {
            double J = (i == 1) ? J_first : J_bulk;
            H.set(l(i),r(i+1),J);
            H.set(l(i+1),r(i),J);
        }
    }
    return H;
}


// homogeneous tridiagonal matrix
ITensor
make_single_particle_hamiltonian(const int N , const vector<double> J, const vector<double> h, const bool spinful)
{
    return make_hopping_matrix(N, h[0], J[0], J[0], spinful);
}


// tridiagonal matrix with a different first hopping
ITensor
make_single_particle_hamiltonian_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful)
{
    return make_hopping_matrix(N, h[0], J[0], J[1], spinful);
}


// Fields divided by the gate counts, plus the hopping written with the operators A, Adag and the
// Jordan-Wigner string F (https://itensor.org/docs.cgi?page=tutorials/fermions); for a mirrored
// chain the order of the two sites in the strings is reversed.

ITensor
make_free_fermion_bond_hamiltonian(const SiteSet sites, const int j, const double J, const vector<double> hup, const vector<double> hdn, const int count_left, const int count_right, const bool mirrored)
{
    int k = j+1;
    ITensor H_S  = (hup[0] * op(sites,"Nup",j) + hdn[0] * op(sites,"Ndn",j)) * op(sites,"Id",k) / double(count_left);
    H_S         += op(sites,"Id",j) * (hup[1] * op(sites,"Nup",k) + hdn[1] * op(sites,"Ndn",k)) / double(count_right);

    ITensor H_SS;
    if(!mirrored)
    {
        H_SS  = J * ( op(sites,"Adagup*F",j) * op(sites,"Aup",k)   - op(sites,"Aup*F",j) * op(sites,"Adagup",k)   );
        H_SS += J * ( op(sites,"Adagdn",j)   * op(sites,"F*Adn",k) - op(sites,"Adn",j)   * op(sites,"F*Adagdn",k) );
    }
    else
    {
        H_SS  = J * ( op(sites,"Aup",j)   * op(sites,"Adagup*F",k) - op(sites,"Adagup",j)   * op(sites,"Aup*F",k) );
        H_SS += J * ( op(sites,"F*Adn",j) * op(sites,"Adagdn",k)   - op(sites,"F*Adagdn",j) * op(sites,"Adn",k)   );
    }
    return H_S + H_SS;
}


// one bond term per bond (fields shared between the neighbouring bonds), exponentiated with
// ITensor BondGate, then the reversed sweep
vector<TebdGate>
make_free_fermion_gates(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{
    int N = length(sites);
    vector<TebdGate> gates;
    for(int j = 1 ; j <= N-1 ; j++)
    {
        ITensor H = make_free_fermion_bond_hamiltonian(sites, j, J[j-1], {hup[j-1], hup[j]}, {hdn[j-1], hdn[j]},
                                                       count_gates_containing(j, 1, 2, N), count_gates_containing(j+1, 1, 2, N));
        gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()));
    }
    return make_symmetric_sweep(gates);
}
