/**
 * @file tight_binding.h
 * @brief Quadratic fermionic models: MPOs, single-particle matrices and TEBD gates.
 */
#ifndef MYTN_MODELS_TIGHT_BINDING_H
#define MYTN_MODELS_TIGHT_BINDING_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/**
 * @brief MPO of spinful fermions on N sites,
 *        H = sum_{j,sigma} J_j (c^dag_{j,sigma} c_{j+1,sigma} + h.c.) + sum_j (hup_j n_{j,up} + hdn_j n_{j,dn}).
 * @param N     Number of sites.
 * @param sites Electron site set.
 * @param J     Hoppings, used cyclically (a single value gives a homogeneous chain).
 * @param hup   Fields on up electrons, used cyclically (default 0).
 * @param hdn   Fields on down electrons, used cyclically (default 0).
 */
MPO
make_tight_binding_mpo(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup = {0.}, const vector<double> hdn = {0.});

/**
 * @brief Single-particle matrix h of the homogeneous chain, H = sum_{ij} c^dag_i h_ij c_j,
 *        with h_ii = h[0] and h_{i,i+1} = h_{i+1,i} = J[0].
 * @param N       Number of sites.
 * @param J       Hopping J[0] (only the first entry is used).
 * @param h       On-site energy h[0] (only the first entry is used).
 * @param spinful If true the matrix has dimension 2N, otherwise N.
 * @return ITensor with indices (r, r').
 */
ITensor
make_single_particle_hamiltonian(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

/**
 * @brief As make_single_particle_hamiltonian, with an impurity on site 1.
 * @param N       Number of sites.
 * @param J       {hopping between sites 1 and 2, hopping elsewhere}.
 * @param h       On-site energy h[0].
 * @param spinful If true the matrix has dimension 2N, otherwise N.
 */
ITensor
make_single_particle_hamiltonian_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

/**
 * @brief Two-site term on the bond (j, j+1) of free spinful fermions (Electron sites): hopping J for
 *        both spins and the fields of sites j and j+1, divided by the number of gates sharing them
 *        (count_left, count_right; see count_gates_containing in mps/gates.h).
 * @param sites       Electron site set.
 * @param j           Left site of the bond.
 * @param J           Hopping.
 * @param hup         Fields on up electrons {site j, site j+1}.
 * @param hdn         Fields on down electrons {site j, site j+1}.
 * @param count_left  Number of gates sharing the field of site j.
 * @param count_right Number of gates sharing the field of site j+1.
 * @param mirrored    Jordan-Wigner strings for a chain stored in reversed order (bra half of a
 *                    purified state, see models/impurity.h).
 */
ITensor
make_free_fermion_bond_hamiltonian(const SiteSet sites, const int j, const double J, const vector<double> hup, const vector<double> hdn, const int count_left, const int count_right, const bool mirrored = false);

/**
 * @brief Second-order Trotter gates of free spinful fermions (Electron sites):
 *        hoppings J_j between j and j+1 for both spins, fields hup_j n_{j,up} + hdn_j n_{j,dn}.
 * @param sites Electron site set of N sites.
 * @param J     Hoppings (N-1 values).
 * @param hup   Fields on up electrons (N values).
 * @param hdn   Fields on down electrons (N values).
 * @param dt    Time step.
 */
vector<TebdGate>
make_free_fermion_gates(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

#endif
