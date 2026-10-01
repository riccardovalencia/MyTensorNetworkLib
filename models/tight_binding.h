/**
 * @file tight_binding.h
 * @brief Quadratic fermionic models: MPOs, single-particle matrices and TEBD gates.
 */
#ifndef MYTN_MODELS_TIGHT_BINDING_H
#define MYTN_MODELS_TIGHT_BINDING_H

#include <itensor/all.h>

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
hamiltonian_tight_binding_electrons(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup = {0.}, const vector<double> hdn = {0.});

/**
 * @brief Single-particle matrix h of the homogeneous chain, H = sum_{ij} c^dag_i h_ij c_j,
 *        with h_ii = h[0] and h_{i,i+1} = h_{i+1,i} = J[0].
 * @param N       Number of sites.
 * @param spinful If true the matrix has dimension 2N, otherwise N.
 * @return ITensor with indices (r, r').
 */
ITensor
hamiltonian_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

/** @brief As hamiltonian_number_conserving_fermions, with hopping J[0] between sites 1 and 2 (impurity) and J[1] elsewhere. */
ITensor
hamiltonian_number_conserving_fermions_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

/**
 * @brief Second-order Trotter gates of free spinful fermions (Electron sites):
 *        hoppings J_j between j and j+1 for both spins, fields hup_j n_{j,up} + hdn_j n_{j,dn}.
 */
vector<BondGate>
gates_free_spinful_fermions(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt);

#endif
