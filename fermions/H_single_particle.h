#ifndef MYTN_FERMIONS_H_SINGLE_PARTICLE_H
#define MYTN_FERMIONS_H_SINGLE_PARTICLE_H

// Single-particle Hamiltonians h of quadratic fermionic models, H = sum_{ij} c^dag_i h_ij c_j,
// returned as an L x L ITensor with indices (r, r') - L = N, or 2N if spinful.

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Homogeneous tight-binding chain: h_ii = h[0], h_{i,i+1} = h_{i+1,i} = J[0].
ITensor
H_number_conserving_fermions(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

// Impurity model: hopping J[0] between site 1 (the impurity) and site 2, J[1] elsewhere; h_ii = h[0].
ITensor
H_number_conserving_fermions_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful = false);

#endif
