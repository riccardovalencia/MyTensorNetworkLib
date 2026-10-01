#ifndef MYTN_BOSONS_SCALAR_PRODUCT_TWO_CUTOFFS_H
#define MYTN_BOSONS_SCALAR_PRODUCT_TWO_CUTOFFS_H
//get_data.h

#include <itensor/all.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;
using namespace itensor;

// Overlap |<psi1|psi2>| between states with different Fock-space cutoffs (psi1 is embedded in the larger space),
// and variance of the bQEM Hamiltonian (parameters n0, symmetry_sector, s, c) on the embedded state.
// Returns {overlap, variance}.
tuple<double, double>
scalar_product_different_cutoff( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c);
#endif
