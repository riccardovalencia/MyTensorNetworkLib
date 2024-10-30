#ifndef SCALAR_PRODUCT_DIFFERENT_CUTOFFS
#define SCALAR_PRODUCT_DIFFERENT_CUTOFFS
//get_data.h

#include <itensor/all.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;
using namespace itensor;

tuple<double, double>
scalar_product_different_cutoff( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c);
#endif
