#ifndef SCALAR_PRODUCT_DIFFERENT_N0
#define SCALAR_PRODUCT_DIFFERENT_N0
//get_data.h

#include <itensor/all.h>
#include <vector>
#include <iomanip>
#include <complex>

using namespace std;
using namespace itensor;

tuple<double, double>
scalar_product_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c);
#endif
