#ifndef MYTN_SPIN_BOSON_MYMPO_H
#define MYTN_SPIN_BOSON_MYMPO_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// MPO of the PXP Hamiltonian H = omega sum_j P_{j-1} X_j P_{j+1}, P = (1+Z)/2.
MPO
mpo_pxp(const SiteSet s, const double omega);

#endif
