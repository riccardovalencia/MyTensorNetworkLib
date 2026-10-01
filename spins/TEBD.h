#ifndef MYTN_SPINS_TEBD_H
#define MYTN_SPINS_TEBD_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Bond term on sites (b, b+1) of the Ising chain in longitudinal (hx) and transverse (hz) fields,
//   H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ],
// with the single-site terms split between neighbouring bonds (full weight on the edges).
// The result is written to *hterm (use it to build an ITensor BondGate).
void
build_single_step( ITensor *hterm , const SpinHalf sites , const int N , const double J , const double hx , const double hz , const int b );

#endif
