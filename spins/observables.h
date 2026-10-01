#ifndef MYTN_SPINS_OBSERVABLES_H
#define MYTN_SPINS_OBSERVABLES_H

#include <itensor/all.h>
#include "../core/entanglement.h"

using namespace std;
using namespace itensor;

// Print <X_j> and <Z_j> for j = 1..N to stdout.
void
measure_mx_mz( const SpinHalf sites , MPS psi , const int N );

#endif
