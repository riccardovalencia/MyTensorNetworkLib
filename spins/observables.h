#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <itensor/all.h>

using namespace itensor;

double
entanglement_entropy( MPS* , int , int );

void 
measure_mx_mz( const SpinHalf , MPS , const int );

#endif
