#ifndef MYTN_SPIN_BOSON_CUSTOM_SITESET_H
#define MYTN_SPIN_BOSON_CUSTOM_SITESET_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// N sites: site 1 is a boson truncated at max_occ, sites 2..N are spin-1/2.
SiteSet
custom_spin_boson(const int N , const int max_occ);

// Doubled (bra-ket) site set of 2N sites for the purified density matrix of custom_spin_boson(N,max_occ):
//   s_{N-1} ... s_1  b  |  b  s_1 ... s_{N-1}
//   ------ bra -------    ------ ket -------
// The bra is mirrored, so that the bosons (where dissipation acts) sit on the central bond.
SiteSet
custom_spin_boson_doubling(const int N , const int max_occ);

#endif
