#ifndef MYTN_CORE_ENTANGLEMENT_H
#define MYTN_CORE_ENTANGLEMENT_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Von Neumann entanglement entropy (base 2) across the bond [site, site+1].
// Moves the orthogonality center of psi to site.
double
entanglement_entropy( MPS* psi, int site );

// Same, but in natural-log units. N is unused (kept for backward compatibility of the spins module).
double
entanglement_entropy( MPS* psi, int N, int site );

#endif
