#ifndef MYTN_CORE_STATE_MANIPULATION_H
#define MYTN_CORE_STATE_MANIPULATION_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Move the site index at position j1 to position j2 through a sequence of SVDs
// ("virtual" swap: the physical state is unchanged, only the ordering of the sites in the MPS).
// Used to apply long-range gates (see e.g. PRR 2, 043255 (2020)).
void
swap_gate( MPS *psi, int j1, int j2, double cut_off, int maxDim);

// Density matrix |psi><psi| of a pure state, as an MPO.
MPO
from_MPS_to_MPDO(MPS psi);

// Same as from_MPS_to_MPDO, but fusing the ket and bra link indices into a single one.
MPO
from_MPS_to_MPDO_v2(MPS psi );

// Reduced density matrix of psi on the sites i..j (inclusive), as an ITensor with
// unprimed (ket) and primed (bra) site indices.
ITensor
exctract_reduced_density_matrix(MPS * psi, int i, int j);

#endif
