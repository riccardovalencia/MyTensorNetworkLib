#ifndef STATE_MANIPULATION
#define STATE_MANIPULATION

#include <itensor/all.h>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
// compute reduced density matrix of a given MPS bewtween sites i and j included
ITensor
exctract_reduced_density_matrix(MPS * psi, int i, int j);

#endif
