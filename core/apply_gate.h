#ifndef MYTN_CORE_APPLY_GATE_H
#define MYTN_CORE_APPLY_GATE_H

#include <itensor/all.h>
#include "MyClasses.h"

using namespace std;
using namespace itensor;

// Apply a 1-, 2- or 3-site gate on the consecutive sites jn of psi and return the new MPS.
// The orthogonality center is moved to jn[0]; the result is truncated with the SVD parameters
// in args ("Cutoff", "MaxDim", both required). The state is not normalized.
// Typical use with MyBondGate:  psi = apply_gate(psi, g.gate(), g.jn(), {"Cutoff=",1E-12,"MaxDim=",64});
MPS
apply_gate(MPS psi,  const ITensor gate, const vector<int> jn, const Args args);

#endif
