#ifndef GENERAL_H
#define GENERAL_H

#include <itensor/all.h>
#include "MyClasses.h"
using namespace std;
using namespace itensor;

MPS
apply_gate(MPS psi,  const ITensor gate, const vector<int> jn, const Args args);

#endif
