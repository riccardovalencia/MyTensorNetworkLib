#ifndef MYTN_SPIN_BOSON_PARAMETERS_H
#define MYTN_SPIN_BOSON_PARAMETERS_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Nearest-neighbour couplings V_j = 1/|r_j - r_{j+1}|^alpha for positions rj = {{x,y,z}, ...}.
// Returns N-1 values (e.g. alpha = 6 for van der Waals interactions between Rydberg atoms).
vector<double>
compute_potential(const vector< vector<double> > rj , const double alpha);

#endif
