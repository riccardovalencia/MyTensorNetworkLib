#ifndef MYTN_FERMIONS_H_FERMIONIC_H
#define MYTN_FERMIONS_H_FERMIONIC_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// MPO of spinful fermions (Electron sites) on N sites:
//   H = sum_{j,sigma} J_j (c^dag_{j,sigma} c_{j+1,sigma} + h.c.) + sum_j (hup_j n_{j,up} + hdn_j n_{j,dn})
// J, hup, hdn are used cyclically, so a single value gives a homogeneous chain.
// Note: the on-site fields are added for sites 1..N-1 only.
MPO
H_tight_binding_electrons(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup = {0.}, const vector<double> hdn = {0.});

#endif
