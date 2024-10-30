#ifndef H_FERMIONIC
#define H_FERMIONIC

#include <itensor/all.h>

using namespace std;
using namespace itensor;

MPO
H_tight_binding_electrons(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup = {0.}, const vector<double> hdn = {0.});

#endif
