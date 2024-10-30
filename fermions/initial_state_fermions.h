#ifndef INITIAL_STATE_FERMIONS
#define INITIAL_STATE_FERMIONS

#include <itensor/all.h>

using namespace std;
using namespace itensor;

MPS 
fermi_sea_electrons(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, const Sweeps, double min_varH);

MPS
initial_computational_electron_state(const SiteSet sites, const int Nupfill, const int Ndnfill);

#endif
