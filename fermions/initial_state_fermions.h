#ifndef MYTN_FERMIONS_INITIAL_STATE_FERMIONS_H
#define MYTN_FERMIONS_INITIAL_STATE_FERMIONS_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Ground state of H (e.g. H_tight_binding_electrons) with Nupfill up and Ndnfill down electrons, via DMRG
// starting from initial_computational_electron_state. Exits if the energy variance exceeds min_varH
// or the filling is not reproduced.
MPS
fermi_sea_electrons(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, const Sweeps sweeps, double min_varH);

// Product state with the first Nupfill sites occupied by an up electron and the first Ndnfill by a down one.
MPS
initial_computational_electron_state(const SiteSet sites, const int Nupfill, const int Ndnfill);

#endif
