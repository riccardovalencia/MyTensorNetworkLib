#ifndef MYTN_SPIN_BOSON_INITIAL_STATE_H
#define MYTN_SPIN_BOSON_INITIAL_STATE_H

// Initial states for mixed spin-boson systems (custom_spin_boson site sets).
// Generic states (computational basis, insert_state, ...) are in core/initial_state.h.

#include <itensor/all.h>
#include "../core/initial_state.h"

using namespace std;
using namespace itensor;

// |n_photon> (x) |theta,phi>^(N-1): Fock state on the boson (site 1) and all spins in the same
// spin-coherent state (polar angle theta, azimuthal angle phi on the Bloch sphere).
MPS
initialize_spin_boson_state(const SiteSet sites , const int n_photon , double theta, double phi);

// As above, with site-dependent angles theta[j], phi[j] for the spins.
MPS
initialize_spin_boson_state(const SiteSet sites , const int n_photon , const vector<double> theta, const vector<double> phi);

#endif
