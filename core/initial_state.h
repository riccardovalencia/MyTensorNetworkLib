#ifndef MYTN_CORE_INITIAL_STATE_H
#define MYTN_CORE_INITIAL_STATE_H

// Model-independent initial states and tools to embed one MPS into another.

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Product state of spin-1/2 sites, one character of config per site (config.size() == N):
//   basis "z" (default): '0' -> |up_z>, '1' -> |down_z>
//   basis "x"          : '0' -> |+x>,   '1' -> |-x>
//   basis "y"          : '0' -> |+y>,   '1' -> |-y>
// e.g. "0000" all up, "1111" all down, "0011" domain wall, "0101" Neel state.
// Throws ITError for a wrong length, characters other than '0'/'1', or non spin-1/2 sites.
MPS
initial_computational_state(const SiteSet sites , const string config , const string basis = "z");

// Overwrite the sites [start, start+length(psi_seed)-1] of *psi with psi_seed.
// inverted : insert psi_seed in reversed site order
// dagger   : conjugate the tensors of psi_seed (used for the bra half of a purified state)
// The site indices of *psi are kept. Use insert_QN_state for MPS with quantum numbers.
void
insert_state(MPS* psi, MPS psi_seed, const int start, bool inverted, bool dagger);

// Same as insert_state, for MPS with conserved quantum numbers (link indices with flux direction).
void
insert_QN_state(MPS* psi, MPS psi_seed, const int start, bool inverted, bool dagger);

// Insert state_to_insert (with site set sites_state_to_insert, L sites) into *psi_t0
// (site set sites, N sites) starting at site start.
void
insert_state(MPS* psi_t0, MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int start, const int L, const int N);

#endif
