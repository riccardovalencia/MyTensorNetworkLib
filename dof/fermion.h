/**
 * @file fermion.h
 * @brief Spinful fermions (ITensor Electron sites): product states.
 */
#ifndef MYTN_DOF_FERMION_H
#define MYTN_DOF_FERMION_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Product state with up electrons on the first Nupfill sites and down electrons on the
 *        first Ndnfill sites (doubly occupied where both apply), as an InitState.
 * @param sites   Electron site set.
 * @param Nupfill Number of up electrons (<= N).
 * @param Ndnfill Number of down electrons (<= N).
 * @throws ITError if Nupfill or Ndnfill exceeds the number of sites.
 */
InitState
make_electron_init_state(const SiteSet sites, const int Nupfill, const int Ndnfill);

/**
 * @brief MPS of make_electron_init_state (same arguments).
 * @param sites   Electron site set.
 * @param Nupfill Number of up electrons (<= N).
 * @param Ndnfill Number of down electrons (<= N).
 */
MPS
make_electron_product_state(const SiteSet sites, const int Nupfill, const int Ndnfill);

#endif
