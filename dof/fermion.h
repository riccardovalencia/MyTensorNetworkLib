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
 *        first Ndnfill sites (doubly occupied where both apply).
 * @param sites   Electron site set.
 * @param Nupfill Number of up electrons (<= N).
 * @param Ndnfill Number of down electrons (<= N).
 */
MPS
initial_computational_electron_state(const SiteSet sites, const int Nupfill, const int Ndnfill);

#endif
