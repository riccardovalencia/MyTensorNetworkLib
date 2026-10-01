/**
 * @file entanglement.h
 * @brief Entanglement entropy of an MPS.
 */
#ifndef MYTN_MPS_ENTANGLEMENT_H
#define MYTN_MPS_ENTANGLEMENT_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Von Neumann entanglement entropy across the bond (site, site+1), in base 2.
 *
 * @param psi  State; its orthogonality center is moved to site.
 * @param site Left site of the bond, 1 <= site < N.
 * @return S = -sum_k p_k log2 p_k, with p_k the squared Schmidt values (p_k < 1e-12 ignored).
 */
double
entanglement_entropy( MPS* psi, int site );

/**
 * @brief Same as entanglement_entropy(MPS*, int), in natural-log units (nats).
 *
 * @param psi  State; its orthogonality center is moved to site.
 * @param N    Unused (kept for compatibility with older drivers).
 * @param site Left site of the bond.
 * @return S = -sum_k p_k ln p_k.
 */
double
entanglement_entropy( MPS* psi, int N, int site );

#endif
