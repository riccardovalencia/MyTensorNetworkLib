/**
 * @file entanglement.h
 * @brief Entanglement entropy of an MPS.
 */
#ifndef MYTN_MPS_ENTANGLEMENT_H
#define MYTN_MPS_ENTANGLEMENT_H

#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;

/**
 * @brief Von Neumann entanglement entropy across the bond (site, site+1).
 *
 * @param psi         State; its orthogonality center is moved to site.
 * @param site        Left site of the bond, 1 <= site < N.
 * @param natural_log Use the natural logarithm (entropy in nats) instead of log2 (default).
 * @return S = -sum_k p_k log p_k, with p_k the squared Schmidt values (p_k < 1e-12 ignored).
 */
double
compute_entanglement_entropy( MPS* psi, int site, bool natural_log = false );

/**
 * @brief Entanglement entropy of every cut: S_j across the bond (j, j+1) for j = 1..N-1.
 * @param psi         State; its orthogonality center is moved.
 * @param natural_log Use the natural logarithm instead of log2 (default).
 * @return N-1 values.
 */
vector<double>
compute_entanglement_entropies( MPS* psi, bool natural_log = false );

#endif
