/**
 * @file spin_half.h
 * @brief Spin-1/2 degrees of freedom: product states and local observables.
 *
 * Conventions: |0> = |up_z>, |1> = |down_z>; X, Y, Z are Pauli matrices;
 * n = (1 - Z)/2 = |down_z><down_z| is the excitation (Rydberg) projector.
 */
#ifndef MYTN_DOF_SPIN_HALF_H
#define MYTN_DOF_SPIN_HALF_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Pauli matrix sigma^direction between the indices in (input) and out (output) of a spin-1/2,
 *        with basis state 1 = |up_z>, 2 = |down_z>; e.g. make_pauli_operator(s, prime(s), "x").
 * @throws ITError for a direction other than "x", "y", "z".
 */
ITensor
make_pauli_operator(const Index& in, const Index& out, const string& direction);

/**
 * @brief Operator measured by measure_magnetization on the site with index s: the Pauli matrix
 *        sigma^direction on a spin-1/2, the number operator on a boson (any direction), zero otherwise.
 */
ITensor
make_magnetization_operator(const Index& s, const string& direction);

/**
 * @brief Product state of spin-1/2 sites from a string of 0 and 1.
 *
 * | basis | '0'   | '1'   |
 * |-------|-------|-------|
 * | "z"   | up_z  | down_z|
 * | "x"   | +x    | -x    |
 * | "y"   | +y    | -y    |
 *
 * @code
 * MPS kink  = make_product_state(sites, "1100000000");      // z basis
 * MPS wall  = make_product_state(sites, "0000011111", "x"); // domain wall along x
 * @endcode
 *
 * @param sites  Spin-1/2 site set of N sites.
 * @param config One character per site, '0' or '1' (config.size() == N).
 * @param basis  "z" (default), "x" or "y".
 * @return Normalized product state with orthogonality center on site 1.
 * @throws ITError on a wrong length, characters other than '0'/'1', an unknown basis
 *         or sites that are not spin-1/2.
 */
MPS
make_product_state(const SiteSet sites , const string config , const string basis = "z");

/**
 * @brief Local magnetization (spins) or occupation (bosons) on every site.
 * @param psi       State; its orthogonality center is moved.
 * @param sites     Site set (spin-1/2 and/or boson sites).
 * @param direction "x", "y" or "z".
 * @return For each site j: <sigma^direction_j> on spin-1/2 sites, <n_j> on boson sites.
 */
vector<double>
measure_magnetization(MPS* psi, const SiteSet sites , string direction);

/**
 * @brief Number of kinks, sum_j <n_j (1 - n_{j+1})> (pairs |down_z>_j |up_z>_{j+1}).
 * @param psi   State; its orthogonality center is moved.
 * @param sites Spin-1/2 site set (N > 1).
 */
double
measure_kink_number( MPS* psi, const SiteSet sites);

/**
 * @brief Density-density correlations <n_start n_j> for j = 1..N.
 * @param psi       State; its orthogonality center is moved.
 * @param sites     Spin-1/2 site set.
 * @param start     Reference site.
 * @param connected If true, return <n_start n_j> - <n_start><n_j>.
 * @return N values, one per site j.
 */
vector<double>
measure_density_correlations(MPS* psi, const SiteSet sites, const int start, const bool connected);

/**
 * @brief Print <X_j> and <Z_j> for j = 1..N to stdout.
 * @param sites Spin-1/2 site set.
 * @param psi   State (copied).
 * @param N     Number of sites.
 */
void
print_magnetization( const SpinHalf sites , MPS psi , const int N );

#endif
