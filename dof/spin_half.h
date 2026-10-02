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
 * @param in        Input (ket) index; dag-ed, so the operator also works with quantum numbers
 *                  (then only "z" has a well-defined flux).
 * @param out       Output index (usually prime(in)).
 * @param direction "x", "y" or "z".
 * @throws ITError for a direction other than "x", "y", "z".
 */
ITensor
make_pauli_operator(const Index& in, const Index& out, const string& direction);

/**
 * @brief Operator measured by measure_magnetization on the site with index s: the Pauli matrix
 *        sigma^direction on a spin-1/2, the number operator on a boson (any direction), zero otherwise.
 * @param s         Site index (the operator has indices s, s').
 * @param direction "x", "y" or "z" (ignored for bosons).
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
 * @brief Configuration string (for make_product_state, make_init_state) of a standard state of N spins.
 * @param N    Number of sites.
 * @param name "up" (00...0), "down" (11...1), "wall" (0...01...1, first N/2 sites '0') or "neel" (0101...).
 * @throws ITError for another name.
 */
string
make_standard_config(const int N, const string& name);

/**
 * @brief InitState of a configuration in the z basis ('0' = "Up", '1' = "Dn"), e.g. to fix the
 *        sector of a search with conserved quantum numbers (find_ground_state).
 * @param sites  Spin-1/2 site set of N sites.
 * @param config One character per site, '0' or '1'.
 * @throws ITError on a wrong length or characters other than '0'/'1'.
 */
InitState
make_init_state(const SiteSet& sites, const string& config);

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
 * @brief As measure_magnetization for a density matrix: Tr(rho sigma^direction_j) on spin-1/2 sites,
 *        Tr(rho n_j) on boson sites.
 * @param rho       Density matrix (MPO with indices s, s'), normalized to Tr(rho) = 1.
 * @param sites     Site set.
 * @param direction "x", "y" or "z".
 */
vector<double>
measure_magnetization(const MPO& rho, const SiteSet& sites, const string& direction);

/**
 * @brief Number of kinks, sum_j <n_j (1 - n_{j+1})> (pairs |down_z>_j |up_z>_{j+1}).
 * @param psi   State; its orthogonality center is moved.
 * @param sites Spin-1/2 site set (N > 1).
 * @return The expected number of kinks.
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
