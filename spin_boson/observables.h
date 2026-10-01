#ifndef MYTN_SPIN_BOSON_OBSERVABLES_H
#define MYTN_SPIN_BOSON_OBSERVABLES_H

// Observables for spin-1/2 / spin-boson chains, for pure states (MPS) and for purified density
// matrices in the impurity geometry (bra on sites 1..N mirrored, ket on N+1..2N; see TEBD_edge_dissipation.h).
// Entanglement entropy is in core/entanglement.h.

#include <itensor/all.h>
#include "../core/entanglement.h"

using namespace std;
using namespace itensor;

// One value per site: <sigma^a_j> (a = "x","y","z") on spin-1/2 sites, <n_j> on boson sites.
vector<double>
measure_magnetization(MPS* psi, const SiteSet sites , string direction);

// Number of kinks sum_j <n_j (1 - n_{j+1})>, i.e. pairs |down_z>_j |up_z>_{j+1}.
double
measure_kink( MPS* psi, const SiteSet sites);

// <n_start n_j> for j = 1..N, n = |down_z><down_z|; connected (minus <n_start><n_j>) if connected = true.
vector<double>
measure_correlations(MPS* psi, const SiteSet sites, const int start, const bool connected);

// Tr(rho sigma^a_j) for the purified impurity state, on physical site q (q = -1: all sites).
// If compute_normalization, divides by Tr(rho).
vector<complex<double> >
measure_magnetization_impurity_first_site(MPS* psi , string direction, bool compute_normalization = true, int q = -1);

// Tr(rho O_q) for a local operator O (given on site 1 indices) on physical site q (q = -1: all sites).
vector<complex<double> >
measure_local_obs_impurity_first_site(MPS *psi , const ITensor O, bool compute_normalization, int q = -1);

// Tr(rho) of the purified impurity state.
double
compute_norm_purifed_impurity(MPS* psi);

// Tr(rho) of the purified impurity state, for MPS with quantum numbers.
double
compute_norm_purifed_impurity_QN(MPS* psi);

// Tr(rho O_q1 O_q2) on physical sites q1, q2; connected correlation if connected = true.
complex<double>
measure_correlation_impurity_first_site(MPS *psi , const ITensor O, bool compute_normalization, int q1, int q2, bool connected = false);

#endif
