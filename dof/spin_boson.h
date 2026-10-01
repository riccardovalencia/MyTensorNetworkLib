/**
 * @file spin_boson.h
 * @brief Mixed spin-boson systems: one bosonic mode (site 1) and N-1 spin-1/2 (sites 2..N),
 *        e.g. atoms in a single-mode cavity.
 */
#ifndef MYTN_DOF_SPIN_BOSON_H
#define MYTN_DOF_SPIN_BOSON_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/**
 * @brief Site set with a truncated boson on site 1 and spin-1/2 on sites 2..N.
 * @param N       Total number of sites.
 * @param max_occ Maximum occupation of the boson.
 */
SiteSet
make_spin_boson_sites(const int N , const int max_occ);

/**
 * @brief Doubled (bra-ket) site set of 2N sites for the purified density matrix of
 *        make_spin_boson_sites(N, max_occ).
 *
 * @verbatim
 *   s_{N-1} ... s_1  b  |  b  s_1 ... s_{N-1}
 *   ------ bra -------     ------ ket -------
 * @endverbatim
 * The bra is mirrored so that the two bosons, where dissipation acts, share the central bond.
 *
 * @param N       Number of physical sites (boson + spins).
 * @param max_occ Maximum occupation of the boson.
 */
SiteSet
make_purified_spin_boson_sites(const int N , const int max_occ);

/**
 * @brief Product state |n_photon> (x) |theta,phi>^(N-1), with the spin-coherent state
 *        |theta,phi> = cos(theta/2)|up_z> + e^{i phi} sin(theta/2)|down_z>.
 * The MPS is real for phi = 0 and complex otherwise.
 * @param sites    Site set from make_spin_boson_sites.
 * @param n_photon Fock state of the boson.
 * @param theta    Polar angle of the spin-coherent state on the Bloch sphere.
 * @param phi      Azimuthal angle.
 */
MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , double theta, double phi);

/**
 * @brief Same as above with site-dependent angles theta[j], phi[j] for the spins.
 */
MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , const vector<double> theta, const vector<double> phi);

#endif
