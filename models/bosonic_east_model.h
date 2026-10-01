/**
 * @file bosonic_east_model.h
 * @brief Bosonic quantum east model (bosonic east model): Hamiltonian MPOs and TEBD gates.
 *
 * Chain of `size` bosonic sites (ITensor Boson) with sigma^x_j = a_j + a^dag_j and
 * U = 1 - 2c. The core kinetically-constrained term is
 *
 *     H = -1/2 sum_j n_j ( e^{-s} sigma^x_{j+1} - U n_{j+1} - 1 ),
 *
 * where an excitation on site j facilitates creation/annihilation on site j+1.
 * Variants fix the occupation n0 of a virtual site 0 (symmetry sector) and add a boundary
 * term (1 - symmetry)/2 n_size. Parameters common to several functions:
 * - size     : number of sites;
 * - n0       : occupation of the virtual site 0, coupled to site 1;
 * - symmetry : symmetry-sector eigenvalue entering the boundary term on site size;
 * - s        : J = e^{-s} is the facilitated hopping amplitude;
 * - c        : density-density interaction, U = 1 - 2c.
 */
#ifndef MYTN_MODELS_BOSONIC_EAST_MODEL_H
#define MYTN_MODELS_BOSONIC_EAST_MODEL_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

/** @name Hamiltonian MPOs */
///@{

/**
 * @brief H = -n0/2 (e^{-s} sigma^x_1 - U n_1 - 1) - 1/2 sum_j n_j (e^{-s} sigma^x_{j+1} - U n_{j+1} - 1)
 *            + (1 - symmetry)/2 n_size.
 */
MPO
make_bosonic_east_model_mpo( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

/** @brief make_bosonic_east_model_mpo plus a drive 0.05 sum_{j=2}^{size-1} sigma^x_j. */
MPO
make_bosonic_east_model_mpo_with_drift( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

/** @brief -make_bosonic_east_model_mpo (to target the highest-energy state with DMRG). */
MPO
make_bosonic_east_model_mpo_minus( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

/**
 * @brief Bulk bosonic east model term plus on-site interaction and hopping (no site 0, no boundary term):
 *        H = -1/2 sum_j n_j (e^{-s} sigma^x_{j+1} - U n_{j+1} - 1) + epsilon/2 sum_j n_j^2
 *            - t/2 sum_j (a^dag_j a_{j+1} + h.c.).
 */
MPO
make_bosonic_east_model_mpo_onsite_hopping( const SiteSet sites, int size , double s, double c , double epsilon, double t);

/**
 * @brief On-site interaction only, H = epsilon/2 sum_j n_j^2.
 * @note n0, symmetry, s and c are not used.
 */
MPO
make_bosonic_east_model_mpo_onsite( const SiteSet sites, int size , int n0, double symmetry , double s, double c , double epsilon);

/**
 * @brief Diagonal part without nearest-neighbour terms:
 *        H = n0/2 + 1/2 sum_{j<size} n_j + c/2 sum_j n_j^2 + (1 - symmetry)/2 n_size.
 * @note s is not used.
 */
MPO
make_bosonic_east_model_mpo_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double s, double c );

/**
 * @brief As make_bosonic_east_model_mpo, but site 1 is left untouched and plays the role of site 0
 *        (the chain starts at site 2). Used to find ground states in the sector fixed by n0
 *        starting from a state |n0>|psi>.
 */
MPO
make_bosonic_east_model_mpo_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double s, double c);

/**
 * @brief Bulk bosonic east model Hamiltonian without the virtual site 0 (n0 is not fixed),
 *        H = -1/2 sum_j n_j (e^{-s} sigma^x_{j+1} - U n_{j+1}) + 1/2 sum_{j<size} n_j + (1 - symmetry)/2 n_size.
 * @note n0 is not a conserved quantity of this Hamiltonian.
 */
MPO
make_bosonic_east_model_mpo_n0_not_fixed( const SiteSet sites, int size , double symmetry , double s, double c);

/**
 * @brief MPO approximation of exp(-i dt H) (ITensor toExpH, first order in dt), with H as in
 *        make_bosonic_east_model_mpo_n0_not_fixed and hopping amplitude J. Used by make_dressed_operator.
 * @note n0 is not used.
 */
MPO
make_bosonic_east_model_evolution_mpo( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt);

/**
 * @brief Energy variance <H^2> - <H>^2 of make_bosonic_east_model_mpo on psi.
 *
 * The symmetry eigenvalue of the sector is read from
 * "<symmetry_sector_dir>/symmetry_sector<symmetry_sector>_maxcutoff30_s<s>_c<c>.dat"
 * (entry number cut_off_fock_space of the file).
 *
 * @param psi                 State.
 * @param sites               Boson site set.
 * @param size                Number of sites.
 * @param cut_off_fock_space  Fock-space cutoff, selects the entry of the data file.
 * @param n0                  Occupation of the virtual site 0.
 * @param symmetry_sector     Index of the symmetry sector.
 * @param s, c                bosonic east model parameters.
 * @param symmetry_sector_dir Folder containing the symmetry-sector data files.
 * @throws ITError if the data file cannot be opened.
 */
double
compute_bosonic_east_model_energy_variance(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c, const string symmetry_sector_dir);
///@}

/** @name TEBD */
///@{

/**
 * @brief Bond term on (j, j+1) of the bulk Hamiltonian with hopping J (the n_j term is assigned to the
 *        bond; the last bond also takes n_size).
 * @param J Facilitated hopping amplitude (e^{-s}).
 * @return Bond Hamiltonian.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian( const SiteSet sites , const int size , const double J , const double c , const int j );

/** @brief As make_bosonic_east_model_bond_hamiltonian, including the coupling to the virtual site 0 with occupation n0 on the first bond. */
ITensor
make_bosonic_east_model_bond_hamiltonian_n0( const SiteSet sites , const int size , const int n0, const double J , const double c , const int j );

/**
 * @brief As make_bosonic_east_model_bond_hamiltonian plus the effective non-hermitian term -i gamma/2 n_j^2
 *        of dephasing L_j = sqrt(gamma) n_j.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian_dephasing( const SiteSet sites , const int size , const double J , const double c , const double gamma, const int j );

/**
 * @brief As make_bosonic_east_model_bond_hamiltonian plus -i/2 L_j^dag L_j for arbitrary local jump operators.
 * @param Lj  Jump operators, one per site.
 * @param Ljd Their hermitian conjugates.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian_jumps( const SiteSet sites , const int size , const double J , const double c , const int j, vector<ITensor> &Lj, vector<ITensor> &Ljd);

/**
 * @brief Second-order Trotter step (BondGate list, for gateTEvol) of the bulk Hamiltonian.
 * @param dynamics "closed" (default) or "open" (adds dephasing with rate gamma).
 * @param gamma    Dephasing rate (only for "open").
 */
vector<BondGate>
make_bosonic_east_model_gates(const SiteSet sites, const int size, const double dt, const double J, const double c, const string dynamics = "closed" , const double gamma = 0.);

/** @brief As make_bosonic_east_model_gates with arbitrary local jump operators (see make_bosonic_east_model_bond_hamiltonian_jumps). */
vector<BondGate>
make_bosonic_east_model_gates_open(const SiteSet sites, const int size, const double dt, const double J, const double c , vector<ITensor> &Lj, vector<ITensor> &Ljd );
///@}

#endif
