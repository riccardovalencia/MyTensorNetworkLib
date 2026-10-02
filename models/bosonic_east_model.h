/**
 * @file bosonic_east_model.h
 * @brief Bosonic quantum east model (bosonic east model): Hamiltonian MPOs and TEBD gates.
 *
 * Chain of `size` bosonic sites (ITensor Boson) with sigma^x_j = a_j + a^dag_j. The core
 * kinetically-constrained term is
 *
 *     H = -1/2 sum_j n_j ( J sigma^x_{j+1} - U n_{j+1} - 1 ),
 *
 * where an excitation on site j facilitates creation/annihilation on site j+1.
 * Variants fix the occupation n0 of a virtual site 0 (symmetry sector) and add a boundary
 * term (1 - symmetry)/2 n_size. Parameters common to several functions, not repeated below:
 * - sites    : Boson site set (at least `size` sites);
 * - size     : number of sites;
 * - n0       : occupation of the virtual site 0, coupled to site 1;
 * - symmetry : symmetry-sector eigenvalue entering the boundary term on site size;
 * - J        : facilitated hopping amplitude (J = e^{-s} in the (s, c) parametrization);
 * - U        : nearest-neighbour density-density coefficient (U = 1 - 2c).
 * The functions that read or name data files (energy variance, overlaps, io/) take s and c instead.
 */
#ifndef MYTN_MODELS_BOSONIC_EAST_MODEL_H
#define MYTN_MODELS_BOSONIC_EAST_MODEL_H

#include <itensor/all.h>
#include "../mps/gates.h"

using namespace std;
using namespace itensor;

/** @name Hamiltonian MPOs */
///@{

/**
 * @brief H = -n0/2 (J sigma^x_1 - U n_1 - 1) - 1/2 sum_j n_j (J sigma^x_{j+1} - U n_{j+1} - 1)
 *            + (1 - symmetry)/2 n_size.
 */
MPO
make_bosonic_east_model_mpo( const SiteSet sites, int size , int n0, double symmetry , double J, double U);

/**
 * @brief make_bosonic_east_model_mpo plus a drive Omega sum_{j=2}^{size-1} sigma^x_j.
 * @param Omega Drive amplitude (default 0.05).
 */
MPO
make_bosonic_east_model_mpo_with_drift( const SiteSet sites, int size , int n0, double symmetry , double J, double U, double Omega = 0.05);

/** @brief -make_bosonic_east_model_mpo (to target the highest-energy state with DMRG). */
MPO
make_bosonic_east_model_mpo_minus( const SiteSet sites, int size , int n0, double symmetry , double J, double U);

/**
 * @brief Bulk bosonic east model term plus on-site interaction and hopping (no site 0, no boundary term):
 *        H = -1/2 sum_j n_j (J sigma^x_{j+1} - U n_{j+1} - 1) + epsilon/2 sum_j n_j^2
 *            - t/2 sum_j (a^dag_j a_{j+1} + h.c.).
 * @param epsilon On-site interaction.
 * @param t       Hopping.
 */
MPO
make_bosonic_east_model_mpo_onsite_hopping( const SiteSet sites, int size , double J, double U , double epsilon, double t);

/**
 * @brief On-site interaction only, H = epsilon/2 sum_j n_j^2.
 * @param epsilon On-site interaction.
 */
MPO
make_bosonic_east_model_mpo_onsite( const SiteSet sites, int size , double epsilon);

/**
 * @brief Diagonal part without nearest-neighbour terms:
 *        H = n0/2 + 1/2 sum_{j<size} n_j + c/2 sum_j n_j^2 + (1 - symmetry)/2 n_size.
 * @param c Coefficient of the on-site term (the c of U = 1 - 2c).
 */
MPO
make_bosonic_east_model_mpo_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double c );

/**
 * @brief As make_bosonic_east_model_mpo, but site 1 is left untouched and plays the role of site 0
 *        (the chain starts at site 2). Used to find ground states in the sector fixed by n0
 *        starting from a state |n0>|psi>.
 */
MPO
make_bosonic_east_model_mpo_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double J, double U);

/**
 * @brief Bulk bosonic east model Hamiltonian without the virtual site 0 (n0 is not fixed),
 *        H = -1/2 sum_j n_j (J sigma^x_{j+1} - U n_{j+1}) + 1/2 sum_{j<size} n_j + (1 - symmetry)/2 n_size.
 * @note n0 is not a conserved quantity of this Hamiltonian.
 */
MPO
make_bosonic_east_model_mpo_n0_not_fixed( const SiteSet sites, int size , double symmetry , double J, double U);

/**
 * @brief MPO approximation of exp(-i dt H) (ITensor toExpH, first order in dt), with H as in
 *        make_bosonic_east_model_mpo_n0_not_fixed. Used by make_dressed_operator.
 * @param dt Time step.
 */
MPO
make_bosonic_east_model_evolution_mpo( const SiteSet sites, int size , double symmetry , double J, double U, double dt);

/**
 * @brief Energy variance <H^2> - <H>^2 of make_bosonic_east_model_mpo (J = e^{-s}, U = 1 - 2c) on psi.
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
 * @param s, c                bosonic east model parameters (they also name the data file).
 * @param symmetry_sector_dir Folder containing the symmetry-sector data files.
 * @throws ITError if the data file cannot be opened.
 */
double
compute_bosonic_east_model_energy_variance(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c, const string symmetry_sector_dir);
///@}

/** @name TEBD */
///@{

/**
 * @brief Bond term on (j, j+1) of the bulk Hamiltonian (the n_j term is assigned to the bond; the
 *        last bond also takes n_size).
 * @param j Left site of the bond (1 <= j < size).
 * @return Bond Hamiltonian.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian( const SiteSet sites , const int size , const double J , const double U , const int j );

/**
 * @brief As make_bosonic_east_model_bond_hamiltonian, including the coupling to the virtual site 0
 *        with occupation n0 on the first bond (the n_j terms are shared between the two bonds of each site).
 * @param j Left site of the bond.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian_n0( const SiteSet sites , const int size , const int n0, const double J , const double U , const int j );

/**
 * @brief As make_bosonic_east_model_bond_hamiltonian plus the effective non-hermitian term -i gamma/2 n_j^2
 *        of dephasing L_j = sqrt(gamma) n_j.
 * @param gamma Dephasing rate.
 * @param j     Left site of the bond.
 */
ITensor
make_bosonic_east_model_bond_hamiltonian_dephasing( const SiteSet sites , const int size , const double J , const double U , const double gamma, const int j );

/**
 * @brief Gates of one second-order Trotter step of the bulk Hamiltonian (apply with tebd_step).
 * @param dt       Time step.
 * @param dynamics "closed" (default) or "open" (adds dephasing with rate gamma).
 * @param gamma    Dephasing rate (only for "open").
 */
vector<TebdGate>
make_bosonic_east_model_gates(const SiteSet sites, const int size, const double dt, const double J, const double U, const string dynamics = "closed" , const double gamma = 0.);

///@}

#endif
