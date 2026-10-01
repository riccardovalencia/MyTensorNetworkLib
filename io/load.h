/**
 * @file load.h
 * @brief Data management: load states written by previous simulations
 *        (bosonic quantum east model folder layout).
 */
#ifndef MYTN_IO_LOAD_H
#define MYTN_IO_LOAD_H

#include <itensor/all.h>
#include <complex>
#include <string>
#include <tuple>

using namespace std;
using namespace itensor;

/**
 * @brief Load the bQEM ground state with the largest available bond dimension.
 *
 * Looks in "<results_dir><size>_cutoff<lambda>/mmGcbQEM_size..._sector<symmetry_sector>_c<c>/
 * mmGcbQEM_size..._s<s>_v<version>/" for "sites_file_n0<n0>" and "ground_state_file_n0<n0>_chi<chi>",
 * with chi = bond_dimension * scaling_bond_dimension^k, and returns the first version whose energy
 * variance is below 1e-8.
 *
 * @param symmetry_sector_dir Folder of the symmetry-sector data files (see compute_variance_hamiltonian_bqem).
 * @return {state, sites, energy variance}.
 */
tuple<MPS, Boson, double>
search_ground_state_max_bond_chi(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension, const string symmetry_sector_dir);

/** @brief As search_ground_state_max_bond_chi, for the folder layout without the "_v<version>" suffix. */
tuple<MPS, Boson, double>
search_ground_state_max_bond_chi_no_v(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension, const string symmetry_sector_dir);

/**
 * @brief Load a state prepared by adiabatic dressing ("sites_size..." and "psi_file_size..." in results_dir).
 * @param state_choice 0: super-coherent state, 3: cat state.
 * @param beta         Duration of the ramp (in the file name).
 * @return {state, sites}. Exits if the files are not found.
 */
tuple<MPS, Boson>
search_state_adiabatic_coherent(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta );

#endif
