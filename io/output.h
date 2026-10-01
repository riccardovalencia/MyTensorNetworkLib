/**
 * @file output.h
 * @brief Data management: file names and output files.
 *
 * All files are written in the working directory.
 */
#ifndef MYTN_IO_OUTPUT_H
#define MYTN_IO_OUTPUT_H

#include <itensor/all.h>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;

/** @name Generic */
///@{

/**
 * @brief Write a matrix to file_name: one row per line, prefixed by the row index.
 * @param set_output_precision Number of digits.
 */
void
print_matrix( string file_name , vector<vector<double> > &matrix , int set_output_precision);
///@}

/** @name Spin chains and full counting statistics */
///@{

/** @brief Set sites_file = "sites_N<N>" and psi_file = "psi_N<N>_nstep". */
void
build_file_tebd( stringstream *sites_file , stringstream *psi_file , const int N );

/** @brief As build_file_tebd, plus "N<N>_GF_real", "N<N>_GF_imag", "N<N>_Moments_real", "N<N>_Moments_imag". */
void
build_file_full_counting( stringstream *sites_file , stringstream *psi_file , stringstream *save_real , stringstream *save_imag , stringstream *saveRealMoments , stringstream *saveImagMoments , const int N );

/** @brief As build_file_tebd, plus save_file = "N<N>_entropy.dat". */
void
build_file_entanglement_entropy( stringstream *sites_file , stringstream *psi_file , stringstream *save_file , const int N );

/** @brief Write the parameters of an Ising-chain simulation to "input.txt" (ITensor InputGroup format). */
void
print_input( int N , double J , double hx , double hz , double ttotal , double tstep , double nmeas , double bonddim );

/**
 * @brief Print time reached, timings, maximum bond dimension and half-chain entanglement entropy.
 * @return The half-chain entanglement entropy (natural log).
 */
double
print_info( time_t time_elapsed_step , time_t time_elapsed_total , MPS *psi , const int N , const int nmeas , const int n , const double tstep );

/**
 * @brief Write a generating function G[size][k] to the file named *save_file: one line per theta,
 *        with columns theta, G(size = 0), G(size = 1), ...
 */
void
printing_generating_function( const stringstream *save_file , const int numberPoints , const int maxLength , vector<vector<double> > &G );
///@}

/** @name Bosonic quantum east model */
///@{

/** @brief Write the DMRG parameters of a bQEM run to "input.txt". */
void
print_input_dmrg(int size , double s , double c , double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

/** @brief As print_input_dmrg, including on-site interaction epsilon and hopping t. */
void
print_input_dmrg_hopping(int size , double s , double c , double epsilon, double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

/**
 * @brief Write the occupations (site, <n_j>) to "occupation_number_size<size>_s<100s>_c<100c>_cutoff<..>_n0<..>_chi<bond_dimension>...".
 */
void
print_occupation_number( vector<double> &occupation_number, int size, double s, double c , int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision);

/** @brief As print_occupation_number for one disorder realization at time `time`; file name includes gamma, t and index. */
void
print_occupation_number_realization_disorder( vector<double> &occupation_number, int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma);

/** @brief As print_occupation_number for a mean-field state ("mean_field_occupation_number_size...dat"). */
void
print_occupation_number_mean_field( vector<double> &occupation_number, int size, double s, double c , int cut_off_fock_space , int n0 , int set_output_precision);

/** @brief As print_occupation_number for the excited state number number_state ("..._E<number_state>..."). */
void
print_occupation_number_excited_states( vector<double> &occupation_number, int size, double s, double c , int cut_off_fock_space , int n0 , int number_state , int set_output_precision);

/**
 * @brief Write the Fock-state probabilities (output of measure_projector_all_sites), one line per site,
 *        to "projector_fock_space_size<size>_s<100s>_c<100c>_cutoff<..>_n0<..>_chi<bond_dimension>...".
 */
void
print_projector_fockspace( vector<vector<double> > &projector_all_sites , int size, double s, double c , int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision);

/** @brief As print_projector_fockspace for one disorder realization at time `time`. */
void
print_projector_fockspace_realization_disorder( vector<vector<double> > &projector_all_sites , int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma);

/** @brief As print_projector_fockspace for a mean-field state. */
void
print_projector_fockspace_mean_field( vector<vector<double> > &projector_all_sites , int size, double s, double c , int cut_off_fock_space , int n0 , int set_output_precision);

/** @brief As print_projector_fockspace for the excited state number number_state. */
void
print_projector_fockspace_excited_states( vector<vector<double> > &projector_all_sites , int size, double s, double c , int cut_off_fock_space , int n0 , int number_state, int set_output_precision);
///@}

#endif
