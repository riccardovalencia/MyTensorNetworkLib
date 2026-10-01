/**
 * @file output.h
 * @brief Data management: output files.
 *
 * The callers choose the file names (e.g. with tinyformat::format, encoding the parameters).
 */
#ifndef MYTN_IO_OUTPUT_H
#define MYTN_IO_OUTPUT_H

#include <itensor/all.h>
#include <complex>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;

/**
 * @brief Write one value per site: lines "j value", j = 1, 2, ...
 * @param precision Number of significant digits.
 */
void
write_site_values( const string file_name, const vector<double>& values, int precision );

/**
 * @brief Write one row per site: lines "j table[j-1][0] table[j-1][1] ...", j = 1, 2, ...
 *        (e.g. the Fock-state probabilities of measure_fock_probabilities).
 * @param precision Number of significant digits.
 */
void
write_site_table( const string file_name, const vector<vector<double> >& table, int precision );

/**
 * @brief Write generating functions G_l(theta) of several blocks: one line per theta,
 *        "theta Re G_1 Im G_1 Re G_2 Im G_2 ...".
 * @param theta Values of theta (e.g. from make_theta_grid).
 * @param G     G[l-1][k] = G_l(theta[k]) (e.g. from compute_generating_function).
 */
void
write_generating_function( const string file_name, const vector<double>& theta, const vector<vector<complex<double> > >& G );

/** @brief Write the DMRG parameters of a bosonic east model run to "input.txt". */
void
write_dmrg_input(int size , double s , double c , double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

/** @brief As write_dmrg_input, including on-site interaction epsilon and hopping t. */
void
write_dmrg_input_hopping(int size , double s , double c , double epsilon, double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

#endif
