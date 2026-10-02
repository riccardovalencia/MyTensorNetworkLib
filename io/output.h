/**
 * @file output.h
 * @brief Data management: output files.
 *
 * The callers choose the file names. The examples write each run into its own folder
 * <folder>/<run name>/ (make_run_directory), the run name encoding the parameters.
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
 * @brief Create the folder <folder>/<run_name>/ (and its parents) if needed, and copy the input
 *        file of the run into it as input.txt.
 * @code
 * string dir = make_run_directory("data", tinyformat::format("rydberg_N%d_D%d", N, max_dim), argv[1]);
 * ofstream out(dir + "observables.txt");
 * @endcode
 * @param folder     Parent folder (e.g. "data").
 * @param run_name   Name of the run, encoding all the parameters that change the results, so that
 *                   different inputs never share a folder.
 * @param input_file Input file of the run (not copied if empty).
 * @return The path of the folder, ending with '/'.
 */
string
make_run_directory(const string& folder, const string& run_name, const string& input_file = "");

/**
 * @brief Write one value per site: lines "j value", j = 1, 2, ...
 * @param file_name File to (over)write.
 * @param values    values[j-1] is written on line j.
 * @param precision Number of significant digits.
 */
void
write_site_values( const string file_name, const vector<double>& values, int precision );

/**
 * @brief Write one row per site: lines "j table[j-1][0] table[j-1][1] ...", j = 1, 2, ...
 *        (e.g. the Fock-state probabilities of measure_fock_probabilities).
 * @param file_name File to (over)write.
 * @param table     One row per site.
 * @param precision Number of significant digits.
 */
void
write_site_table( const string file_name, const vector<vector<double> >& table, int precision );

/**
 * @brief Write generating functions G_l(theta) of several blocks: one line per theta,
 *        "theta Re G_1 Im G_1 Re G_2 Im G_2 ...".
 * @param file_name File to (over)write (10 significant digits, with a header line).
 * @param theta Values of theta (e.g. from make_theta_grid).
 * @param G     G[l-1][k] = G_l(theta[k]) (e.g. from compute_generating_function).
 */
void
write_generating_function( const string file_name, const vector<double>& theta, const vector<vector<complex<double> > >& G );

/**
 * @brief Write the DMRG parameters of a bosonic east model run to "input.txt" (and to stderr).
 * @param size                   Number of sites.
 * @param s, c                   bosonic east model parameters.
 * @param simmetry_sector        Symmetry sector.
 * @param cut_off_fock_space     Fock-space cutoff.
 * @param scaling_bond_dimension Factor of the bond dimension between DMRG rounds.
 * @param bond_dimension         Initial bond dimension.
 * @param precision_dmrg         Convergence threshold.
 */
void
write_dmrg_input(int size , double s , double c , double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

/**
 * @brief As write_dmrg_input (same arguments), including the on-site interaction epsilon and the
 *        hopping t of make_bosonic_east_model_mpo_onsite_hopping.
 */
void
write_dmrg_input_hopping(int size , double s , double c , double epsilon, double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

#endif
