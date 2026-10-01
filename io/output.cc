/**
 * @file output.cc
 * @brief Implementation of output.h (interfaces documented in the header, logic commented here).
 */
#include "output.h"
#include <itensor/all.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
using namespace itensor;


//----------------------------------------------------------------------
// Summary of the DMRG parameters of a bosonic quantum east model run; extra_physical lists
// additional physical parameters (name, value). It is written to "input.txt" and to stderr.

static void
write_dmrg_summary(int size , double s , double c , const vector<pair<string,double> >& extra_physical, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
{
	ostringstream text;
	text << "Input: \n \n"
	     << "Physical quantities.\n"
	     << "size : " << size << "\n"
	     << "s : " << s  << "\n"
	     << "c : " << c << "\n";
	for(const auto& [name, value] : extra_physical) text << name << " : " << value << "\n";
	text << "simmetry_sector : " << simmetry_sector << "\n \n"
	     << "Numerical quantities.\n"
	     << "cut_off_fock_space : " << cut_off_fock_space << "\n"
	     << "starting bond dimension (then it is increased of a factor " << scaling_bond_dimension << " in the folowing DMRG) : " << bond_dimension << "\n"
	     << "precision dmrg : " << precision_dmrg << "\n";

	string input_file = "input.txt";
	cout << "FileName = " << input_file << endl;
	ofstream(input_file) << text.str();
	cerr << "\n\n" << text.str();
}


// summary without extra physical parameters
void
write_dmrg_input(int size , double s , double c ,double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
{
	write_dmrg_summary(size, s, c, {}, simmetry_sector, cut_off_fock_space, scaling_bond_dimension, bond_dimension, precision_dmrg);
}


// summary with epsilon and t
void
write_dmrg_input_hopping(int size , double s , double c , double epsilon , double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
{
	write_dmrg_summary(size, s, c, {{"epsilon", epsilon}, {"t", t}}, simmetry_sector, cut_off_fock_space, scaling_bond_dimension, bond_dimension, precision_dmrg);
}


// one line per entry, prefixed with its 1-based site
void
write_site_values( const string file_name, const vector<double>& values, int precision )
{
    ofstream out(file_name);
    out << setprecision(precision);
    for(size_t j = 0 ; j < values.size() ; j++) out << j+1 << " " << values[j] << "\n";
}


// one line per row, prefixed with its 1-based site
void
write_site_table( const string file_name, const vector<vector<double> >& table, int precision )
{
    ofstream out(file_name);
    out << setprecision(precision);
    for(size_t j = 0 ; j < table.size() ; j++)
    {
        out << j+1;
        for(double x : table[j]) out << " " << x;
        out << "\n";
    }
}


// header, then one line per theta with real and imaginary parts of every block
void
write_generating_function( const string file_name, const vector<double>& theta, const vector<vector<complex<double> > >& G )
{
    ofstream out(file_name);
    out << setprecision(10) << "# theta . Re G_l . Im G_l  (l = 1, 2, ...)\n";
    for(size_t k = 0 ; k < theta.size() ; k++)
    {
        out << theta[k];
        for(const auto& G_l : G) out << " " << G_l[k].real() << " " << G_l[k].imag();
        out << "\n";
    }
}
