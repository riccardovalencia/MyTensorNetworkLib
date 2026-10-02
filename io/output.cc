/**
 * @file output.cc
 * @brief Implementation of output.h (interfaces documented in the header, logic commented here).
 */
#include "output.h"
#include <itensor/all.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;


// create_directories also creates the missing parents and does nothing if the folder exists
string
make_run_directory(const string& folder, const string& run_name, const string& input_file)
{
    string dir = folder + "/" + run_name + "/";
    std::filesystem::create_directories(dir);
    if(!input_file.empty())
        std::filesystem::copy_file(input_file, dir + "input.txt", std::filesystem::copy_options::overwrite_existing);
    return dir;
}


// two small tables: the result of the returned run, the energy of every restart
void
write_ground_state( const string& dir, const GroundState& ground_state )
{
    ofstream out(dir + "ground_state.txt");
    out << setprecision(14) << "# E . <H^2>-<H>^2 . converged\n";
    out << ground_state.energy << " " << ground_state.variance << " " << ground_state.converged << "\n";

    ofstream out_restarts(dir + "ground_state_restarts.txt");
    out_restarts << setprecision(14) << "# restart . E\n";
    for(size_t r = 0 ; r < ground_state.restart_energies.size() ; r++)
        out_restarts << r+1 << " " << ground_state.restart_energies[r] << "\n";
}


// x, then the values separated by spaces, then a newline (flushed, so that the file can be read during the run)
void
write_row( ostream& out, const double x, const vector<double>& values )
{
    out << x;
    for(double v : values) out << " " << v;
    out << endl;
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
        for(const vector<complex<double> >& G_l : G) out << " " << G_l[k].real() << " " << G_l[k].imag();
        out << "\n";
    }
}
