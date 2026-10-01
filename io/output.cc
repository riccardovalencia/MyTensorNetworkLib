/**
 * @file output.cc
 * @brief Implementation of output.h (the functions are documented in the header).
 */
#include "output.h"
#include <itensor/all.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;


//----------------------------------------------------------------------
//print input DMRG in Bosonic Quantum East Model
void
write_dmrg_input(int size , double s , double c ,double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
	{
	string input_file = "input.txt";
	cout << "FileName = " << input_file << endl;
	ofstream SaveInput( input_file.c_str() );

	SaveInput << "Input: \n \n"
			<< "Physical quantities.\n"
			<< "size : " << size << "\n"
			<< "s : " << s  << "\n"
			<< "c : " << c << "\n"
			<< "simmetry_sector : " << simmetry_sector << "\n \n"	
			<< "Numerical quantities.\n"
			<< "cut_off_fock_space : " << cut_off_fock_space << "\n"
			<< "starting bond dimension (then it is increased of a factor " << scaling_bond_dimension << " in the folowing DMRG) : " << bond_dimension << "\n"
			<< "precision dmrg : " << precision_dmrg << endl;

	SaveInput.close();
	
	cerr << "\n\nInput: \n \n"
		<< "Physical quantities.\n"
		<< "size : " << size << "\n"
		<< "s : " << s  << "\n"
		<< "c : " << c << "\n"
		<< "simmetry_sector : " << simmetry_sector << "\n \n"	
		<< "Numerical quantities.\n"
		<< "cut_off_fock_space : " << cut_off_fock_space << "\n"
        << "starting bond dimension (then it is increased of a factor " << scaling_bond_dimension << " in the folowing DMRG) : " << bond_dimension << "\n"
		<< "precision dmrg : " << precision_dmrg << endl;
	
	}


//----------------------------------------------------------------------
//print input DMRG in Bosonic Quantum East Model with hopping
void
write_dmrg_input_hopping(int size , double s , double c , double epsilon , double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
	{
	string input_file = "input.txt";
	cout << "FileName = " << input_file << endl;
	ofstream SaveInput( input_file.c_str() );

	SaveInput << "Input: \n \n"
			<< "Physical quantities.\n"
			<< "size : " << size << "\n"
			<< "s : " << s  << "\n"
			<< "c : " << c << "\n"
			<< "epsilon : " << epsilon << "\n"
			<< "t : " << t << "\n"
			<< "simmetry_sector : " << simmetry_sector << "\n \n"	
			<< "Numerical quantities.\n"
			<< "cut_off_fock_space : " << cut_off_fock_space << "\n"
			<< "starting bond dimension (then it is increased of a factor " << scaling_bond_dimension << " in the folowing DMRG) : " << bond_dimension << "\n"
			<< "precision dmrg : " << precision_dmrg << endl;

	SaveInput.close();
	
	cerr << "\n\nInput: \n \n"
		<< "Physical quantities.\n"
		<< "size : " << size << "\n"
		<< "s : " << s  << "\n"
		<< "c : " << c << "\n"
		<< "epsilon : " << epsilon << "\n"
		<< "t : " << t << "\n"
		<< "simmetry_sector : " << simmetry_sector << "\n \n"	
		<< "Numerical quantities.\n"
		<< "cut_off_fock_space : " << cut_off_fock_space << "\n"
        << "starting bond dimension (then it is increased of a factor " << scaling_bond_dimension << " in the folowing DMRG) : " << bond_dimension << "\n"
		<< "precision dmrg : " << precision_dmrg << endl;
	
	}


void
write_site_values( const string file_name, const vector<double>& values, int precision )
{
    ofstream out(file_name);
    out << setprecision(precision);
    for(size_t j = 0 ; j < values.size() ; j++) out << j+1 << " " << values[j] << "\n";
}


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
