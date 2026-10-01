#include "external_file.h"
#include "observables.h"
#include <itensor/all.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;



//----------------------------------------------------------------------
//print input DMRG in Bosonic Quantum East Model
void
print_input_DMRG(int size , double s , double c ,double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
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
print_input_DMRG_hopping(int size , double s , double c , double epsilon , double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg )
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


//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model) of a single realization of disorder

void
print_occupation_number_realization_disorder( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma)
{
	// stringstream file_occupation_number ; 
	// file_occupation_number << "occupation_number_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_t" << int(time) << "_index" << index << ".dat";
	string file_occupation_number ; 
	file_occupation_number = tinyformat::format("occupation_number_size%d_s%.2f_c%.2f_cutoff%d_n0%d_gamma%.3f_t%.3f_index%d.dat", size, s, c, cut_off_fock_space, n0, gamma, time, index);
	
	ofstream save_file_occupation_number( file_occupation_number );
	save_file_occupation_number << setprecision(set_output_precision);
	for(int j = 1 ; j <= size ; j++)
		{
		save_file_occupation_number << j << " " << occupation_number[j-1] << endl;
		}
	save_file_occupation_number.close();
}

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model) of a single realization of disorder

void
print_projector_fockspace_realization_disorder( vector<vector<double> > &projector_all_sites ,  int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma)
{
	// stringstream file_projector_fock;
	// file_projector_fock << "projector_fock_space_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_chi" << bond_dimension << "_index" << index << ".dat";
	string file_projector_fock ; 
	file_projector_fock = tinyformat::format("projector_fock_space_size%d_s%.2f_c%.2f_cutoff%d_n0%d_gamma%.3f_t%.3f_index%d.dat", size, s, c, cut_off_fock_space, n0, gamma, time, index);

	ofstream save_file_projector_fock( file_projector_fock );
	save_file_projector_fock << setprecision(set_output_precision);	

	for(int j = 1 ; j <= size ; j++) 
		{
		save_file_projector_fock << j ;
		for(int n = 0 ; n <= cut_off_fock_space ; n++)
			{
			save_file_projector_fock <<	" " << projector_all_sites[j-1][n];
			}
		save_file_projector_fock << "\n";
		}

	save_file_projector_fock.close();
		
}

void 
print_matrix( string file_name , vector<vector<double> > &matrix ,int set_output_precision)
{	
	ofstream save_file( file_name );
	save_file << setprecision(set_output_precision);	

	int nrow = matrix.size();
	int ncol = matrix[0].size();

	for(int j = 0 ; j < nrow ; j++) 
		{
		save_file << j ;
		for(int n = 0 ; n < ncol ; n++)
			{
			save_file <<	" " << matrix[j][n];
			}
		save_file << "\n";
		}

	save_file.close();

}





//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model)

void
print_occupation_number( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision)
{
	stringstream file_occupation_number ; 
	file_occupation_number << "occupation_number_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_chi" << bond_dimension << ".dat";
	ofstream save_file_occupation_number( file_occupation_number.str() );
	save_file_occupation_number << setprecision(set_output_precision);
	for(int j = 1 ; j <= size ; j++)
		{
		save_file_occupation_number << j << " " << occupation_number[j-1] << endl;
		}
	save_file_occupation_number.close();
}

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model)

void
print_projector_fockspace( vector<vector<double> > &projector_all_sites ,  int size, double s, double c , int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision)
{
	stringstream file_projector_fock;
	file_projector_fock << "projector_fock_space_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_chi" << bond_dimension << ".dat";
	ofstream save_file_projector_fock( file_projector_fock.str() );
	save_file_projector_fock << setprecision(set_output_precision);	

	for(int j = 1 ; j <= size ; j++) 
		{
		save_file_projector_fock << j ;
		for(int n = 0 ; n <= cut_off_fock_space ; n++)
			{
			save_file_projector_fock <<	" " << projector_all_sites[j-1][n];
			}
		save_file_projector_fock << "\n";
		}

	save_file_projector_fock.close();
		
}


//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model)											MEAN FIELD

void
print_occupation_number_mean_field( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int set_output_precision)
{
	stringstream file_occupation_number ; 
	// file_occupation_number << "mean_field_occupation_number_size" << size << "_s" << int(s*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
	file_occupation_number << "mean_field_occupation_number_size" << size << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
	ofstream save_file_occupation_number( file_occupation_number.str() );
	save_file_occupation_number << setprecision(set_output_precision);
	for(int j = 1 ; j <= size ; j++)
		{
		save_file_occupation_number << j << " " << occupation_number[j-1] << endl;
		}
	save_file_occupation_number.close();
}

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model)			MEAN FIELD

void
print_projector_fockspace_mean_field( vector<vector<double> > &projector_all_sites ,  int size, double s,  double c , int cut_off_fock_space , int n0 , int set_output_precision)
{
	stringstream file_projector_fock;
	// file_projector_fock << "mean_field_projector_fock_space_size" << size << "_s" << int(s*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
	file_projector_fock << "mean_field_projector_fock_space_size" << size <<  "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
	ofstream save_file_projector_fock( file_projector_fock.str() );
	save_file_projector_fock << setprecision(set_output_precision);	

	for(int j = 1 ; j <= size ; j++) 
		{
		save_file_projector_fock << j ;
		for(int n = 0 ; n <= cut_off_fock_space ; n++)
			{
			save_file_projector_fock <<	" " << projector_all_sites[j-1][n];
			}
		save_file_projector_fock << "\n";
		}

	save_file_projector_fock.close();
		
}

void
print_occupation_number_excited_states( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int number_state, int set_output_precision)
{
	stringstream file_occupation_number ; 
	file_occupation_number << "occupation_number_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_E" << number_state << ".dat";
	ofstream save_file_occupation_number( file_occupation_number.str() );
	save_file_occupation_number << setprecision(set_output_precision);
	for(int j = 1 ; j <= size ; j++)
		{
		save_file_occupation_number << j << " " << occupation_number[j-1] << endl;
		}
	save_file_occupation_number.close();
}

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model)

void
print_projector_fockspace_excited_states( vector<vector<double> > &projector_all_sites ,  int size, double s, double c , int cut_off_fock_space , int n0 , int number_state, int set_output_precision)
{
	stringstream file_projector_fock;
	file_projector_fock << "projector_fock_space_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << "_E" << number_state << ".dat";
	ofstream save_file_projector_fock( file_projector_fock.str() );
	save_file_projector_fock << setprecision(set_output_precision);	

	for(int j = 1 ; j <= size ; j++) 
		{
		save_file_projector_fock << j ;
		for(int n = 0 ; n <= cut_off_fock_space ; n++)
			{
			save_file_projector_fock <<	" " << projector_all_sites[j-1][n];
			}
		save_file_projector_fock << "\n";
		}

	save_file_projector_fock.close();
		
}



//
//
//
//
