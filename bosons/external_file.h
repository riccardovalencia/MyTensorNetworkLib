#ifndef EXTERNAL_FILE_H
#define EXTERNAL_FILE_H


#include <itensor/all.h>
#include <iostream>
#include <sstream> // for ostringstream
#include <vector>
#include <string>
#include <iomanip>
#include <complex>
#include <ctime>

using namespace std;
using namespace itensor;


//----------------------------------------------------------------------
//print input DMRG in Bosonic Quantum East Model
void
print_input_DMRG(int size , double s , double c ,double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

//----------------------------------------------------------------------
//print input DMRG in Bosonic Quantum East Model with hopping
void
print_input_DMRG_hopping(int size , double s , double c , double epsilon, double t, double simmetry_sector , int cut_off_fock_space , int scaling_bond_dimension , int bond_dimension , double precision_dmrg );

//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model) of a single realization of disorder
void
print_occupation_number_realization_disorder( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma);

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model) of a single realization of disorder
void
print_projector_fockspace_realization_disorder( vector<vector<double> > &projector_all_sites ,  int size, double s, double c , int cut_off_fock_space , int n0 , double time, int set_output_precision, int index, double gamma);

//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model)
void
print_occupation_number( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision);

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model)
void
print_projector_fockspace( vector<vector<double> > &projector_all_sites ,  int size, double s, double c ,int cut_off_fock_space , int n0 , int bond_dimension, int set_output_precision);

//
void 
print_matrix( string file_name , vector<vector<double> > &matrix ,int set_output_precision);

//----------------------------------------------------------------------
//print occupation_number (bosonic quantum east model) MEAN FIELD
void
print_occupation_number_mean_field( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int set_output_precision);

//----------------------------------------------------------------------
//print projector over fock space for each physical site (bosonic quantum east model) MEAN FIELD
void
print_projector_fockspace_mean_field( vector<vector<double> > &projector_all_sites ,  int size, double s, double c , int cut_off_fock_space , int n0 , int set_output_precision);

//print occupation_number excited states (bosonic quantum east model)
void
print_occupation_number_excited_states( vector<double> &occupation_number,  int size, double s, double c , int cut_off_fock_space , int n0 , int number_state , int set_output_precision);

//----------------------------------------------------------------------
//print projector over fock space for each physical site of excited states (bosonic quantum east model)
void
print_projector_fockspace_excited_states( vector<vector<double> > &projector_all_sites ,  int size, double s, double c ,int cut_off_fock_space , int n0 , int number_state, int set_output_precision);


//----------------------------------------------------------------------
//files to save the information about sites and the state
void 
build_file_TEBD( stringstream * ,  stringstream * , const int  );

//----------------------------------------------------------------------
//files to save information for measuring generating function and moments of a certain full counting statistics
void
build_file_full_counting( stringstream * ,  stringstream * , stringstream * , stringstream * ,  stringstream * ,  stringstream * , const int );

//----------------------------------------------------------------------
//file to save entanglement entropy of a one dimensional system
void
build_file_entanglement_entropy( stringstream * ,  stringstream * , stringstream * , const int );
	
//----------------------------------------------------------------------
//print input for time evolution in Ising model with longitudinal and transversal magnetic fields
void
print_input(int , double , double , double , double , double , double , double );

//----------------------------------------------------------------------
//print information during time evolution: time reached, time needed to do a single step, time needed in total, max bond dimension, entanglement entropy
//double
//print_info( time_t , time_t  , MPS * , const int , const int , const int , const double );

//----------------------------------------------------------------------
//print at time fixed the generating function of a certain probability distribution


#endif
