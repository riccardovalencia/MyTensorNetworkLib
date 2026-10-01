#ifndef MYTN_BOSONS_EXTERNAL_FILE_H
#define MYTN_BOSONS_EXTERNAL_FILE_H

// Console/file output for the bosonic quantum east model (bQEM) simulations.


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

#endif
