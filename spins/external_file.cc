#include "external_file.h"
#include "observables.h"
#include "full_counting_statistics.h"
#include <itensor/all.h>
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

// Contains functions for handling specific-purpose input/output data files.


//----------------------------------------------------------------------
//files to save the information about sites and the state
void 
build_file_TEBD( stringstream *sites_file ,  stringstream *psi_file , const int N )
	{
	*sites_file << "sites_N" << N ;
	*psi_file << "psi_N" << N << "_nstep";
	}
	
//----------------------------------------------------------------------
//files to save information for measuring generating function and moments of a certain full counting statistics
void
build_file_full_counting( stringstream *sites_file ,  stringstream *psi_file , stringstream *save_real , stringstream *save_imag ,  stringstream *saveRealMoments ,  stringstream *saveImagMoments , const int N )
	{
		*sites_file << "sites_N" << N ;
		*psi_file << "psi_N" << N << "_nstep";
		*save_real << "N" << N <<  "_GF_real";
		*save_imag << "N" << N <<  "_GF_imag";
		*saveRealMoments << "N" << N <<  "_Moments_real";
		*saveImagMoments << "N" << N <<  "_Moments_imag";
	}	
	
//----------------------------------------------------------------------
//file to save entanglement entropy of a one dimensional system
void
build_file_entanglement_entropy( stringstream *sites_file ,  stringstream *psi_file , stringstream *save_file , const int N )
	{
		*sites_file << "sites_N" << N ;
		*psi_file << "psi_N" << N << "_nstep";
		*save_file << "N" << N << "_entropy.dat";
	}	
	
//----------------------------------------------------------------------
//print input for time evolution in Ising model with longitudinal and transversal magnetic fields
void
print_input(int N , double J , double hx , double hz , double ttotal , double tstep , double nmeas , double bonddim )
	{
	string input_file = "input.txt";
	cout << "FileName = " << input_file << endl;
	ofstream SaveInput( input_file.c_str() );
	SaveInput << "input" << "\n"
			  << "{\n"
			  << "N = " << N << ";\n"
			  << "J = " << J << ";\n"
			  << "hx = " << hx << ";\n"
			  << "hz = " << hz << ";\n"
			  << "ttotal = " << ttotal << ";\n"
			  << "tstep = " << tstep << ";\n"
			  << "nmeas = " << nmeas << ";\n"
			  << "bond_dim = " << bonddim << ";\n"
			  << "maxLength = " << N/2 << ";\n"
			  << "numberPoints = " << 200 << ";\n"
			  << "}" << endl;
	SaveInput.close();
	}

//----------------------------------------------------------------------
//print information during time evolution: time reached, time needed to do a single step, time needed in total, max bond dimension, entanglement entropy
double
print_info( time_t time_elapsed_step , time_t  time_elapsed_total , MPS *psi , const int N , const int nmeas , const int n , const double tstep )
	{
	double entropy = entanglement_entropy( psi , N , N/2 );
	cout << "Time evolution : " << n * nmeas * tstep << "\n"
		 << "Single Step Time : " << time_elapsed_step << "\n"
		 << "Total Time : " <<  time_elapsed_total << "\n"
		 << "Max bond dimension : " << maxM( *(psi) ) << "\n"
		 << "Entanglement entropy : " << entropy << "\n" << endl;
	return entropy;
	}
	
//----------------------------------------------------------------------
//print at time fixed the generating function of a certain probability distribution
//nb function put in full_counting_statistics.* because otherwise there are compiling problems


