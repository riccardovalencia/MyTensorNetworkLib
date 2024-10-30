// #include "external_file.h"
#include "full_counting_statistics_termal.h"
#include <itensor/all.h>
#include <vector>
#include <string>
#include <complex>
#include <iomanip>

#include <itensor/all.h>
#include <iostream>
#include <sstream> // for ostringstream
#include <vector>
#include <string>
#include <iomanip>
#include <complex>


using namespace std;
using namespace itensor;

// Functions for computing and saving the generating function of the magnetization on M sites (\sum_{j=1}^M s_j^x)


void 
printing_generating_function( const stringstream *save_file ,  const int numberPoints , const int maxLength , vector<vector<double> > &G )
	{
	ofstream SaveFile( (*save_file).str() );		
	
	SaveFile << setprecision(10) << fixed;

	int col , row;
	
	double theta = -M_PI;
	
	for( col = 0 ; col <= numberPoints ; col++ )
		{
		SaveFile << theta << " ";
		for( row = 0 ; row < maxLength ; row++ ) SaveFile << G[row][col] << " ";
		if( col != (numberPoints - 1) ) SaveFile << "\n";	
		theta += theta_step( col , numberPoints);
		}
	SaveFile.close();	
	}

//----------------------------------------------------------------------
//return the value of theta where to evaluate generating function
double 
theta_step( int col , int numberPoints )
	{
	double thetaStep;
	if ( col <= numberPoints / 4 || col >= (3 * numberPoints / 4 - 1) ) thetaStep = ( M_PI - 1. ) / ( numberPoints / 4 );
	else thetaStep = 2. / (numberPoints / 2.);	
	return thetaStep;
	}
	
//----------------------------------------------------------------------
// measure of the generating function of probability distribution function of total magnetization in a certain
// subsystem at time fixed and size of subsystem on an input MPO (representing a mixed state).
// The subsystem is centered along the finite chain

void
generaring_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPO* psi , const SpinHalf sites )
	{
	int start;
	
	if( size % 2 == 0)	start = ( N/2 - size / 2 );					//if size is even we go to the left of the center
	else start = ( N/2 - (size + 1) / 2 ) ;
	// (*psi).position( start );	
	cerr << "(Simmetric) Measuring size = " << size + 1 << "\t"
		 << "start = " << start << endl;

	double theta = -M_PI ;
    
	for( int col = 0 ; col < numberPoints ; col++ ) 					//particular value of theta
		{

		ITensor Sx = op(sites, "Sx", start );
		ITensor Obs = expHermitian(Sx , theta * 1_i  ); 

		ITensor Meas;
		
        Meas = (*psi)(1);
        Meas *= op(sites,"Id",1);

        for( int j = 2 ; j < start ; j++) Meas *= (*psi)(j) * op(sites,"Id",j);
	
		if( size == 0 ) Meas *= (*psi)(start) * Obs ;
					
		else
			{
			// Meas *= (*psi)(start) * op(sites,"Id",start) ;
			for( int row = 0 ; row <= size ; row++ )
				{
				Sx = op(sites,"Sx", start + row);
				Obs = expHermitian(Sx , theta * 1_i  ); 
				Meas *= (*psi)(start + row) ;
				Meas *= Obs ;
				}

			}
        for( int j = start + size + 1 ; j <= N ; j++) Meas *= (*psi)(j) * op(sites,"Id",j);
		
		theta += theta_step( col , numberPoints);

		complex<double> SingleMeasure = eltC(Meas);
		singleGreal.push_back( SingleMeasure.real() );
		singleGimag.push_back( SingleMeasure.imag() );
		}
	}
	
//----------------------------------------------------------------------
//measure of generating function and saving of the information	
void	
measure_generating_function( MPO* rho , const SpinHalf sites , int N , int maxLength , int numberPoints , double hx , double hz)
	{
	vector<vector<double> > Greal;
	vector<vector<double> > Gimag;
	
	stringstream save_real_two , save_imag_two;
	
	save_real_two << "Termal_N" << N << "_hx_" << hx << "_hz_" << hz << "_GF_real.dat";
	save_imag_two << "Termal_N" << N << "_hx_" << hx << "_hz_" << hz << "_GF_imag.dat";

	if( fileExists( save_real_two.str() ) == false)
		{
		for( int size = 0 ; size < maxLength ; size++ )
			{			
			vector<double> singleGreal;
			vector<double> singleGimag;			
			generaring_function_sim_size( singleGreal , singleGimag , size , N , numberPoints , rho , sites );
			Greal.push_back( singleGreal );
			Gimag.push_back( singleGimag );
			}
		
		printing_generating_function( &save_real_two , numberPoints , maxLength , Greal );  
		printing_generating_function( &save_imag_two , numberPoints , maxLength , Gimag );  
		}
	else cout << "The files " << save_real_two.str() << " and "
			  << save_imag_two.str()
			  << " already exist." << endl;

	}	
