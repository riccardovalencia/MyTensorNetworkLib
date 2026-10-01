#include "full_counting_statistics.h"
#include <itensor/all.h>
#include <vector>
#include <string>
#include <complex>
#include <iomanip>
#include <iostream>
#include <sstream> // for ostringstream


using namespace std;
using namespace itensor;

void 
printing_generating_function( const stringstream *save_file ,  const int numberPoints , const int maxLength , vector<vector<double> > &G )
	{
	ofstream SaveFile( (*save_file).str() );		
	
	SaveFile << setprecision(10) << fixed;

	int col , row;
	
	double theta = -M_PI;
	
	for( col = 0 ; col < numberPoints ; col++ )
		{
		SaveFile << theta << " ";
		for( row = 0 ; row < maxLength ; row++ ) SaveFile << G[row][col] << " ";
		if( col != (numberPoints - 1) ) SaveFile << "\n";	
		theta += theta_step( col , numberPoints);
		}
	SaveFile.close();	
	}

//----------------------------------------------------------------------
//return theta where to evaluate generating function
double 
theta_step( int col , int numberPoints )
	{
	double thetaStep;
	if ( col <= numberPoints / 4 || col >= (3 * numberPoints / 4 - 1) ) thetaStep = ( M_PI - 1. ) / ( numberPoints / 4 );
	else thetaStep = 2. / (numberPoints / 2.);	
	return thetaStep;
	}
	
//----------------------------------------------------------------------
//----------------------------------------------------------------------
// measure of the generating function of probability distribution function of total magnetization in a certain
// subsystem at time fixed and size of subsystem on an input MPS (representing a pure state).
// The subsystem is centered along the finite chain

void
generaring_function_sim_size( vector<double> &singleGreal , vector<double> &singleGimag , int size , int N , int numberPoints , MPS* psi , const SpinHalf sites )
	{
	int start;
	
	if( size % 2 == 0)	start = ( N/2 - size / 2 );					//if size is even we go to the left of the center
	else start = ( N/2 - (size + 1) / 2 ) ;
	(*psi).position( start );	
	cout << "(Simmetric) Measuring size = " << size + 1 << "\t"
		 << "start = " << start << endl;

	double theta = -M_PI ;

	for( int col = 0 ; col < numberPoints ; col++ ) 					//particular value of theta
		{
		ITensor Sx = sites.op( "Sx", start );
		ITensor Obs = expHermitian(Sx , theta * 1_i  ); 
		ITensor Meas;
		
		if( size == 0 )
			{
			Meas = (*psi)(start) * Obs * dag( prime( (*psi)(start) , "Site" ) );
			}		
		else
			{
			Index ir = commonIndex( (*psi)(start) , (*psi)(start + 1) , "Link");
			Meas = (*psi)(start) * Obs * dag( prime( prime( (*psi)(start) , "Site" ) , ir ) );
			
			for( int row = 1 ; row < size ; row++ )
				{
				Sx = sites.op("Sx", start + row);
				Obs = expHermitian(Sx , theta * 1_i  ); 
				Meas *= (*psi)(start + row) ;
				Meas *= Obs ;
				Meas *= dag( prime( (*psi)(start + row) ) );
				}
			
			Meas *= (*psi)( start + size );
			Index il = commonIndex( (*psi)( start + size ), (*psi)( start + size -1 ), "Link");
	
			Sx = sites.op("Sx", start + size );
			Obs = expHermitian(Sx , theta * 1_i  ); 
			Meas *= Obs;
			Meas *= dag( prime( prime( (*psi)( start + size ), il ) , "Site") );
			}
			
		theta += theta_step( col , numberPoints);
	
		complex<double> SingleMeasure = eltC(Meas);
		singleGreal.push_back( SingleMeasure.real() );
		singleGimag.push_back( SingleMeasure.imag() );
		}
	}
	
//----------------------------------------------------------------------
//measure of generating function and saving of the information	
void	
measure_generating_function( MPS* psi , const SpinHalf sites , int N , int n , int maxLength , int numberPoints , stringstream* save_real , stringstream* save_imag )
	{
	vector<vector<double> > Greal;
	vector<vector<double> > Gimag;
	
	stringstream save_real_two , save_imag_two;
	save_real_two << (*save_real).str();
	save_imag_two << (*save_imag).str();
	
	save_real_two << n << ".dat";
	save_imag_two << n << ".dat";				
		
	if( fileExists( save_real_two.str() ) == false)
		{
		for( int size = 0 ; size < maxLength ; size++ )
			{			
			vector<double> singleGreal;
			vector<double> singleGimag;			
			generaring_function_sim_size( singleGreal , singleGimag , size , N , numberPoints , psi , sites );
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

//----------------------------------------------------------------------
//return the MPO of the total magnetization of a spin-1/2 system of a subsystem of size "size"
MPO
build_totalSx( const SpinHalf sites , const int start ,  const int size )
	{
	
	AutoMPO ampo(sites);
		
	for(int j = start ; j <= start + size ; j++) ampo += 1. , "Sx" , j ;

	MPO totalSx = toMPO(ampo);	
	
	return totalSx;
	}	
	
//----------------------------------------------------------------------
// measure the first 4 moments of the full counting statistics of the total magnetization of a system of N spin-1/2 system
// on a subsystem of size maxLength centered along the chain of size 

void
measuring_moments( MPS *psi , const SpinHalf sites , const int N , const int n , const int maxLength , stringstream* saveRealMoments , stringstream* saveImagMoments )
	{
	int start;
	
	vector<vector<double> > MReal;
	vector<vector<double> > MImag; 
	
	stringstream saveRealMomentsTwo , saveImagMomentsTwo;
	saveRealMomentsTwo << (*saveRealMoments).str() << n << ".dat";
	saveImagMomentsTwo << (*saveImagMoments).str() << n << ".dat";	
	
	
	if( fileExists( saveRealMomentsTwo.str() ) == false)
		{
		for( int size = 0 ; size < maxLength ; size ++) 
			{
			cout << "nmeas : " << n << "\tSize : " << size << endl;
			if( size % 2 == 0)	start = ( N/2 - size / 2 );					//if size is even we go to the left of the center
			else start = ( N/2 - (size + 1) / 2 ) ;
		
			(*psi).position( start );
			MPO totalSx = build_totalSx( sites , start , size );
	
			// powers of the MPO: in ITensor v3, A*B = nmultMPO(A,prime(B)) with primes 2 -> 1
			Args args_mult = {"MaxDim",500,"Cutoff",1E-16};
			MPO Moment2 = nmultMPO(totalSx,prime(totalSx),args_mult); Moment2.mapPrime(2,1);
			MPO Moment3 = nmultMPO(Moment2,prime(totalSx),args_mult); Moment3.mapPrime(2,1);
			MPO Moment4 = nmultMPO(Moment2,prime(Moment2),args_mult); Moment4.mapPrime(2,1);
	
			complex<double> M1 = innerC( *psi , totalSx , *psi );
			complex<double> M2 = innerC( *psi , Moment2 , *psi );
			complex<double> M3 = innerC( *psi , Moment3 , *psi );
			complex<double> M4 = innerC( *psi , Moment4 , *psi );
		
			complex<double> C1 = M1;
			complex<double> C2 = M2 - M1*M1;
			complex<double> C3 = M3 - 3*M2*M1 + 2*M1*M1*M1;
			complex<double> C4 = M4 - 4*M3*M1 - 3*M2*M2 + 12*M2*M1*M1 - 6*M1*M1*M1*M1;
			
			vector<double> MSingleReal;
			vector<double> MSingleImag;
			
			MSingleReal.push_back( C1.real() );
			MSingleReal.push_back( C2.real() );
			MSingleReal.push_back( C3.real() );
			MSingleReal.push_back( C4.real() );
					
			MSingleImag.push_back( C1.imag() );
			MSingleImag.push_back( C2.imag() );
			MSingleImag.push_back( C3.imag() );
			MSingleImag.push_back( C4.imag() );
			
			MReal.push_back( MSingleReal );
			MImag.push_back( MSingleImag );
			}
					
			ofstream SaveFileReal( saveRealMomentsTwo.str().c_str() );
			ofstream SaveFileImag( saveImagMomentsTwo.str().c_str() ); 
	
			SaveFileReal << setprecision(10) << fixed;
			SaveFileImag << setprecision(10) << fixed;
	
			int col , row;
			
			for( row = 0 ; row < maxLength ; row++ )
				{
				for( col = 0 ; col <= 3 ; col++ )
					{
					SaveFileReal << MReal[row][col] << " ";
					SaveFileImag << MImag[row][col] << " ";
					}
				if( row != (maxLength-1) ) 
					{
					SaveFileReal << "\n";	
					SaveFileImag << "\n";				
					}
				}
				
			SaveFileReal.close();	
			SaveFileImag.close();	
		}
	
	else cout << "The files " << saveRealMomentsTwo.str() << " and "
			  << saveImagMomentsTwo.str()
			  << " already exist.\n" << endl;		
	
	}


// ======================================================================
// Mixed states (density matrix rho given as an MPO)
// ======================================================================

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
