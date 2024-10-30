#include "external_file.h"
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
			Meas = (*psi).A(start) * Obs * dag( prime( (*psi).A(start) , Site ) );
			}		
		else
			{
			Index ir = commonIndex( (*psi).A(start) , (*psi).A(start + 1) , Link);
			Meas = (*psi).A(start) * Obs * dag( prime( (*psi).A(start) , Site , ir ) );
			
			for( int row = 1 ; row < size ; row++ )
				{
				Sx = sites.op("Sx", start + row);
				Obs = expHermitian(Sx , theta * 1_i  ); 
				Meas *= (*psi).A(start + row) ;
				Meas *= Obs ;
				Meas *= dag( prime( (*psi).A(start + row) , Site , Link) );
				}
			
			Meas *= (*psi).A( start + size );
			Index il = commonIndex( (*psi).A( start + size ), (*psi).A( start + size -1 ), Link);
	
			Sx = sites.op("Sx", start + size );
			Obs = expHermitian(Sx , theta * 1_i  ); 
			Meas *= Obs;
			Meas *= dag( prime( (*psi).A( start + size ), il , Site) );
			}
			
		theta += theta_step( col , numberPoints);
	
		complex<double> SingleMeasure = Meas.cplx();
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

	MPO totalSx(ampo);	
	
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
	
			MPO Moment2;
			MPO Moment3;
			MPO Moment4;
		
			nmultMPO(totalSx,totalSx,Moment2,{"Maxm",500,"Cutoff",1E-16});
			nmultMPO(Moment2,totalSx,Moment3,{"Maxm",500,"Cutoff",1E-16});
			nmultMPO(Moment2,Moment2,Moment4,{"Maxm",500,"Cutoff",1E-16});
	
			complex<double> M1 = overlap( *psi , totalSx , *psi );
			complex<double> M2 = overlap( *psi , Moment2 , *psi );
			complex<double> M3 = overlap( *psi , Moment3 , *psi );
			complex<double> M4 = overlap( *psi , Moment4 , *psi );
		
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
	
