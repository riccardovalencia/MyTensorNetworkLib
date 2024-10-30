#include "TEBD.h"
#include <itensor/all.h>
#include <iostream>
#include <math.h>       

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------

//single time step for time evolution in bosonic quantum east model chain with losses
 
void
build_single_step_jumps( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const double gamma, const int j )
{

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id",j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Xj_plus_1  = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	*hterm = - 0.5 * Nj * ( J * Xj_plus_1 - U * Nj_plus_1) ;

	if( j==size-1) *hterm += 0.50 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   *hterm += 0.50 * Nj * Idj_plus_1 ;

	// *hterm  += -0.5 * gamma * Cplx_i * Nj * Idj_plus_1;
	// if(j==size -1) *hterm += - 0.5 * gamma * Cplx_i * Idj * Nj_plus_1; 

	// TESTING
//	*hterm  			   = - 0.5 * gamma * Cplx_i * Nj  * Idj_plus_1;
//	if(j==size -1) *hterm += - 0.5 * gamma * Cplx_i * Idj * Nj_plus_1;

	// density noise
	
	ITensor Nj_square = op(sites, "N" , j) * prime(op(sites, "N" , j),"Site");
	Nj_square.mapPrime(2,1);

	*hterm += -0.5 * gamma * Cplx_i * Nj_square * Idj_plus_1;

	Nj_square = op(sites, "N" , j+1) * prime(op(sites, "N" , j+1),"Site");
	Nj_square.mapPrime(2,1);

	if(j==size - 1) *hterm += -0.5 * gamma * Cplx_i * Idj * Nj_square; 
	
	}

// single time step for the time evoluton of the open bosonic quantum east model - the jump operators are in the vector<ITensor> Lj and Ljd

//single time step for time evolution in bosonic quantum east model chain with losses
 
void
build_single_step_jumps_v2( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const int j , vector<ITensor> &Lj, vector<ITensor> &Ljd)
{

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id",j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Xj_plus_1  = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	// TESTING H = 0

	*hterm = - 0.5 * Nj * ( J * Xj_plus_1 - U * Nj_plus_1) ;
	*hterm += 0.50 * Nj * Idj_plus_1 ;
	if( j==size-1) *hterm += 0.5  * Idj * Nj_plus_1;

	ITensor LdL_j = prime(Ljd[j-1],"Site") * Lj[j-1];
	LdL_j.mapPrime(2,1);

	// cerr << "Inside TEBD.cc site : " << j << endl;
	// PrintData(LdL_j);


	*hterm += -0.5 * Cplx_i * LdL_j * Idj_plus_1;

	if( j == size-1 )
	{
		LdL_j =  prime(Ljd[j],"Site") * Lj[j];
		LdL_j.mapPrime(2,1);
		*hterm += -0.5 * Cplx_i * Idj * LdL_j;
	}

	
	}



//single time step for time evolution in bosonic quantum east model chain

void
build_single_step( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const int j )
	{

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id",j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Sxj_plus_1 = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	*hterm = - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;
/*
	if( j==1 )     *hterm += 0.5  * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;
	if( j==size-1) *hterm += 0.25 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   *hterm += 0.25 * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;
*/

	if( j==size-1) *hterm += 0.50 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   *hterm += 0.50 * Nj * Idj_plus_1 ;


	}


//----------------------------------------------------------------------

//single time step for time evolution in bosonic quantum east model chain, with symmetry on n0

void
build_single_step_n0( ITensor *hterm , const SiteSet sites , const int size , const int n0, const double J , const double c , const int j )
	{

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id", j);
	ITensor Sxj = op(sites, "A" , j) + op(sites, "Adag" , j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Sxj_plus_1 = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	if( j== 1) 
	{	
		*hterm =  - 0.5 * n0 * ( J * Sxj * Idj_plus_1 - U * Nj * Idj_plus_1) ;
		*hterm += - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;
	}
	else *hterm = - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;

	if( j==1 )     *hterm += 0.5  * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;
	if( j==size-1) *hterm += 0.25 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   *hterm += 0.25 * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;

	}


// Build a full TEBD time step under the bosonic quantum east hamiltonian with only next-neighbour density-density interaction
// i.e. H = - 0.5 \sum_i n_i (exp(-s)\sigma_i+1^x - U n_i+1 - 1) where we do not fix any symmetry sector. 
// deprecated the usage of open function since 21.03.22 -> use 
void
build_TEBD_dt_step_H(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c , const string dynamics ,const double gamma)
{

	for(int j = 1; j <= size-1; j++)
		{
		ITensor hterm;
		if(dynamics == "closed")
			{
			// build_single_step_n0( &hterm , sites , size , n0, J , c , j );
			build_single_step( &hterm , sites , size , J , c , j );
			}
		if(dynamics == "open")
			{
			build_single_step_jumps( &hterm , sites , size , J , c , gamma, j );
			}
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm); 
		gates.push_back(g);
		}

	for(int j = size-1; j >= 1; j-=1)
		{
		ITensor hterm;
		if(dynamics == "closed")
			{
			// build_single_step_n0( &hterm , sites , size , n0, J , c , j );
				build_single_step( &hterm , sites , size , J , c , j );
			}
		if(dynamics == "open")
			{
			build_single_step_jumps( &hterm , sites , size , J , c , gamma, j );
			}
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm); 
		gates.push_back(g);
		}
}


// Perform adiabatic transformation from s=infty to s via a linear protocol


void
adiabatic_transformation_linear_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double beta)
{

	// beta controls the slope of the linear ramping. The greater is beta the slower is the protocol
	int L = length(*psi_start);
	double J_target = exp(-s);
	double t = 0;
	double J = 0;
	double tolerance = 1E-8;
	// auto TEBD_args = Args("Cutoff=",1E-16,"Verbose=",false,"MaxDim=",1000 );	
	auto TEBD_args = Args("Cutoff=",1E-10,"Verbose=",false,"MaxDim=",50 );	


	cerr << "dt : " << dt << endl;
	cerr << "beta : " << beta << endl;
	do
	{
		t += dt;
		J = J_target * t / beta;
		cerr << J << endl;
		auto gates = vector<BondGate>();
		build_TEBD_dt_step_H(gates, sites, L, dt, J, c);
		gateTEvol( gates , dt , dt , *psi_start , TEBD_args); 
	}while( J_target > J );




}


// Perform adiabatic transformation from s=infty to s via aa tanh(x) protocol


void
adiabatic_transformation_tanh_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double beta)
{

	// beta controls the slope of the tanh. The greater is beta the slower is the protocol
	int L = length(*psi_start);
	double J_target = exp(-s);
	double t = 0;
	double J = 0;
	double tolerance = 1E-6;
	auto TEBD_args = Args("Cutoff=",1E-16,"Verbose=",false,"MaxDim=",1000 );		
	
	do
	{
		t += dt;
		J = J_target * tanh(t/beta);
		cerr << J << endl;
		auto gates = vector<BondGate>();
		build_TEBD_dt_step_H(gates, sites, L, dt, J, c);
		gateTEvol( gates , dt , dt , *psi_start , TEBD_args); 
	}while((J_target - J)/(J_target + J) > tolerance );





	
	// int n/umber_steps = 5*int(T);

	// 

	// for(int step=1; step<=number_steps; step++)
	// {
	// 	double J = J_target * tanh(step/T);
	// 	cerr << J << endl;
	// 	auto gates = vector<BondGate>();
	// 	build_TEBD_dt_step_H(gates, sites, L, dt, J, c);
	// 	gateTEvol( gates , dt , dt , *psi_start , TEBD_args); 
	// }
	// *psi_start /= norm(*psi_start);


}


void
build_TEBD_dt_step_H_open(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c , vector<ITensor> &Lj, vector<ITensor> &Ljd )
{
	for(int j = 1; j <= size-1; j++)
		{
		ITensor hterm;
		build_single_step_jumps_v2( &hterm , sites , size , J , c , j , Lj, Ljd);
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm); 
		gates.push_back(g);
		}

	for(int j = size-1; j >= 1; j-=1)
		{
		ITensor hterm;
		build_single_step_jumps_v2( &hterm , sites , size , J , c , j ,Lj, Ljd);
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm); 
		gates.push_back(g);
		}
}
