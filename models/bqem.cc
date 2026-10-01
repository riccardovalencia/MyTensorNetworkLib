/**
 * @file bqem.cc
 * @brief Implementation of bqem.h (the functions are documented in the header).
 */
#include "bqem.h"
#include "../models/spin_chain.h"
#include <itensor/all.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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
build_tebd_dt_step(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c , const string dynamics ,const double gamma)
{

	for(int j = 1; j <= size-1; j++)
		{
		ITensor hterm;
		if(dynamics == "closed")
			{
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


void
build_tebd_dt_step_open(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c , vector<ITensor> &Lj, vector<ITensor> &Ljd )
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


//----------------------------------------------------------------------
// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector + 0.5 * n_L

MPO
hamiltonian_bqem( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	double U = 1-2*c;

	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo += 0.5 * n0 * U , "N", 1 ;
	ampo += n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;
		}

	double eigenvalue_symmetry = -1*symmetry;


	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size; 

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


//----------------------------------------------------------------------
// Hamiltonian H = - 0.5 n_0 ( exp(-s) \sigma_1^x -1) - 0.5 \sum_{j=1}^{L-1} n_j (exp(-s) \sigma_{j+1}^x - (1-2c) n_{j+1} -1) + 0.5 n_L * sector + 0.5 * n_L

MPO
hamiltonian_bqem_with_drift( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{

	double U = 1-2*c;
	double Omega = 0.05;
	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo += 0.5 * n0 * U , "N", 1 ;
	ampo += n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;
		}

	for(int j = 2 ; j <= size-1 ; j++)
		{
		ampo += Omega , "A" , j ;
		ampo += Omega , "Adag" , j;
		}

	double eigenvalue_symmetry = -1*symmetry;


	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size; 

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


MPO
hamiltonian_bqem_minus( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{


	double U = 1-2*c;

	auto ampo_minus = AutoMPO(sites);

	ampo_minus += n0 * exp(-s) * 0.5  , "A" , 1 ;
	ampo_minus += n0 * exp(-s) * 0.5  , "Adag" , 1 ;
	ampo_minus += -0.5 * n0 * U , "N", 1 ;
	ampo_minus += -n0 * 0.5 , "Id", 1;

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo_minus += exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo_minus += exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo_minus += -0.5 * U , "N", j , "N" , j+1 ;
		ampo_minus += -0.5 , "N", j , "Id", j+1;
		}


	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = symmetry;

	ampo_minus += -0.5 * eigenvalue_symmetry , "N", size; 
	ampo_minus += -0.5, "N", size;

	MPO H_minus = toMPO(ampo_minus,{"Exact=",true});

	return H_minus;

}


MPO
hamiltonian_bqem_onsite_hopping( const SiteSet sites, int size , double s, double c, double epsilon, double t)
{


	double U = 1-2*c;
	auto ampo = AutoMPO(sites);


	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - exp(-s) * 0.5 , "N" , j , "A" , j+1;
		ampo += - exp(-s) * 0.5 , "N" , j , "Adag" , j+1;
		ampo += 0.5 * U , "N", j , "N" , j+1 ;
		ampo += 0.5 , "N", j , "Id", j+1;

		// on site density-density
		ampo += 0.5 * epsilon , "N", j , "N", j;

		ampo += -t * 0.5, "Adag", j , "A", j+1;
		ampo += -t * 0.5, "A", j+1 , "Adag", j;

		}

	// on site density-density
	ampo += 0.5 * epsilon , "N", size , "N", size;

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


MPO
hamiltonian_bqem_onsite( const SiteSet sites, int size , int n0, double symmetry , double s, double c , double epsilon)
{
	 
	double U = 1-2*c;
	auto ampo = AutoMPO(sites);


	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += 0.5 * epsilon, "N", j , "N", j ;
		}

	ampo += 0.5 * epsilon , "N", size , "N", size;

	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}




MPO
hamiltonian_bqem_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double s, double c )
{

	auto ampo = AutoMPO(sites);

	ampo +=   n0 * 0.5 , "Id", 1;

	// non serve mettere n_0^2 dato che e' una costante.

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += 0.5 , "N", j , "Id", j+1;

		// on site density-density
		ampo += 0.5 * c , "N", j , "N", j;

		}

	// on site density-density
	ampo += 0.5 * c , "N", size , "N", size;
	ampo += - symmetry * 0.5 , "N", size; 
	ampo += 0.5, "N", size;

	MPO H = toMPO(ampo,{"Exact=",true});


	return H;
}


// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// site 0 is not touched. This is used to find the ground state within a certain symmetry sector fixed by
// n0 from a site of the form |n0>|psi>.

MPO
hamiltonian_bqem_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{


	double U  = 1 - 2 * c;
	auto ampo = AutoMPO(sites);

	ampo += - n0 * exp(-s) * 0.5  , "A" , 2 ;
	ampo += - n0 * exp(-s) * 0.5  , "Adag" , 2 ;
	ampo +=  0.5 * n0 * U , "N", 2 ;
	ampo +=   n0 * 0.5 , "Id", 2;

	for(int j = 2 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * exp(-s) , "N" , j , "A"    , j+1;
		ampo += - 0.5 * exp(-s) , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}


	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;
 
	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;


	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


// ------------------------------------------------------------
// Build Hamiltonian of Bosonic quantum east model with next-neighbour density-density interaction.
// No on-site density-density interaction
// Notice that n0 is not a number, but its defined in the Hamiltonian itself. So, it is not a conserved quantity
// and we do not expect DMRG to conserve it (it will try to minimize the Hamiltonian and that's it. At most it 
// will be stuck for a long time if you initialize DMRG with a state |n_0> \otimes |\psi> since |n_0> is an eigenstate
// of the operators acting on the 0-th site). 

MPO
hamiltonian_bqem_n0_not_fixed( const SiteSet sites, int size , double symmetry , double s, double c)
{


	double U = 1-2*c;
	auto ampo = AutoMPO(sites);

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * exp(-s) , "N" , j , "A"    , j+1;
		ampo += - 0.5 * exp(-s) , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}


	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;

	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;


	MPO H = toMPO(ampo,{"Exact=",true});

	return H;
}


// Exponential of the Hamiltonian. This is done in order to apply it to operators

MPO
exp_hamiltonian_bqem_n0_not_fixed( const SiteSet sites, int size , int n0, double symmetry , double J, double c, double dt)
{
	 
	double U  = 1 - 2 * c;
	auto ampo = AutoMPO(sites);

	for(int j = 1 ; j <= size-1 ; j++)
		{
		ampo += - 0.5 * J , "N" , j , "A"    , j+1;
		ampo += - 0.5 * J , "N" , j , "Adag" , j+1;
		ampo +=   0.5 * U , "N", j , "N" , j+1 ;
		ampo +=   0.5     , "N", j ;
		}

	// eigenvalues for the finite cut-off system
	
	double eigenvalue_symmetry = -1*symmetry;

	ampo += 0.5 * eigenvalue_symmetry , "N", size; 
	ampo += 0.5, "N", size;	

	MPO expH = toExpH( ampo , Cplx_i * dt );
	return expH;
}


double 
compute_variance_hamiltonian_bqem(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c, const string symmetry_sector_dir)
{
	double symmetry;
	ifstream symmetry_sector_file;
	string name_symmetry_sector_file = tinyformat::format("%s/symmetry_sector%d_maxcutoff30_s%.2f_c%.2f.dat",symmetry_sector_dir,int(symmetry_sector),s,c);
	cerr << name_symmetry_sector_file << endl;
	symmetry_sector_file.open(name_symmetry_sector_file);
	if(!symmetry_sector_file.is_open())
		throw ITError("compute_variance_hamiltonian_bqem: cannot open " + name_symmetry_sector_file);

	for (int i = 1; i < cut_off_fock_space; i++)
	{
		float tmp;
		symmetry_sector_file >> tmp;
	}

	symmetry_sector_file >> symmetry;

	MPO H = hamiltonian_bqem( sites, size , n0, symmetry , s, c) ; 


	double variance = inner((*psi),H,H,(*psi))  -   inner((*psi),H,(*psi)) *  inner((*psi),H,(*psi));
	cerr << "I have computed variance " << endl;
	return variance;
}
