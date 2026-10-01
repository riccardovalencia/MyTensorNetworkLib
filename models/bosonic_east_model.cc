/**
 * @file bosonic_east_model.cc
 * @brief Implementation of bosonic_east_model.h (the functions are documented in the header).
 */
#include "bosonic_east_model.h"
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
// bond term on (j, j+1) of the bosonic quantum east model chain

ITensor
make_bosonic_east_model_bond_hamiltonian( const SiteSet sites , const int size , const double J , const double c , const int j )
	{
    ITensor hterm;

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id",j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Sxj_plus_1 = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	hterm = - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;

	if( j==size-1) hterm += 0.50 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   hterm += 0.50 * Nj * Idj_plus_1 ;
    return hterm;
}


// n^2 on site j

static ITensor
make_number_squared(const SiteSet& sites, const int j)
{
	ITensor N2 = op(sites, "N" , j) * prime(op(sites, "N" , j),"Site");
	N2.mapPrime(2,1);
	return N2;
}


// bond term plus the non-hermitian density-noise term -i gamma/2 n_j^2 (and n_{j+1}^2 on the last bond)

ITensor
make_bosonic_east_model_bond_hamiltonian_dephasing( const SiteSet sites , const int size , const double J , const double c , const double gamma, const int j )
{
	ITensor hterm = make_bosonic_east_model_bond_hamiltonian(sites, size, J, c, j);
	hterm += -0.5 * gamma * Cplx_i * make_number_squared(sites, j) * op(sites, "Id", j+1);
	if(j==size - 1) hterm += -0.5 * gamma * Cplx_i * op(sites, "Id", j) * make_number_squared(sites, j+1);
    return hterm;
}


//----------------------------------------------------------------------

//single time step for time evolution in bosonic quantum east model chain, with symmetry on n0

ITensor
make_bosonic_east_model_bond_hamiltonian_n0( const SiteSet sites , const int size , const int n0, const double J , const double c , const int j )
	{
    ITensor hterm;

	double U = 1 - 2*c;		
	ITensor Nj  = op(sites, "N" , j);					
	ITensor Idj = op(sites, "Id", j);
	ITensor Sxj = op(sites, "A" , j) + op(sites, "Adag" , j);

	ITensor Nj_plus_1  = op(sites, "N" , j+1);		
	ITensor Sxj_plus_1 = op(sites, "A" , j+1) + op(sites, "Adag" , j+1);
	ITensor Idj_plus_1 = op(sites, "Id" , j+1);

	if( j== 1) 
	{	
		hterm =  - 0.5 * n0 * ( J * Sxj * Idj_plus_1 - U * Nj * Idj_plus_1) ;
		hterm += - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;
	}
	else hterm = - 0.5 * Nj * ( J * Sxj_plus_1 - U * Nj_plus_1) ;

	if( j==1 )     hterm += 0.5  * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;
	if( j==size-1) hterm += 0.25 * Nj * Idj_plus_1 + 0.5  * Idj * Nj_plus_1;
	else 		   hterm += 0.25 * Nj * Idj_plus_1 + 0.25 * Idj * Nj_plus_1;
    return hterm;
}


// Full TEBD time step under the bosonic quantum east Hamiltonian with only next-neighbour density-density
// interaction, H = - 0.5 \sum_i n_i (exp(-s)\sigma_i+1^x - U n_i+1 - 1), without fixing a symmetry sector.

vector<TebdGate>
make_bosonic_east_model_gates(const SiteSet sites, const int size, const double dt, const double J, const double c , const string dynamics ,const double gamma)
{
    if(dynamics != "closed" && dynamics != "open")
        throw ITError("make_bosonic_east_model_gates: dynamics must be \"closed\" or \"open\"");

    vector<TebdGate> gates;
	for(int j = 1; j <= size-1; j++)
	{
		ITensor hterm = (dynamics == "closed") ? make_bosonic_east_model_bond_hamiltonian( sites , size , J , c , j )
		                                       : make_bosonic_east_model_bond_hamiltonian_dephasing( sites , size , J , c , gamma, j );
		gates.push_back(TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hterm).gate()));
	}
    return make_symmetric_sweep(gates);
}


//----------------------------------------------------------------------
// Terms of the bosonic quantum east model, multiplied by prefactor, on the sites first..size:
//   H = - 0.5 n0 (J sigma^x_first - U n_first - 1)                      (only if n0 != 0)
//       - 0.5 sum_{j=first}^{size-1} n_j (J sigma^x_{j+1} - U n_{j+1} - 1)
//       + 0.5 (1 - symmetry) n_size
// with sigma^x = a + a^dag and U = 1 - 2c. The virtual site 0 with occupation n0 couples to site first.

static AutoMPO
make_bosonic_east_model_terms( const SiteSet& sites, int size, int first, int n0, double symmetry, double J, double c, double prefactor = 1. )
{
	double U  = 1 - 2*c;
	auto ampo = AutoMPO(sites);

	if(n0 != 0)
		{
		ampo += - prefactor * n0 * J * 0.5 , "A"    , first;
		ampo += - prefactor * n0 * J * 0.5 , "Adag" , first;
		ampo +=   prefactor * 0.5 * n0 * U , "N"    , first;
		ampo +=   prefactor * n0 * 0.5     , "Id"   , first;
		}

	for(int j = first ; j <= size-1 ; j++)
		{
		ampo += - prefactor * 0.5 * J , "N" , j , "A"    , j+1;
		ampo += - prefactor * 0.5 * J , "N" , j , "Adag" , j+1;
		ampo +=   prefactor * 0.5 * U , "N" , j , "N"    , j+1;
		ampo +=   prefactor * 0.5     , "N" , j ;
		}

	ampo += - prefactor * 0.5 * symmetry , "N", size;
	ampo +=   prefactor * 0.5            , "N", size;
	return ampo;
}


MPO
make_bosonic_east_model_mpo( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{
	return toMPO(make_bosonic_east_model_terms(sites, size, 1, n0, symmetry, exp(-s), c), {"Exact=",true});
}


MPO
make_bosonic_east_model_mpo_with_drift( const SiteSet sites, int size , int n0, double symmetry , double s, double c, double Omega)
{
	auto ampo = make_bosonic_east_model_terms(sites, size, 1, n0, symmetry, exp(-s), c);
	for(int j = 2 ; j <= size-1 ; j++)
		{
		ampo += Omega , "A"    , j;
		ampo += Omega , "Adag" , j;
		}
	return toMPO(ampo, {"Exact=",true});
}


MPO
make_bosonic_east_model_mpo_minus( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{
	return toMPO(make_bosonic_east_model_terms(sites, size, 1, n0, symmetry, exp(-s), c, -1.), {"Exact=",true});
}


MPO
make_bosonic_east_model_mpo_onsite_hopping( const SiteSet sites, int size , double s, double c, double epsilon, double t)
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
make_bosonic_east_model_mpo_onsite( const SiteSet sites, int size , double epsilon)
{
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
make_bosonic_east_model_mpo_onsite_nonext( const SiteSet sites, int size , int n0, double symmetry , double c )
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
// Site 1 is not touched: it plays the role of site 0, to find the ground state within the symmetry
// sector fixed by n0 from a state of the form |n0>|psi>.

MPO
make_bosonic_east_model_mpo_n0_untouched( const SiteSet sites, int size , int n0, double symmetry , double s, double c)
{
	return toMPO(make_bosonic_east_model_terms(sites, size, 2, n0, symmetry, exp(-s), c), {"Exact=",true});
}


// ------------------------------------------------------------
// No virtual site 0: n_1 is not a conserved quantity, and DMRG is not expected to conserve it (when
// initialized in |n_0> (x) |psi> it may stay stuck for a long time, since |n_0> is an eigenstate of the
// operators acting on site 1).

MPO
make_bosonic_east_model_mpo_n0_not_fixed( const SiteSet sites, int size , double symmetry , double s, double c)
{
	return toMPO(make_bosonic_east_model_terms(sites, size, 1, 0, symmetry, exp(-s), c), {"Exact=",true});
}


// exp(i dt H) of make_bosonic_east_model_mpo_n0_not_fixed with hopping J, to evolve operators

MPO
make_bosonic_east_model_evolution_mpo( const SiteSet sites, int size , double symmetry , double J, double c, double dt)
{
	return toExpH(make_bosonic_east_model_terms(sites, size, 1, 0, symmetry, J, c), Cplx_i * dt);
}


// eigenvalue of the symmetry sector for the given Fock-space cutoff: entry number cut_off_fock_space of
// "<symmetry_sector_dir>/symmetry_sector<sector>_maxcutoff30_s<s>_c<c>.dat"

static double
read_symmetry_eigenvalue(const string& symmetry_sector_dir, int symmetry_sector, int cut_off_fock_space, double s, double c)
{
	string file_name = tinyformat::format("%s/symmetry_sector%d_maxcutoff30_s%.2f_c%.2f.dat",symmetry_sector_dir,symmetry_sector,s,c);
	ifstream file(file_name);
	if(!file.is_open()) throw ITError("compute_bosonic_east_model_energy_variance: cannot open " + file_name);

	float skipped;
	for(int i = 1 ; i < cut_off_fock_space ; i++) file >> skipped;
	double symmetry;
	file >> symmetry;
	return symmetry;
}


double 
compute_bosonic_east_model_energy_variance(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c, const string symmetry_sector_dir)
{
	double symmetry = read_symmetry_eigenvalue(symmetry_sector_dir, symmetry_sector, cut_off_fock_space, s, c);
	MPO H = make_bosonic_east_model_mpo( sites, size , n0, symmetry , s, c);
	double energy = inner(*psi,H,*psi);
	return inner(*psi,H,H,*psi) - energy * energy;
}
