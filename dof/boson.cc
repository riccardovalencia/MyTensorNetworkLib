/**
 * @file boson.cc
 * @brief Implementation of boson.h (the functions are documented in the header).
 */
#include "boson.h"
#include "../mps/mps_tools.h"
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


// ----------------------------------------------------------
// amplitudes of single-site states in the Fock basis |0>, ..., |dim-1>

static vector<Cplx>
fock_amplitudes( int dim, int n )
{
    vector<Cplx> a(dim, 0.);
    a[n] = 1.;
    return a;
}

static vector<Cplx>
coherent_amplitudes( int dim, Cplx alpha )
{
    vector<Cplx> a;
    for(int k = 0 ; k < dim ; k++) a.push_back(coherent_state_amplitude(alpha, k));
    return a;
}

static vector<Cplx>
squeezed_amplitudes( int dim, double r )
{
    vector<Cplx> a;
    for(int k = 0 ; k < dim ; k++) a.push_back(k % 2 == 0 ? squeezed_state_amplitude(r, k) : 0.);
    return a;
}


// ----------------------------------------------------------
// product states

void
set_site_occupation( MPS* psi, const SiteSet sites, int site, int n )
{
    set_site_tensor(psi, sites, site, fock_amplitudes(dim(sites(site)), n));
}

void
set_vacuum_state( MPS* psi, const SiteSet sites )
{
    for(int j = 1 ; j <= length(*psi) ; j++) set_site_occupation(psi, sites, j, 0);
}

void
set_unit_filling_state( MPS* psi, const SiteSet sites )
{
    for(int j = 1 ; j <= length(*psi) ; j++) set_site_occupation(psi, sites, j, 1);
}

void
set_fock_excitation( MPS* psi, const SiteSet sites, int n0, int site )
{
    set_vacuum_state(psi, sites);
    set_site_occupation(psi, sites, site, n0);
}

void
set_kink_state( MPS* psi, const SiteSet sites, int number_ones )
{
    for(int j = 1 ; j <= length(*psi) ; j++) set_site_occupation(psi, sites, j, j <= number_ones ? 1 : 0);
}


// ----------------------------------------------------------
// coherent, squeezed and cat states

void
set_coherent_state_on_site( MPS* psi, const SiteSet sites, int site, Cplx alpha )
{
    set_vacuum_state(psi, sites);
    set_site_tensor(psi, sites, site, coherent_amplitudes(dim(sites(site)), alpha));
}

void
set_coherent_state_all_sites( MPS* psi, const SiteSet sites, Cplx alpha )
{
    for(int j = 1 ; j <= length(*psi) ; j++) set_site_tensor(psi, sites, j, coherent_amplitudes(dim(sites(j)), alpha));
}

void
set_squeezed_state_on_site( MPS* psi, const SiteSet sites, int site, double r )
{
    set_vacuum_state(psi, sites);
    set_site_tensor(psi, sites, site, squeezed_amplitudes(dim(sites(site)), r));
}

void
set_cat_state_on_site( MPS* psi, const SiteSet sites, int site, Cplx alpha )
{
    // (|alpha> + |-alpha>) normalized
    MPS plus  = randomMPS(sites);
    MPS minus = randomMPS(sites);
    set_coherent_state_on_site(&plus,  sites, site,  alpha);
    set_coherent_state_on_site(&minus, sites, site, -alpha);
    *psi = sqrt(0.5) * sum(plus, minus);
    (*psi).position(1);
    (*psi).normalize();
}


// ----------------------------------------------------------------------
// Initial bosonic state of the form |n0> |0000...0>.




// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.


// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.


// ----------------------------------------------------------------------
// Initial all one bosonic state |1111...1>.


// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.




// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.








// factorial n!

double
factorial(int n)
{
	double factorial = 1;
	for(int j=1; j<=n ; j++) factorial *= j;
	return factorial;	
}


// weight coherent state

complex<double> 
coherent_state_amplitude( const complex<double> alpha, const int k)
{
	complex<double> alpha_power_k = pow(alpha,k);
	return exp(-abs(alpha)*abs(alpha) / 2 ) * (alpha_power_k.real() + 1i*alpha_power_k.imag()) / sqrt(factorial(k));
}


// weight squeezed state

double 
squeezed_state_amplitude( const double r, const int k)
{
	double result;
	result  = pow((-1 * tanh(r)),int(k/2));
	result *= sqrt(factorial(k));
	result /= (pow(2,k/2) * factorial(int(k/2)));
	result /= sqrt(cosh(r));
	return result;
}


//----------------------------------------------------------------------
// expectation value: <\sigma_j^x^2>
double 
measure_sigma_x_squared( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);
	ITensor observable2 = prime(observable);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}


//----------------------------------------------------------------------
// expectation value: <\sigma_j^x>
double 
measure_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(ket,"Site"));
	double expectation_value = elt(bra * observable * ket);

	return expectation_value;

}


//----------------------------------------------------------------------
// expectation value: <\sigma_j^x n_j>
double 
measure_sigma_x_n( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites, "A", j);
	observable += op(sites,"Adag",j);
	ITensor observable2 = prime(op(sites,"N",j));

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}


//----------------------------------------------------------------------
// expectation value: <n_j \sigma_j^x>
double 
measure_n_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	ITensor observable = op(sites,"N",j);
	ITensor observable2 = op(sites, "A", j);
	observable2 += op(sites,"Adag",j);
	observable2 = prime(observable2);

	(*state).position(j);
	
	auto ket = (*state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);

	return expectation_value;

}


// measure occupation number along the 1D chain

void 
measure_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &occupation_number )
{
for(int j = 1 ; j <= size ; j++)
	{
	ITensor observable = op(sites, "N" , j);
	(*ground_state).position(j);
	auto ket = (*ground_state)(j);
	auto bra = dag(prime(ket,"Site"));
	complex<double> expectation_value = eltC(bra * observable * ket);
	occupation_number.push_back( expectation_value.real() );
	}
}


double
compute_imbalance( vector<double> &occupation_number, int k)
{

	double nk = occupation_number[k];
	occupation_number.erase(occupation_number.begin()+k);
	double nmax = *max_element(occupation_number.begin(), occupation_number.end());
	
	return (nk-nmax)/(nk+nmax);

}


double compute_max_cutoff_probability(vector<vector<double> > &projector_all_sites,const int size,const int cut_off)
{
	double maximum = -1;
	for(int j = 1 ; j <= size ; j++) maximum = std::max(projector_all_sites[j-1][cut_off-1], maximum);

	return maximum;

}


// measure squareoccupation number along the 1D chain

void 
measure_occupation_number_squared( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &square_occupation_number )
{
for(int j = 1 ; j <= size ; j++)
	{
	ITensor observable = op(sites, "N" , j);
	ITensor observable2 = prime(op(sites, "N" , j));

	(*ground_state).position(j);
	auto ket = (*ground_state)(j);
	auto bra = dag(prime(prime(ket,"Site"),"Site"));
	double expectation_value = elt(bra * observable * observable2 * ket);
	square_occupation_number.push_back( expectation_value );
	}
	
}


//---------------------------------------------------------------------

// measure of the projector along all the sites and all the Fock space

void measure_fock_probabilities( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space ,  vector<vector<double> > &projector_all_sites ,  vector<double> &occupation_number)
{
	for(int j = 1 ; j <= size ; j++)
	{
	vector<double> projector_single_site;
	(*ground_state).position(j);


	double check_number_particles = 0.;
	double check_normalization_projector = 0. ;

	for( int n = 0 ; n <= cut_off_fock_space ; n++ )
		{
		Index physical_index = sites(j);
		Index physical_index_prime = prime(sites(j));
		ITensor projector = ITensor(physical_index, physical_index_prime); //initialized with all the elements equal to zero.
	
		projector.set(physical_index(n+1),physical_index_prime(n+1),1.);

		ITensor ket = (*ground_state)(j);
		ITensor bra = dag(prime((*ground_state)(j),"Site"));
		
		complex<double> expectation_value_c = eltC(bra * projector * ket);
		double expectation_value = expectation_value_c.real();
		projector_single_site.push_back(expectation_value);

		check_number_particles += n * expectation_value ;
		check_normalization_projector += expectation_value;
	
		}

	cerr << "Site : " << j << " Difference (occupation_number - projector) " << (occupation_number[j-1]-check_number_particles) << endl;	
	cerr << "Site : " << j << " Normalization projector (sum of projectors should be 1) : " << check_normalization_projector << endl;	
	projector_all_sites.push_back(projector_single_site);
	}
}


// Measure the covariance matrix size X size. Each element is <N_i N_j>_c - 
void
measure_number_covariance( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN_system, vector<vector<double> > &covariance_matrix_NN_gaussian, vector<vector<double> > &relative_error)
{
	int size = length(*psi);

	(*psi).position(1);
	(*psi).normalize();
	
	for(int i = 1; i <= size; i++)
	{
		vector<double> covariance_ij;
		vector<double> covariance_ij_gaussian;
		vector<double> relative_error_ij;

		for(int j = 1; j <= size ; j++)
		{
			
			double NiNj;
			complex<double> Aid_Ajd ;
			complex<double> Aid_Aj  ;
			complex<double> Ai_Aj   ; 
			complex<double> Ai_Ajd  ;

			// define A 
			ITensor Ai_op  = op(sites, "A" , i);
			ITensor Aj_op  = op(sites, "A" , j);
			ITensor Adi_op  = op(sites, "Adag" , i);
			ITensor Adj_op  = op(sites, "Adag" , j);

			(*psi).position(i);
			auto ket = (*psi)(i);
			auto bra = dag(prime(ket,"Site"));
			complex<double> Ai = eltC( bra * Ai_op * ket);

			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			complex<double> Aj = eltC( bra * Aj_op * ket);

			// attempt - 14.12.21
			ITensor Ni_op  = op(sites, "N" , i);
			ITensor Nj_op  = op(sites, "N" , j);


			if( i != j)
			{
				NiNj    = abs(measure_two_point_function(psi,sites ,Ni_op,Nj_op , i , j));
				Aid_Ajd = measure_two_point_function(psi,sites,Adi_op,Adj_op , i , j);
				Ai_Aj   = measure_two_point_function(psi,sites,Ai_op,Aj_op , i , j);
				Aid_Aj  = measure_two_point_function(psi,sites,Adi_op,Aj_op  , i , j);
				Ai_Ajd  = measure_two_point_function(psi,sites,Ai_op,Adj_op  , i , j);	
			}


			else
			{
				(*psi).position(i);
				auto ket = (*psi)(i);
				auto bra =  dag(prime(prime(ket,"Site"),"Site"));
				NiNj = elt(bra * prime(Ni_op) * Ni_op * ket);
				Aid_Ajd = eltC(bra * prime(Adi_op) * Adi_op * ket);
				Ai_Aj   = eltC(bra * prime(Ai_op) * Ai_op * ket);
				Aid_Aj = eltC(bra * prime(Adi_op) * Ai_op * ket);
				Ai_Ajd = eltC(bra * prime(Ai_op) * Adi_op * ket);

			}
		
			(*psi).position(i);
			ket = (*psi)(i);
			bra = dag(prime(ket,"Site"));
			double Ni = elt( bra * Ni_op * ket);
			Ai = eltC( bra * Ai_op * ket);
 
			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			double Nj = elt( bra * Nj_op * ket);
			Aj = eltC( bra * Aj_op * ket);

			if( i == j && i==1)
			{
			cerr << " Observables " << endl;
			cerr << "NN : " << NiNj << endl;
			cerr << "N : " << Ni << endl;
			cerr << "A : " << Ai << endl; 
			cerr << "AdA : " << Aid_Aj << endl;
			cerr << "AdAd : " << Aid_Ajd << endl;
			}


			// measure average occupation number
			ITensor observable = op(sites, "N" , i);
			(*psi).position(i);
			ket = (*psi)(i);
			bra = dag(prime(ket,"Site"));
			double ni = eltC(bra * observable * ket).real();

			observable = op(sites, "N" , j);
			(*psi).position(j);
			ket = (*psi)(j);
			bra = dag(prime(ket,"Site"));
			double nj = eltC(bra * observable * ket).real();

			if( i == j && i==1)
			{
			cerr << "N_" << i << " " << ni << " " << Ni << endl;


			}


			complex<double> NiNj_gaussian_approx_connected = (Aid_Ajd * Ai_Aj + Aid_Aj * Ai_Ajd) - 2 * abs(Ai)*abs(Ai) * abs(Aj)*abs(Aj);
			double NiNj_connected = (NiNj - Ni*Nj);
		
			if( i == j && i==1)
			{
			cerr << "N_" << i << "N_" << j << " : " << NiNj_gaussian_approx_connected << " " << NiNj_connected << " " << (abs(NiNj_gaussian_approx_connected)-NiNj_connected)/abs(NiNj_gaussian_approx_connected) << endl;
			}
			covariance_ij.push_back(NiNj_connected);
			covariance_ij_gaussian.push_back(NiNj_gaussian_approx_connected.real());
			relative_error_ij.push_back((NiNj_gaussian_approx_connected.real()-NiNj_connected)/NiNj_connected);

		}
		covariance_matrix_NN_system.push_back(covariance_ij);
		covariance_matrix_NN_gaussian.push_back(covariance_ij_gaussian);
		relative_error.push_back(relative_error_ij);
		
	}
}


double
measure_variance_x(MPS *psi, MPO *A, MPO *Adag )
{

	complex<double> a_expectation_value     = innerC(*psi, *A   , *psi);
	complex<double> aa_expectation_value    = innerC(*psi, *A   , *A, *psi);
	complex<double> adaga_expectation_value = innerC(*psi, *Adag, *A, *psi);

	complex<double> delta_x = 1 + 2 * adaga_expectation_value + 2 * aa_expectation_value.real() - 4 * a_expectation_value.real() * a_expectation_value.real();

	return delta_x.real();
}


double
measure_variance_p(MPS *psi, MPO *A, MPO *Adag )
{

	complex<double> a_expectation_value     = innerC(*psi, *A   , *psi);
	complex<double> aa_expectation_value    = innerC(*psi, *A   , *A, *psi);
	complex<double> adaga_expectation_value = innerC(*psi, *Adag, *A, *psi);

	complex<double> delta_p = 1 + 2 * adaga_expectation_value - 2 * aa_expectation_value.real() - 4 * a_expectation_value.imag() * a_expectation_value.imag();

	return delta_p.real();

}


double
measure_squeezing(MPS *psi, const SiteSet sites, const int j)
{
	if( j <= length(*psi) )
	{
		ITensor Nj = op(sites, "N" , j);
		ITensor adag2 = op(sites, "Adag", j) * prime(op(sites, "Adag",j));
		adag2.mapPrime(2,1);


		(*psi).position(j);
		auto ket = (*psi)(j);
		auto bra = dag(prime(ket,"Site"));
		complex<double> Nj_exp = eltC(bra * Nj * ket);
		complex<double> adag2_exp = eltC(bra * adag2 * ket);
		double squeezing = 	1 + 2 * Nj_exp.real() - 2 * abs(adag2_exp);
		return squeezing;
	}
	else return -1;
}


double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N)
{
	complex<double> N_measured = innerC((*psi),N,(*psi));
	complex<double> A2_measured = innerC((*psi),A,A,(*psi));

	double squeezing_dressed = 1 + 2 * N_measured.real() - 2 * abs(A2_measured);
	return squeezing_dressed;
}
