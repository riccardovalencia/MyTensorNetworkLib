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
// local observables (sigma^x = a + a^dag)

static ITensor
make_sigma_x(const SiteSet& sites, const int j)
{
	return op(sites, "A", j) + op(sites, "Adag", j);
}


double 
measure_sigma_x_squared( MPS *state , const SiteSet sites , const int j )
{
	ITensor X = make_sigma_x(sites, j);
	return real(measure_local_operator(state, multSiteOps(X, X), j));
}


double 
measure_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	return real(measure_local_operator(state, make_sigma_x(sites, j), j));
}


double 
measure_sigma_x_n( MPS *state , const SiteSet sites , const int j )
{
	return real(measure_local_operator(state, multSiteOps(make_sigma_x(sites, j), op(sites, "N", j)), j));
}


double 
measure_n_sigma_x( MPS *state , const SiteSet sites , const int j )
{
	return real(measure_local_operator(state, multSiteOps(op(sites, "N", j), make_sigma_x(sites, j)), j));
}


void 
measure_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &occupation_number )
{
	for(int j = 1 ; j <= size ; j++)
		occupation_number.push_back( real(measure_local_operator(ground_state, op(sites, "N", j), j)) );
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


void 
measure_occupation_number_squared( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &square_occupation_number )
{
	for(int j = 1 ; j <= size ; j++)
	{
		ITensor Nj = op(sites, "N", j);
		square_occupation_number.push_back( real(measure_local_operator(ground_state, multSiteOps(Nj, Nj), j)) );
	}
}


void
measure_fock_probabilities( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space ,  vector<vector<double> > &projector_all_sites )
{
	for(int j = 1 ; j <= size ; j++)
	{
		Index s = sites(j);
		vector<double> probabilities;
		for( int n = 0 ; n <= cut_off_fock_space ; n++ )
		{
			ITensor projector = ITensor(s, prime(s));   // |n><n|
			projector.set(s(n+1), prime(s)(n+1), 1.);
			probabilities.push_back( real(measure_local_operator(ground_state, projector, j)) );
		}
		projector_all_sites.push_back(probabilities);
	}
}


// <O_i P_j>: two-point function for i != j, product of the operators on the same site otherwise
static Cplx
measure_pair(MPS *psi, const SiteSet& sites, const ITensor& O_i, const ITensor& P_j, const int i, const int j)
{
	if(i != j) return measure_two_point_function(psi, sites, O_i, P_j, i, j);
	return measure_local_operator(psi, multSiteOps(O_i, P_j), i);
}


void
measure_number_covariance( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN_system, vector<vector<double> > &covariance_matrix_NN_gaussian, vector<vector<double> > &relative_error)
{
	int size = length(*psi);

	(*psi).position(1);
	(*psi).normalize();

	for(int i = 1; i <= size; i++)
	{
		vector<double> covariance_ij, covariance_ij_gaussian, relative_error_ij;
		for(int j = 1; j <= size ; j++)
		{
			ITensor Ai  = op(sites, "A" , i),    Aj  = op(sites, "A" , j);
			ITensor Adi = op(sites, "Adag" , i), Adj = op(sites, "Adag" , j);
			ITensor Ni  = op(sites, "N" , i),    Nj  = op(sites, "N" , j);

			Cplx NiNj_c  = measure_pair(psi, sites, Ni,  Nj,  i, j);
			double NiNj  = (i != j) ? abs(NiNj_c) : real(NiNj_c);
			Cplx Aid_Ajd = measure_pair(psi, sites, Adi, Adj, i, j);
			Cplx Ai_Aj   = measure_pair(psi, sites, Ai,  Aj,  i, j);
			Cplx Aid_Aj  = measure_pair(psi, sites, Adi, Aj,  i, j);
			Cplx Ai_Ajd  = measure_pair(psi, sites, Ai,  Adj, i, j);

			double ni = real(measure_local_operator(psi, Ni, i));
			double nj = real(measure_local_operator(psi, Nj, j));
			Cplx   ai = measure_local_operator(psi, Ai, i);
			Cplx   aj = measure_local_operator(psi, Aj, j);

			// Wick factorization of <n_i n_j>_c for a Gaussian state
			Cplx NiNj_gaussian_connected = (Aid_Ajd * Ai_Aj + Aid_Aj * Ai_Ajd) - 2 * abs(ai)*abs(ai) * abs(aj)*abs(aj);
			double NiNj_connected = NiNj - ni*nj;

			covariance_ij.push_back(NiNj_connected);
			covariance_ij_gaussian.push_back(NiNj_gaussian_connected.real());
			relative_error_ij.push_back((NiNj_gaussian_connected.real()-NiNj_connected)/NiNj_connected);
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
	if( j > length(*psi) ) return -1;
	ITensor Adag = op(sites, "Adag", j);
	Cplx Nj_exp    = measure_local_operator(psi, op(sites, "N", j), j);
	Cplx adag2_exp = measure_local_operator(psi, multSiteOps(Adag, Adag), j);
	return 1 + 2 * Nj_exp.real() - 2 * abs(adag2_exp);
}


double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N)
{
	complex<double> N_measured = innerC((*psi),N,(*psi));
	complex<double> A2_measured = innerC((*psi),A,A,(*psi));

	double squeezing_dressed = 1 + 2 * N_measured.real() - 2 * abs(A2_measured);
	return squeezing_dressed;
}
