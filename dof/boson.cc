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


// ----------------------------------------------------------------------
// Initial bosonic state of the form |n0> |0000...0>.
void
initial_state_n0_excitation( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position)
{
	
	// SITE 1

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(excitation_position == 1)
	{
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	else
	{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	(*psi).set(1,wf);

	// SITE 2 TO SIZE-1

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);

		if(j==excitation_position)
		{
			for( int d=1; d <= n0; d++) wf.set(sj(d),li(1),ri(1), 0);
			wf.set(sj(n0+1),li(1),ri(1),1);
			for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		}
		else
		{
			wf.set(sj(1),li(1),ri(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		}
		(*psi).set(j,wf);
	}

	// SITE SIZE (LAST ONE)
	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);

	if(excitation_position == size)
	{
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	else
	{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}
	(*psi).set(size,wf);

}


void
initial_state_n0_excitation_pinned( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position)
{
	
	// SITE 1

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(excitation_position == 1)
	{
		sj = sites(1);
		li = commonIndex((*psi)(1),(*psi)(2));
		wf = ITensor(sj,li);
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		(*psi).set(1,wf);
	}


	// SITE 2 TO SIZE-1

	else if(excitation_position > 1 && excitation_position < size)
	{
		sj = sites(excitation_position);
		li = commonIndex((*psi)(excitation_position-1),(*psi)(excitation_position));
		Index ri = commonIndex((*psi)(excitation_position),(*psi)(excitation_position+1));
		wf = ITensor(sj,li,ri);


		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1),ri(1), 0);
		wf.set(sj(n0+1),li(1),ri(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		
		(*psi).set(excitation_position,wf);
	}

	// SITE SIZE (LAST ONE)
	
	if(excitation_position == size)
	{
		sj  = sites(size);
		li = commonIndex((*psi)(size-1),(*psi)(size));
		wf = ITensor(sj,li);
		for( int d=1; d <= n0; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n0+1),li(1),1);
		for( int d=n0+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		(*psi).set(size,wf);
	}


}


// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.


// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.
void
initial_state_vacuum_state_correct_link( MPS* psi, const SiteSet sites, const int size )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 1);
	for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		wf.set(sj(1),li(1),ri(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 1);
	for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(size,wf);
}


// ----------------------------------------------------------------------
// Initial all one bosonic state |1111...1>.
void
initial_state_all_one_state_correct_link( MPS* psi, const SiteSet sites, const int size )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 0);
	wf.set(sj(2),li(1), 1);
	for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		wf.set(sj(1),li(1),ri(1), 0);
		wf.set(sj(2),li(1),ri(1), 1);
		for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	wf.set(sj(1),li(1), 0);
	wf.set(sj(2),li(1), 1);
	for( int d=3; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	(*psi).set(size,wf);
}


// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.
void
coherent_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);

	if(site==1)
		{
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weight_coherent_state(alpha, d-1));
		}
	else
		{
		wf.set(sj(1), li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		}
	(*psi).set(1,wf);

	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		if(j==site)
			{
			for( int d=1; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), weight_coherent_state(alpha, d-1));
			}
		else
			{
			wf.set(sj(1),li(1),ri(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), 0);
			}
		(*psi).set(j,wf);
	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);
	
	if(site==size)
		{
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weight_coherent_state(alpha, d-1));
		}
	else
		{
		wf.set(sj(1),li(1), 1);
		for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
		}
	
	(*psi).set(size,wf);


}


void
coherent_state_all_sites( MPS* psi, const SiteSet sites, const int size , complex<double> alpha )
{

	Index sj = sites(1);
	Index li = commonIndex((*psi)(1),(*psi)(2));
	ITensor wf = ITensor(sj,li);


	for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weight_coherent_state(alpha, d-1));
	(*psi).set(1,wf);


	for( int j = 2 ; j < size; j++)
	{
		sj = sites(j);
		li = commonIndex((*psi)(j-1),(*psi)(j));
		Index ri = commonIndex((*psi)(j),(*psi)(j+1));
		wf = ITensor(sj,li,ri);
		for( int d=1; d <= dim(sj); d++) wf.set(sj(d), li(1), ri(1), weight_coherent_state(alpha, d-1));
		(*psi).set(j,wf);

	}

	sj  = sites(size);
	li = commonIndex((*psi)(size-1),(*psi)(size));
	wf = ITensor(sj,li);


	for( int d=1; d <= dim(sj); d++) wf.set(sj(d),li(1), weight_coherent_state(alpha, d-1));
	(*psi).set(size,wf);
}


// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.
void
squeezed_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const double r )
{

	Index sj = sites(1);


	if( size > 1)
	{
		Index li = commonIndex((*psi)(1),(*psi)(2));
		ITensor wf = ITensor(sj,li);
		if(site==1)
			{
			for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  weight_squeezed_state(r, d-1));
			for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  0);
			}
		else
			{
			wf.set(sj(1),li(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
			}
		(*psi).set(1,wf);

		for( int j = 2 ; j < size; j++)
		{
			sj = sites(j);
			li = commonIndex((*psi)(j-1),(*psi)(j));
			Index ri = commonIndex((*psi)(j),(*psi)(j+1));
			wf = ITensor(sj,li,ri);
			if(j==site)
				{
				for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1),ri(1), weight_squeezed_state(r, d-1));
				for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),ri(1),  0);
				}
			else
				{
				wf.set(sj(1),li(1),ri(1), 1);
				for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
				}
			(*psi).set(j,wf);
		}

		sj  = sites(size);
		li = commonIndex((*psi)(size-1),(*psi)(size));
		wf = ITensor(sj,li);
		
		if(site==size)
			{
			for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),li(1), weight_squeezed_state(r, d-1));
			for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),li(1),  0);
			}
		else
			{
			wf.set(sj(1),li(1), 1);
			for( int d=2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
			}
		
		(*psi).set(size,wf);
	}

	else
	{
		ITensor wf = ITensor(sj);
		for( int d=1; d <= dim(sj); d+=2) wf.set(sj(d),  weight_squeezed_state(r, d-1));
		for( int d=2; d <= dim(sj); d+=2) wf.set(sj(d),  0);
		(*psi).set(size,wf);
	}


}


void
initial_state_cat_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha )
{
	MPS psi_t_1 = randomMPS(sites);
	MPS psi_t_2 = randomMPS(sites);
	coherent_state_site_j( &psi_t_1, sites, size , site, alpha );
	coherent_state_site_j( &psi_t_2, sites, size , site, -alpha );
	*psi = sqrt(0.5) * sum(psi_t_1,psi_t_2);
	(*psi).position(1);
	(*psi).normalize();
}


void
kink_state(MPS *psi, const SiteSet sites, const int number_ones)
{
	int L =  length(*psi);
	
	for(int j=1 ; j<= number_ones; j++)     put_occupation(psi,sites,j,1);
	for(int j=number_ones + 1 ; j<= L; j++) put_occupation(psi,sites,j,0);

}


void
put_occupation(MPS *psi, const SiteSet sites, const int position, const int n)
{
	Index sj;
	Index li;
	Index ri;
	ITensor wf;
	int L = length(*psi);
	if( position == 1)
	{
		sj = sites(1);
		li = commonIndex((*psi)(1),(*psi)(2));
		wf = ITensor(sj,li);
		for( int d=1; d <= n; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n+1),li(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}

	else if( position == L)
	{
		sj  = sites(L);
		li = commonIndex((*psi)(L-1),(*psi)(L));
		wf = ITensor(sj,li);

		for( int d=1; d <= n; d++) wf.set(sj(d),li(1), 0);
		wf.set(sj(n+1),li(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1), 0);
	}

	else
	{
		sj = sites(position);
		li = commonIndex((*psi)(position-1),(*psi)(position));
		ri = commonIndex((*psi)(position),(*psi)(position+1));
		wf = ITensor(sj,li,ri);

		for( int d=1; d <= n; d++) wf.set(sj(d),li(1),ri(1), 0);
		wf.set(sj(n+1),li(1),ri(1),1);
		for( int d=n+2; d <= dim(sj); d++) wf.set(sj(d),li(1),ri(1), 0);
		
	}

	(*psi).set(position,wf);

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
weight_coherent_state( const complex<double> alpha, const int k)
{
	complex<double> alpha_power_k = pow(alpha,k);
	return exp(-abs(alpha)*abs(alpha) / 2 ) * (alpha_power_k.real() + 1i*alpha_power_k.imag()) / sqrt(factorial(k));
}


// weight squeezed state

double 
weight_squeezed_state( const double r, const int k)
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
expectation_value_sigma_x_square( MPS *state , const SiteSet sites , const int j )
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
expectation_value_sigma_x( MPS *state , const SiteSet sites , const int j )
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
expectation_value_sigma_x_n( MPS *state , const SiteSet sites , const int j )
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
expectation_value_n_sigma_x( MPS *state , const SiteSet sites , const int j )
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
measure_imbalance( vector<double> &occupation_number, int k)
{

	double nk = occupation_number[k];
	occupation_number.erase(occupation_number.begin()+k);
	double nmax = *max_element(occupation_number.begin(), occupation_number.end());
	
	return (nk-nmax)/(nk+nmax);

}


double max_projector_at_cutoff(vector<vector<double> > &projector_all_sites,const int size,const int cut_off)
{
	double maximum = -1;
	for(int j = 1 ; j <= size ; j++) maximum = std::max(projector_all_sites[j-1][cut_off-1], maximum);

	return maximum;

}


// measure squareoccupation number along the 1D chain

void 
measure_square_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &square_occupation_number )
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

void measure_projector_all_sites( MPS *ground_state , const SiteSet sites , const int size , const int cut_off_fock_space ,  vector<vector<double> > &projector_all_sites ,  vector<double> &occupation_number)
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
measure_covariance_matrix_number_operator( MPS *psi , const SiteSet sites , vector<vector<double> > &covariance_matrix_NN_system, vector<vector<double> > &covariance_matrix_NN_gaussian, vector<vector<double> > &relative_error)
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
				NiNj    = abs(compute_two_point(psi,sites ,Ni_op,Nj_op , i , j));
				Aid_Ajd = compute_two_point(psi,sites,Adi_op,Adj_op , i , j);
				Ai_Aj   = compute_two_point(psi,sites,Ai_op,Aj_op , i , j);
				Aid_Aj  = compute_two_point(psi,sites,Adi_op,Aj_op  , i , j);
				Ai_Ajd  = compute_two_point(psi,sites,Ai_op,Adj_op  , i , j);	
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
measure_delta_x(MPS *psi, MPO *A, MPO *Adag )
{

	complex<double> a_expectation_value     = innerC(*psi, *A   , *psi);
	complex<double> aa_expectation_value    = innerC(*psi, *A   , *A, *psi);
	complex<double> adaga_expectation_value = innerC(*psi, *Adag, *A, *psi);

	complex<double> delta_x = 1 + 2 * adaga_expectation_value + 2 * aa_expectation_value.real() - 4 * a_expectation_value.real() * a_expectation_value.real();

	return delta_x.real();
}


double
measure_delta_p(MPS *psi, MPO *A, MPO *Adag )
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
