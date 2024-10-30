#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <itensor/all.h>
#include <vector>

using namespace itensor;
using namespace std;

double
entanglement_entropy( MPS* , int  );

void 
measure_mx_mz( const SpinHalf , MPS , const int );

// expectation value: <\sigma_j^x^2>
double 
expectation_value_sigma_x_square( MPS *state , const SiteSet sites , const int j );

// expectation value: <\sigma_j^x>
double 
expectation_value_sigma_x( MPS *state , const SiteSet sites , const int j );

// expectation value: <\sigma_j^x n_j>
double 
expectation_value_sigma_x_n( MPS *state , const SiteSet sites , const int j );

// expectation value: <\sigma_j^x n_j>
double 
expectation_value_n_sigma_x( MPS *state , const SiteSet sites , const int j );

void 
measure_occupation_number( MPS * , const SiteSet , const int size ,  vector<double> &occupation_number);

double
measure_imbalance( vector<double> &occupation_number, int k);

double
max_projector_at_cutoff(vector<vector<double> > &projector_all_sites,const int size,const int cut_off);

void
measure_covariance_matrix_number_operator( MPS * , const SiteSet , vector<vector<double> > &, vector<vector<double> > &, vector<vector<double> > &);


complex<double> 
compute_two_point( MPS *, const SiteSet , ITensor op_i, ITensor op_j, int i, int j);


void 
measure_square_occupation_number( MPS *ground_state , const SiteSet sites , const int size ,  vector<double> &square_occupation_number );


void 
measure_projector_all_sites( MPS * , const SiteSet , const int size , const int cut_off_fock_space , vector<vector<double> > &projector_all_sites ,  vector<double> &occupation_number);

double 
compute_variance_H_mmGcbQEM(MPS *psi , const SiteSet sites, int size , int cut_off_fock_space, int n0, int symmetry_sector, double s, double c);


// compute variance of bosonic x = a + adag
double
measure_delta_x(MPS *psi, MPO *A, MPO *Adag );

// compute variance of bosonic p = -i(a + adag)

double
measure_delta_p(MPS *psi, MPO *A, MPO *Adag );

// compute squeezing on site j of a bosonic system
double
measure_squeezing(MPS *psi, const SiteSet sites, const int j);


// compute squeezing using the dressed annihilation and creation operators
double
measure_dressed_squeezing(MPS *psi, const SiteSet sites, MPO A, MPO N);
#endif
