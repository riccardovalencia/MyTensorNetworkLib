#ifndef INITIAL_STATE
#define INITIAL_STATE

#include <itensor/all.h>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
//all spins UP along X
void 
initial_state_all_UP( const SpinHalf , MPS* , const int );

//----------------------------------------------------------------------
//all spins DOWN along X
void
initial_state_all_DOWN( const SpinHalf , MPS* , const int );

//----------------------------------------------------------------------
//half chain DOWN along X, half chain UP along X (domain wall with single kink)
void
initial_state_DOMAIN_WALL( const SpinHalf , MPS* , const int );

//----------------------------------------------------------------------
//Initial bosonic state of the form |00..0> |n0> |0000...0>, where n0 is in the position `excitation_position`
void
initial_state_n0_excitation( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position);

//----------------------------------------------------------------------
// Insert n0 excitation on site excitation_position without changing the other sites
void
initial_state_n0_excitation_pinned( MPS* psi, const SiteSet sites, const int size , const int n0 , const int excitation_position);


// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.
void
initial_state_vacuum_state( MPS* psi, const SiteSet sites, const int size );

// ----------------------------------------------------------------------
// Initial vacuum bosonic state |0000...0>.
void
initial_state_vacuum_state_correct_link( MPS* psi, const SiteSet sites, const int size  );

// ----------------------------------------------------------------------
// Initial vacuum bosonic state |1111...1>.
void
initial_state_all_one_state_correct_link( MPS* psi, const SiteSet sites, const int size) ;

// Initial bosonic state of the form |000..0> |alpha>_j |0000...0>, such that a|alpha> = alpha|alpha>.
void
coherent_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha);

// initial bosonic state: product state of coherent states |\alpha>
void
coherent_state_all_sites( MPS* psi, const SiteSet sites, const int size , const complex<double> alpha );

// Initial bosonic state of the form |000..0> |r>_j |0000...0>, such that |r> is squeezed.
void
squeezed_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const double r );

// Initial single-body even cat state on site 'site'

void
initial_state_cat_state_site_j( MPS* psi, const SiteSet sites, const int size , const int site, const complex<double> alpha );


// Initial state |0>^k \otimes |n_0> \otimes |GS(n_0)_L> \otimes |0>^(N-L-k-1)
// It's a state with k 0's, n_0, the ground state of size L in this symmetry sector, and the other 0's.
// The total size of the system is L

MPS
super_bosonic_state(MPS ground_state, const SiteSet sites, const SiteSet ground_sites, const int L, const int N, const int n0, const int k);

// Insert a state within another state, such that you have a state |ground_state> that you want to put in another state |psi_t0> from site start to start+L

void
insert_state(MPS* psi_t0, MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int start, const int L, const int N);


// Return a super-bosonic coherent state via dressing of a single site coherent state
MPS
super_bosonic_coherent_state( MPS *psi_coherent, const SiteSet sites_coherent, const double alpha, const double s, const double c, double dt, double T);


// return the super-bosonic coherent state built as linear combinations of DMRG obtained ground state
MPS
super_bosonic_coherent_state_from_ground_states( MPS *psi_coherent, const SiteSet sites, const complex<double> alpha, const double s, const double c, double dt, double T,Args const& physical_args, Args const& numerical_args);


// Return a super-bosonic squeezed state
MPS
super_bosonic_squeezed_state( MPS *psi_coherent, const SiteSet sites_coherent, const double alpha, const double s, const double c, double dt, double T);


// create a state of the form |1111>|0000> with number_ones 1.
void
kink_state(MPS *psi, const SiteSet sites, const int number_ones);

// put n particles on site position
void
put_occupation(MPS *psi, const SiteSet sites, const int position, const int n);


// Return n!
double
factorial(int n);

// return the weight exp(-alpha^2/2)alpha^k/sqrt(k!) for coherent states

complex<double> 
weigth_coherent_state( const complex<double> alpha, const int k);

// return weight squeezed state
double 
weigth_squeezed_state( const double r, const int k);


// return the dressed version of a given operator O. The dressing is performed via the adiabatic time evolution via the Hamiltonian of the
// bosonic quantum east model

MPO
bQEM_dressed_operator(MPO O,  const SiteSet sites, const double s, const double c, double dt, double T);
#endif
