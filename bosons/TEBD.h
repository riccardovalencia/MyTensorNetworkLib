#ifndef TEBD_H
#define TEBD_H

#include <itensor/all.h>

using namespace std;
using namespace itensor;

//single time step for time evolution in bosonic quantum east model chain

void
build_single_step( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const int j );


//single time step for time evolution in bosonic quantum east model chain, where we specify n0 as a number

void
build_single_step_n0( ITensor *hterm , const SiteSet sites , const int size , const int n0, const double J , const double c , const int j );

// Build a full TEBD time step under the bosonic quantum east hamiltonian with only next-neighbour density-density interaction
// i.e. H = - 0.5 \sum_i n_i (exp(-s)\sigma_i+1^x - U n_i+1 - 1) where we do not fix any symmetry sector. 

void
build_TEBD_dt_step_H(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c, const string dynamics = "closed" , const double gamma = 0.);


// Perform adiabatic transformation from s=infty to s via a linear protocol
void
adiabatic_transformation_linear_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T);


// Perform adiabatic transformation from s=infty to s via aa tanh(x) protocol
void
adiabatic_transformation_tanh_protocol( MPS *psi_start, const SiteSet sites, const double s, const double c, const double dt, const double T);


//single time step for time evolution in bosonic quantum east model chain with losses
 
void
build_single_step_jumps( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const double gamma, const int j);

//generalization of the time evolution of the open bosonic quantum east model with any jump operator list.
 
void
build_single_step_jumps_v2( ITensor *hterm , const SiteSet sites , const int size , const double J , const double c , const int j, vector<ITensor> &Lj, vector<ITensor> &Ljd);

// build the TEBD gates for the open bosonic quantum east model with any list of jump operators

void
build_TEBD_dt_step_H_open(vector<BondGate> &gates, const SiteSet sites, const int size, const double dt, const double J, const double c , vector<ITensor> &Lj, vector<ITensor> &Ljd );


#endif
