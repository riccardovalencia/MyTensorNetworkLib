/**
 * @file adiabatic.cc
 * @brief Implementation of adiabatic.h (the functions are documented in the header).
 */
#include "adiabatic.h"
#include "../models/bqem.h"
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


}
