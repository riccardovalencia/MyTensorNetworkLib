/**
 * @file adiabatic.cc
 * @brief Implementation of adiabatic.h (the functions are documented in the header).
 */
#include "adiabatic.h"
#include "../mps/gates.h"
#include "../models/bosonic_east_model.h"
#include <itensor/all.h>
#include <cmath>
#include <vector>

using namespace std;
using namespace itensor;


// One TEBD step of length dt of the bosonic east model with hopping J, followed by a normalization
static void
apply_ramp_step( MPS *psi, const SiteSet& sites, const double J, const double U, const double dt, const Args& args )
{
	*psi = apply_gates(*psi, make_bosonic_east_model_gates(sites, length(*psi), dt, J, U), args);
	(*psi).position(1);
	(*psi).normalize();
}


// linear protocol J(t) = J_target t / T: the larger T, the slower the ramp

void
evolve_adiabatic_linear_ramp( MPS *psi_start, const SiteSet sites, const double J_target, const double U, const double dt, const double T, const Args& args)
{
	double t = 0;
	double J = 0;
	do
	{
		t += dt;
		J = J_target * t / T;
		apply_ramp_step(psi_start, sites, J, U, dt, args);
	} while( J_target > J );
}


// tanh protocol J(t) = J_target tanh(t / T), until J is within a relative tolerance of J_target

void
evolve_adiabatic_tanh_ramp( MPS *psi_start, const SiteSet sites, const double J_target, const double U, const double dt, const double T, const Args& args)
{
	const double tolerance = 1E-6;
	double t = 0;
	double J = 0;
	do
	{
		t += dt;
		J = J_target * tanh(t/T);
		apply_ramp_step(psi_start, sites, J, U, dt, args);
	} while( (J_target - J)/(J_target + J) > tolerance );
}
