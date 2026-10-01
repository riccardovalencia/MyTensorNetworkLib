/**
 * @file entanglement.cc
 * @brief Implementation of entanglement.h (the functions are documented in the header).
 */
#include "entanglement.h"
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
// Compute entanglement entropy along the bond [site,site+1]

double
entanglement_entropy( MPS* psi , int site)
	{
	(*psi).position(site); 
	ITensor wf = (*psi)(site) * (*psi)(site+1);
	ITensor U  = (*psi)(site);
	ITensor S,V;
	auto spectrum = svd(wf,U,S,V);
	
	double SvN = 0.;
	for(auto p : spectrum.eigs())
		{
		if(p > 1E-12) SvN += -p*log2(p);
		}
	return SvN;
	}


// ----------------------------------------------------------
// Same as above, but with the natural logarithm (entropy in nats).
// Kept for the transverse-field Ising (spins) module; N is unused.

double
entanglement_entropy( MPS* psi , int N , int site)
	{
	(*psi).position(site); 
	ITensor wf = (*psi)(site) * (*psi)(site+1);
	ITensor U  = (*psi)(site);
	ITensor S,V;
	auto spectrum = svd(wf,U,S,V);

	double SvN = 0.;
	for(auto p : spectrum.eigs())
		{
		if(p > 1E-12) SvN += -p*log(p);
		}
	return SvN;
	}
