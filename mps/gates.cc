/**
 * @file gates.cc
 * @brief Implementation of gates.h (the functions are documented in the header).
 */
#include "gates.h"
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


// My classes

TebdGate::TebdGate(const SiteSet sites, vector<int> j, const double dt, const ITensor h)
{
	jn_ = j;
	gate_ = expHermitian(h,-1_i * dt);
	sites_ = sites;
}


ITensor TebdGate::gate()
{
	return gate_;
}


vector<int> TebdGate::sites()
{
	return jn_;
}


void TebdGate::set_gate(ITensor new_gate)
{
	gate_ = new_gate;
}


DissipativeGate::DissipativeGate(const SiteSet sites, vector<int> jket, vector<int> jbra, double dt, ITensor h)
{
	jnket_ = jket;
	jnbra_ = jbra;
	sites_ = sites;

	// linear approximation - tested and works well for our purposes
	gate_ = h * dt;


}


ITensor DissipativeGate::gate()
{
	return gate_;
}


vector<int> DissipativeGate::ket_sites()
{
	return jnket_;
}


vector<int> DissipativeGate::bra_sites()
{
	return jnbra_;
}


OperatorPair::OperatorPair(vector<ITensor> T, vector<int> j, double gamma)
{
	Ti_   = T[0];
    Tj_   = T[1];
    i_    = j[0];
    j_    = j[1];
	gamma_ = gamma; 
}


ITensor OperatorPair::op_i()
{
	return Ti_;
}


ITensor OperatorPair::op_j()
{
	return Tj_;
}


int OperatorPair::site_i()
{
	return i_;
}


int OperatorPair::site_j()
{
	return j_;
}


double OperatorPair::rate()
{
	return gamma_;
}


// ----------------------------------------------------------
// Apply a gate on the jn sites of an MPS psi 
// It is possible to apply up to 3-sites gates.

MPS
apply_gate(MPS psi, const ITensor gate, const vector<int> jn, const Args args)
{
    double cut_off = args.getReal("Cutoff");
    int maxDim  = args.getInt("MaxDim");
	int j = jn[0];

	if( j < 0 || j > length(psi))
	{
		cerr << "Site not valid. Cannot apply gate \n";
		exit(0);
	}

	ITensor AA = gate;
	psi.position(j);
           
	for(int q : jn) AA *= psi(q);  
	AA.mapPrime(1,0);

	if(jn.size() == 1)
	{
		psi.set(j,AA);
	}

	else if(jn.size() == 2)
	{
		auto [U,S,V] = svd(AA,inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
		psi.set(j,U);
		psi.set(j+1,S*V);
	}

	else if(jn.size() == 3)
	{
		auto [U,S,V] = svd(AA,inds(psi(j)),{"Cutoff=",cut_off,"MaxDim=",maxDim});
		Index l =  commonIndex(U,S);
		psi.set(j,U);

		auto [U1,S1,V1] = svd(S*V,{inds(psi(j+1)),l},{"Cutoff=",cut_off,"MaxDim=",maxDim});
		psi.set(j+1,U1);
		psi.set(j+2,S1*V1);
	}

	else
	{
	cerr << jn.size() <<"-gate not implemented (yet)!\n";
	exit(0);
	}
    
	return psi;
}
