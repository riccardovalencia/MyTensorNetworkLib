/**
 * @file tight_binding.cc
 * @brief Implementation of tight_binding.h (the functions are documented in the header).
 */
#include "tight_binding.h"
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


MPO
make_tight_binding_mpo(const int N , const SiteSet sites, const vector<double> J, const vector<double> hup, const vector<double> hdn)
{
    int size_J = J.size();
    int size_hup = hup.size();
    int size_hdn = hdn.size();


    auto ampo = AutoMPO(sites);

    for(int j : range1(N-1))
    {
        
        ampo += J[(j-1)%size_J] , "Cdagup" , j   , "Cup" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagup" , j+1 , "Cup" , j  ;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j   , "Cdn" , j+1;
        ampo += J[(j-1)%size_J] , "Cdagdn" , j+1 , "Cdn" , j  ;
        
    }

    // on-site fields on every site
    for(int j : range1(N))
    {
        ampo += hup[(j-1)%size_hup] , "Nup" , j;
        ampo += hdn[(j-1)%size_hdn] , "Ndn" , j;
    }


    return toMPO(ampo);
}


ITensor
make_single_particle_hamiltonian(const int N , const vector<double> J, const vector<double> h, const bool spinful)
{
    int L = N;
    vector<double> Jp = J;
    vector<double> hp = h; 

    if(spinful) L = 2 * N;

    Index r = Index(L);
    Index l = prime(r);

    ITensor H = ITensor(r,l);

    for(int i=1 ; i <= dim(l) ; i++)
    {
        for(int j=1 ; j <= dim(r) ; j++)
        {   
            if(i==j) H.set(l(i),r(j),hp[0]);
            if(abs(i-j)==1)  H.set(l(i),r(j), Jp[0]);
        }
    } 

    return H;

}


ITensor
make_single_particle_hamiltonian_impurity(const int N , const vector<double> J, const vector<double> h, const bool spinful)
{
    int L = N;
    vector<double> Jp = J;
    vector<double> hp = h; 

    if(spinful) L = 2 * N;

    Index r = Index(L);
    Index l = prime(r);

    ITensor H = ITensor(r,l);

    for(int i=1 ; i <= dim(l) ; i++)
    {
        for(int j=1 ; j <= dim(r) ; j++)
        {   
            if(i==j) H.set(l(i),r(j),hp[0]);
            if(abs(i-j)==1)  
            {
                if(j==1 || i==1) H.set(l(i),r(j), Jp[0]);
                else H.set(l(i),r(j), Jp[1]);
            }

        }
    } 

    return H;

}


// Simulation of free spinful fermions 

vector<TebdGate>
make_free_fermion_gates(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{ 


	int N = length(sites);

    vector<TebdGate> gates;
	
	for(int j=1 ; j <= N-1 ; j+=1)
	{


		vector<ITensor> Adagup;
		vector<ITensor> Aup;

		vector<ITensor> Adagdn;
		vector<ITensor> Adn;

		vector<ITensor> Nup;
		vector<ITensor> Ndn;

		vector<ITensor> Id;

		vector<ITensor> AdagupFi;
		vector<ITensor>	AupFi;
		vector<ITensor>	FiAdn;
		vector<ITensor>	FiAdagdn;


		for(int q=j ; q<=j+1; q++)
		{
			Adagup.push_back(op(sites,"Adagup",q) );
			Aup.push_back(   op(sites,"Aup",q) );
			Adagdn.push_back(op(sites,"Adagdn",q) );
			Adn.push_back(   op(sites,"Adn",q) );
			Nup.push_back(   op(sites,"Nup",q) );
			Ndn.push_back(   op(sites,"Ndn",q) );
			Id.push_back(op(sites,"Id",q));

			// see https://itensor.org/docs.cgi?page=tutorials/fermions
		
			AdagupFi.push_back(op(sites,"Adagup*F",q));
			AupFi.push_back(op(sites,"Aup*F",q));
			FiAdn.push_back(op(sites,"F*Adn",q));
			FiAdagdn.push_back(op(sites,"F*Adagdn",q));


		}

		ITensor H, H_S , H_SS;

		// build single site
	
		if(j ==1)    H_S = (hup[j-1] * Nup[0] + hdn[j-1] * Ndn[0]) * Id[1] ;
		else         H_S = (hup[j-1] * Nup[0] + hdn[j-1] * Ndn[0]) * Id[1] / 2.;

		if(j == N-1) H_S += Id[0] * (hup[j] * Nup[1] + hdn[j] * Ndn[1]) ;
		else         H_S += Id[0] * (hup[j] * Nup[1] + hdn[j] * Ndn[1]) / 2.      ;

		// two sites hopping

		H_SS   = J[j-1] * (  AdagupFi[0] * Aup[1]   - AupFi[0] * Adagup[1] );
		H_SS  += J[j-1] * (  Adagdn[0]   * FiAdn[1] - Adn[0] * FiAdagdn[1] );


		H = H_S + H_SS ;
		TebdGate g = TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()); 
		gates.push_back(g);
			
	}


	vector<TebdGate> gates_reversed = gates;
	
	reverse(gates_reversed.begin(), gates_reversed.end());

	for(TebdGate g : gates_reversed) gates.push_back(g);

	return gates;
}
