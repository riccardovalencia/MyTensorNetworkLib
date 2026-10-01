/**
 * @file spin_chain.cc
 * @brief Implementation of spin_chain.h (the functions are documented in the header).
 */
#include "spin_chain.h"
#include "../models/bosonic_east_model.h"
#include "../mps/gates.h"
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


// -----------------------------------------------------------------
// Gates of exp(-i dt H) with spin Hamiltonian H
//  H  =  - hx \sum_j X_j - Jxx \sum_j X_j X_{j+1} 
//        - hy \sum_j X_j - Jyy \sum_j Y_j Y_{j+1}  
//        - hz \sum_j X_j - Jzz \sum_j Z_j Z_{j+1}
// where J = [Jxx,Jyy,Jzz] and h = [hx,hy,hz]

vector<TebdGate>
make_spin_chain_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt)
{

	int N = length(sites);

    vector<TebdGate> gates;

	cerr << "vector J = (J_xx, J_yy , J_zz)\n";
	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];
	cerr << Jxx << "\n" << Jyy << "\n" << Jzz << "\n";

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];

	for(int j=1 ; j <= N-1 ; j+=1)
	{
		vector<ITensor> X;
		vector<ITensor> Y;
		vector<ITensor> Z;
		vector<ITensor> Id;

		for(int q=j ; q<=j+1; q++)
		{
			Id.push_back(      op(sites,"Id",q) );
			X.push_back(  2 * op(sites,"Sx",q) );
			Y.push_back(  2 * op(sites,"Sy",q) );
			Z.push_back(  2 * op(sites,"Sz",q) );
		}

		ITensor H_S , H_SS;

		if(j==1) H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] ;
		else     H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] / 2.;

		if(j <  N-1) H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) / 2. ;
		else         H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1])      ;

		H_SS  = Jxx * X[0]*X[1] + Jyy * Y[0]*Y[1] + Jzz * Z[0]*Z[1];		
		
		ITensor H = H_S + H_SS;

		vector<int> jn = {j,j+1};
		TebdGate g = TebdGate(jn,dt/2.,H);
		gates.push_back(g);
	}
	
	vector<TebdGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());
	for(TebdGate gate : gates_) gates.push_back(gate);

	return gates;
}


// -----------------------------------------------------------------
// Gates of exp(-i dt H) with spin Hamiltonian H
//  H  =  - hx \sum_j X_j - Jxx \sum_j X_j X_{j+1} 
//        - hy \sum_j X_j - Jyy \sum_j Y_j Y_{j+1}  
//        - hz \sum_j X_j - Jzz \sum_j Z_j Z_{j+1}
// where J = [Jxx,Jyy,Jzz] and h = [hx,hy,hz]
// Same as the one retuning <TebdGate>: overload of the function. 
// Depending on the degree of flexibility and control needed could be better to use one over the other



// effective Hamiltonian of a local dissipative process
// The coherent dynamics is given by the generic short range spin model

vector<TebdGate>
make_spin_chain_effective_gates(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const vector<int> Lj_sites, const vector<double> gamma, const double dt)
{

	int N = length(sites);

    vector<TebdGate> gates;

	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];
	for(int j=1 ; j <= N-1 ; j+=1)
	{
		vector<ITensor> X;
		vector<ITensor> Y;
		vector<ITensor> Z;
		vector<ITensor> Id;

		for(int q=j ; q<=j+1; q++)
		{
			Id.push_back(      op(sites,"Id",q) );
			X.push_back(  2 * op(sites,"Sx",q) );
			Y.push_back(  2 * op(sites,"Sy",q) );
			Z.push_back(  2 * op(sites,"Sz",q) );
		}

		ITensor H_S , H_SS;

		if(j==1) H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] ;
		else     H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] / 2.;

		if(j <  N-1) H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) / 2. ;
		else         H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1])      ;

		for(int k=0 ; k < (int)Lj_sites.size(); k++)
		{
			if(Lj_sites[k]==j)
			{
				ITensor lj  = Lj[k];
				ITensor ljd = conj(lj);
				ljd.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)
				ITensor ljdlj = ljd * lj; 
				ljdlj.mapPrime(2,1);
				H_S -= 0.5 * gamma[k] * Cplx_i * ljdlj * Id[1]; 
			}
			if(j==N-1 && Lj_sites[k] == N)
			{
				ITensor lj = Lj[k];
				ITensor ljd = conj(lj);
				ljd.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)
				ITensor ljdlj = ljd * lj; 
				ljdlj.mapPrime(2,1);
				H_S -= 0.5 * gamma[k] * Cplx_i * Id[0] * ljdlj;
			}
		}

		H_SS  = Jxx * X[0] * X[1];
		H_SS += Jyy * Y[0] * Y[1];
		H_SS += Jzz * Z[0] * Z[1];		
		
		ITensor H = H_S + H_SS;

		vector<int> jn = {j,j+1};
		TebdGate g = TebdGate({j,j+1}, BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H).gate()); 
		gates.push_back(g);
	}
	
	vector<TebdGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());
	for(TebdGate gate : gates_) gates.push_back(gate);

	return gates;
}


vector<TebdGate>
make_local_field_gates(const SiteSet sites , vector<double> omegaj, const double dt)
{

	int N = length(sites);
	vector<TebdGate> gates;

	// H = \sum_j omveja
	for(int j=1 ; j <= N; j++)
	{
		ITensor Sx = 2*op(sites,"Sx",j);
		ITensor Sy = 2*op(sites,"Sy",j);
		ITensor Sz = 2*op(sites,"Sz",j);

		ITensor hj = omegaj[0] * Sx + omegaj[1] * Sy + omegaj[2] * Sz;

		vector<int> jn = {j};

		TebdGate g = TebdGate(jn,dt/2.,hj);
		gates.push_back(g);
	}


	vector<TebdGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(TebdGate gate : gates_) gates.push_back(gate);
	
	return gates;
}


//----------------------------------------------------------------------

//single gate acting on sites [b,b+1] of the Ising model with longitudinal (hx) and transversal (hz) magnetic fields

ITensor
make_ising_bond_hamiltonian( const SpinHalf sites , const int N , const double J , const double hx , const double hz , const int b )
	{
    ITensor hterm;
	ITensor Sx1 = sites.op("Sx",b);					
	ITensor Sx2 = sites.op("Sx",b+1);
	ITensor Sz1 = sites.op("Sz",b);
	ITensor Sz2 = sites.op("Sz",b+1);
	ITensor Id1 = sites.op("Id",b);
	ITensor Id2 = sites.op("Id",b+1);
		
	hterm = - 4 * J * Sx1 * Sx2;
		
	if( b == 1 )
		{
		hterm +=  - 2 * J * hx * ( Sx1 * Id2 + Id1 * Sx2 / 2. ); 									
		hterm +=  - 2 * J * hz * ( Sz1 * Id2 + Id1 * Sz2 / 2. );	
		}
	else if( b == N-1)
		{
		hterm +=  - 2 * J * hx * ( Sx1 * Id2 / 2. + Id1 * Sx2 ); 									
		hterm +=  - 2 * J * hz * ( Sz1 * Id2 / 2. + Id1 * Sz2 );		
		}
	else{
		hterm +=  - 2 * J * hx * ( Sx1 * Id2 + Id1 * Sx2 ) / 2.; 									
		hterm +=  - 2 * J * hz * ( Sz1 * Id2 + Id1 * Sz2 ) / 2.;	
		}
    return hterm;
}


// ----------------------------------------------------------
// second-order Trotter step of the Ising chain: forward sweep with dt/2, then the reversed sweep

vector<TebdGate>
make_ising_gates( const SpinHalf sites , const int N , const double J , const double hx , const double hz , const double dt )
{
    vector<TebdGate> gates;
    for(int b = 1 ; b <= N-1 ; b++)
    {
        ITensor hterm = make_ising_bond_hamiltonian(sites, N, J, hx, hz, b);
        gates.push_back(TebdGate({b, b+1}, BondGate(sites, b, b+1, BondGate::tReal, dt/2., hterm).gate()));
    }
    for(int b = N-1 ; b >= 1 ; b--) gates.push_back(gates[b-1]);
    return gates;
}
