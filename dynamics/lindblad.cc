/**
 * @file lindblad.cc
 * @brief Implementation of lindblad.h (the functions are documented in the header).
 */
#include "lindblad.h"
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


// using this, we have the more standard usage of indices, but we have to use the ITensor convention
// for S^- , which corresponds to the standard S^+ convention. 
// CHECK: 03.09.23 IT COULD BE THAT I AM MESSING UP INDICES TECHNICALLY, SINCE I AM ASSOCIATING THE INDEX 0 TO THE KET,
// AND PRIME TO THE BRA. BUT FROM ITENSOR DEFAULT CONVENTION IT COULD BE THAT THEY ARE SWAPPED. THIS IS WHY IT LOOKS LIKE
// MY CONVENTION OF S^- IS THE OPPOSITE OF THE ONE OF ITENSOR (IN REALITY THEY ARE NOT DIFFERENT). I SHOULD CHECK THIS
// REWRITING A PIECE OF CODE CONCERNING THIS AND TESTING WITH A NON-HERMITIAN JUMP.
vector<MyBondGateDiss>
gates_local_lindblad(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt)
{

	int N = length(sites);
	vector<MyBondGateDiss> gates;

	// ket has index sj
	// bra has index sj'
	// we need a gate  input (sj,sj') -> (sj'',sj''') output
	// 		   _	
	// 	sj''- | | - sj
	// 		  |	|
	// sj'''- | | - sj'
	// 	

	//  The procedure is: (0,1) -> (2,3) (input-output)
	//   sj''
	//   |
	//   Lj - 1/2 (Ljdag \otimes Lj)
	//   | sj
	//   o-
	//   | sj'
	//   LjdagT - 1/2 (LjdagT \otimes LjT)
	//   |
	//   sj'''

	for(int j : lj_sites)
	{
		Index sj  = sites(j);
		Index sj1 = prime(sj);
		Index sj2 = prime(sj,2);
		Index sj3 = prime(sj,3);

		ITensor Idket = ITensor(sj,sj2);
		ITensor Idbra = ITensor(sj1,sj3);

		for(int q=1 ; q<=dim(sj) ; q++)
		{
			Idket.set(sj(q),sj2(q),1.);
			Idbra.set(sj1(q),sj3(q),1.);
		}


		ITensor lj = Lj[j-1];
		ITensor ljd = conj(lj);
		ITensor lj_ = lj;
		ITensor ljd_ = ljd;

		// non-hermitian Hamiltonian

		// ket 
		ljd_.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)
		ITensor ljdlj_I = ljd_ * lj_ * Idbra; 

		// reset
		lj_ = lj;
		ljd_ = ljd;
	
		// bra
		ljd_.mapPrime(0,3); 
		ITensor I_ljdlj = ljd_ * lj_; // acts on bra  - (3,0)
		I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
		I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
		
		// reset
		lj_ = lj;
		ljd_ = ljd;

		// jumps 
		lj_.mapPrime(1,2);
		ljd_.mapPrime(1,3);
		ljd_.mapPrime(0,1);
		ITensor lj_ljd = lj_ * ljd_;


		ITensor Dj = gammaj[j-1] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

		vector<int> jnket = {j};
		vector<int> jnbra = {j};

		MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,Dj);
	
		gates.push_back(g);

	}

	return gates;
}


// using this, we have the more standard usage of indices, but we have to use the ITensor convention
// for S^- , which corresponds to the standard S^+ convention. 

vector<MyBondGateDiss>
gates_local_lindblad(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
{

	int N = length(sites);
	vector<MyBondGateDiss> gates;

	// ket has index sj
	// bra has index sj'
	// we need a gate  input (sj,sj') -> (sj'',sj''') output
	// 		   _	
	// 	sj''- | | - sj
	// 		  |	|
	// sj'''- | | - sj'
	// 	

	//  The procedure is: (0,1) -> (2,3) (input-output)
	//   sj''
	//   |
	//   Lj - 1/2 (Ljdag \otimes Lj)
	//   | sj
	//   o-
	//   | sj'
	//   LjdagT - 1/2 (LjdagT \otimes LjT)
	//   |
	//   sj'''


	for(int j=1 ; j <= N; j++)
	{
		
		Index sj  = sites(j);
		Index sj1 = prime(sj);
		Index sj2 = prime(sj,2);
		Index sj3 = prime(sj,3);

		ITensor Idket = ITensor(sj,sj2);
		ITensor Idbra = ITensor(sj1,sj3);

		for(int q=1 ; q<=dim(sj) ; q++)
		{
			Idket.set(sj(q),sj2(q),1.);
			Idbra.set(sj1(q),sj3(q),1.);
		}


		ITensor lj = Lj[j-1];
		ITensor ljd = conj(lj);
		ITensor lj_ = lj;
		ITensor ljd_ = ljd;

		// non-hermitian Hamiltonian

		// ket 
		ljd_.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)
		ITensor ljdlj_I = ljd_ * lj_ * Idbra; 

		// reset
		lj_ = lj;
		ljd_ = ljd;
	
		// bra
		ljd_.mapPrime(0,3); 
		ITensor I_ljdlj = ljd_ * lj_; // acts on bra  - (3,0)
		I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
		I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
		
		// reset
		lj_ = lj;
		ljd_ = ljd;

		// jumps 
		lj_.mapPrime(1,2);
		ljd_.mapPrime(1,3);
		ljd_.mapPrime(0,1);
		ITensor lj_ljd = lj_ * ljd_;


		ITensor Dj = gammaj[j-1] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

		vector<int> jnket = {j};
		vector<int> jnbra = {j};

		MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,Dj);
	
		gates.push_back(g);

	}


	return gates;
}


// keeping as backup - 4.05.23
// vector<MyBondGateDiss>
// gates_local_lindblad(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
// {

// 	int N = length(sites);
// 	vector<MyBondGateDiss> gates;

// 	// ket has index sj
// 	// bra has index sj'
// 	// we need a gate  input (sj,sj') -> (sj'',sj''') output
// 	// 		   _	
// 	// 	sj''- | | - sj
// 	// 		  |	|
// 	// sj'''- | | - sj'
// 	// 	

// 	//  The procedure is: (0,1) -> (2,3) (input-output)
// 	//   sj''
// 	//   |
// 	//   Lj - 1/2 (Ljdag \otimes Lj)
// 	//   | sj
// 	//   o-
// 	//   | sj'
// 	//   LjdagT - 1/2 (LjdagT \otimes LjT)
// 	//   |
// 	//   sj'''


// 	// vector<MyBondGate> gates_ = gates;
// 	// reverse(gates_.begin(), gates_.end());

// 	// for(MyBondGate gate : gates_) gates.push_back(gate);
	
// 	return gates;
// }

// Non-diagonal local Lindland 
// We consider Lindbland of the form: Li \rho L_{i+1}^\dagger + 0.5 * {Li L_{i+1}, \rho}
// It appears as a 4-sites gate if we apply the jumps and the non-hermitian part at the same time.
// A possibility is to split differently: the non hermitiain part in the hermitian one
// In this way I have at most 2-sites gates. The number of gates that have to be applied is the same.  <- POSSIBLE EFFICIENCY GAIN(?)
// Drawback: we would have long-range interactions in the final case study both in the Hamiltonian part 
// and jump part -> MULTIPLE LOOPS NECESSARY                                                           -> HUGE INEFFICENCY FROM LOOPING

// MyTrainITensor is a personalized class containing ITensors which have to act either on 


// It is a 4-sites object
// // We apply Li \rho L_j^\dagger + L_j \rho L_i^\dagger - 1/2( {L_i^\dagger L_j , \rho} + {L_j^\dagger L_i,\rho} ) (OR SIMILAR)


vector<MyBondGateDiss>
gates_nearest_neighbour_local_lindblad(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt)
{

	int N = length(sites);
	vector<MyBondGateDiss> gates;
	// check size of the two containers

	for(MyTrainITensor T : TTrain)
	{
		
		int i = T.i();
		int j = T.j();

		ITensor li = T.Ti();
		ITensor lj = T.Tj();
		ITensor lid = dag(li);
		ITensor ljd = dag(lj);

		double gamma = T.gamma();

		// site index

		cerr << "Sites : " << i << " " << j << "\n";

		if( abs(i-j)>1 )
		{
			cerr << "non-local dissipation still not implemented! Returning empty set of gates.\n";
			return gates;
		}

		if(abs(i-j)==0)
		{

			ITensor lj_ = lj;
			ITensor ljd_ = ljd;

			Index sj  = sites(j);
			Index sj1 = prime(sj);
			Index sj2 = prime(sj,2);
			Index sj3 = prime(sj,3);

			ITensor Idket = ITensor(sj,sj2);
			ITensor Idbra = ITensor(sj1,sj3);

			for(int q=1 ; q<=dim(sj) ; q++)
			{
				Idket.set(sj(q),sj2(q),1.);
				Idbra.set(sj1(q),sj3(q),1.);
			}

			// non-hermitian Hamiltonian

			// ket
			ljd_.mapPrime(0,2);
			ITensor ljdlj_I = ljd_ * lj_ * Idbra; // acts on ket (sj,sj') -> (sj'',sj''') as desired
		
			// reset
			lj_ = lj;
			ljd_ = ljd;

			// bra
			ljd_.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
			ITensor I_ljdlj = ljd_ * lj_; // acts on bra  - (3,0)
			I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
			I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
			I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
			
			// reset
			lj_ = lj;
			ljd_ = ljd;
			
			// jumps 
			lj_.mapPrime(1,2);
			ljd_.mapPrime(1,3);
			ljd_.mapPrime(0,1);
			ITensor lj_ljd = lj_ * ljd_;

			// all together

			ITensor Dj = gamma * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);
			vector<int> jnket = {j};
			vector<int> jnbra = {j};

			MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dj);
		
			gates.push_back(g);

		}

		if(abs(i-j)==1)
		{
			Index si  = sites(i);
			Index si1 = prime(si);
			Index si2 = prime(si,2);
			Index si3 = prime(si,3);

			Index sj  = sites(j);
			Index sj1 = prime(sj);
			Index sj2 = prime(sj,2);
			Index sj3 = prime(sj,3);

			ITensor Idket = ITensor(si,sj,si2,sj2);
			ITensor Idbra = ITensor(si1,sj1,si3,sj3);

			ITensor Id_iket_jbra = ITensor(si,si2,sj1,sj3);
			ITensor Id_ibra_jket = ITensor(sj,sj2,si1,si3);
						
			for(int q=1 ; q<=dim(sj) ; q++)
			{
				Idket.set( si(q) , si2(q) , sj(q)  , sj2(q) , 1.);
				Idbra.set(si1(q) , si3(q) , sj1(q) , sj3(q) , 1.);
				Id_iket_jbra.set(si(q),si2(q),sj1(q),sj3(q) , 1.);
				Id_ibra_jket.set(sj(q),sj2(q),si1(q),si3(q) , 1.);
			}

			// effective Hamiltonian part - jumps act either on the ket or bra, but not both

			// acts on ket

			ITensor lid_ = lid;
			ITensor ljd_ = ljd;
			ITensor li_  = li;
			ITensor lj_  = lj;

			// ket

			lid_.mapPrime(0,2);
			lid_.mapPrime(1,0);
			ljd_.mapPrime(0,2);
			ljd_.mapPrime(1,0);
			li_.mapPrime(1,2);
			lj_.mapPrime(1,2);
			
			ITensor ljd_li_I = ljd_ * li_ ;
			ITensor lid_lj_I = lid_ * lj_ ;
			ITensor LdL_I = (ljd_li_I + lid_lj_I ) * Idbra ; // acts on ket (sj,sj') -> (sj'',sj''') as desired


			// reset
			lid_ = lid;
			ljd_ = ljd;
			li_  = li;
			lj_  = lj;

			// bra
			ljd_.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
			lid_.mapPrime(0,3);
			li_.mapPrime(1,3);
			lj_.mapPrime(1,3);
			li_.mapPrime(0,1);
			lj_.mapPrime(0,1);


			ITensor I_LdL = ljd_ * li_ + lid_ * lj_; // acts on bra  - (3,0)
			I_LdL = swapPrime(I_LdL,1,3); // transposition
			I_LdL *= Idket;			// acts on bra from (0,1) to (2,3) as desired
			
			// reset

			lid_ = lid;
			ljd_ = ljd;
			li_  = li;
			lj_  = lj;

			// jumps

			lj_.mapPrime(1,2);
			ljd_.mapPrime(1,3);
			ljd_.mapPrime(0,1);

			li_.mapPrime(1,2);
			lid_.mapPrime(1,3);
			lid_.mapPrime(0,1);

			ITensor L_Ld = li_ * ljd_ * Id_ibra_jket + lj_ * lid_ * Id_iket_jbra;
 
			ITensor Dij = gamma * (L_Ld - 0.5 * LdL_I - 0.5 * I_LdL);


			vector<int> jnket = {i,j};
			vector<int> jnbra = {i,j};

			MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dij);
		
			gates.push_back(g);

		}


	}


	vector<MyBondGateDiss> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGateDiss gate : gates_) gates.push_back(gate);
	
	return gates;
}


// Dissipative gate : we have a list of tensors which act 
// We have a Lidbland of the form L_{i,j} \rho L_{i,j}^\dagger + ...

// NOT USEFUL AT THE MOMENT - NOT TESTED (SHOULD WORK)
vector<MyBondGateDiss>
gates_local_nsites_lindblad(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt)
{

	int N = length(sites);
	vector<MyBondGateDiss> gates;


	// check size of the two containers

	if(Lj_sites.size() != Lij_list.size()){
		cerr << "vectors containing jumps and sites where they act have different length!\n";
		cerr << "Lj_sites has length " << Lj_sites.size() << "\n";
		cerr << "Lj has length " << Lij_list.size() << "\n";
		cerr << "Returning empty gates\n";
		return gates;
	}


	int M = Lj_sites.size();

	for(int k=0 ; k < M; k++)
	{
		vector<int> jn = Lj_sites[k];
		if(jn.size() > 2){
			cerr << "lindlbland acting on more than two sites not yet implemented!\n Returning empty gates";
			return gates;
		} 

		// sites where it acts
		int i = jn[0];
		int j = jn[1];

		// list of jump operators

		ITensor lij  = Lij_list[k];
		ITensor lijd = dag(lij);  // it is equal to complex conjugation - it does not swap indices to make the transpose

		// site index

		Index si  = sites(i);
		Index si1 = prime(si);
		Index si2 = prime(si,2);
		Index si3 = prime(si,3);

		Index sj  = sites(j);
		Index sj1 = prime(sj);
		Index sj2 = prime(sj,2);
		Index sj3 = prime(sj,3);

		// Identity for the non-hermitian Hamiltonian part
		
		ITensor Idket = ITensor(si,sj,si2,sj2);
		ITensor Idbra = ITensor(si1,sj1,si3,sj3);
		
		for(int q=1 ; q<=dim(sj) ; q++)
		{
			Idket.set( si(q) , si2(q) , sj(q)  , sj2(q) , 1.);
			Idbra.set(si1(q) , si3(q) , sj1(q) , sj3(q) , 1.);
		}


		// v2: I think correct version

		lijd.mapPrime(0,2);
		ITensor ljdlj_I = lijd * lij * Idbra; // acts on ket
	
		lijd.mapPrime(2,0);
		lijd.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
		ITensor I_ljdlj = lijd * lij; // acts on bra  - (3,0)
		I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
		I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3)
		lijd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
		
		// start v1: I think it does not do correctly the transposition

		// jumps 
		lij.mapPrime(1,2);
		lijd.mapPrime(1,3);
		lijd.mapPrime(0,1);
		ITensor lj_ljd = lij * lijd;


		ITensor Dj = gammaj[k] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

		vector<int> jnket = {i,j};
		vector<int> jnbra = {i,j};

		MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dj);
	
		gates.push_back(g);

	}


	vector<MyBondGateDiss> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGateDiss gate : gates_) gates.push_back(gate);
	
	return gates;
}


// Dissipative impurity acting on the unfolded density matrix in an impurity problem.
// the first half sites represent the bra and evolve via -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place
// Here we apply the jump/nonunitary part on the bond in between

// There was an error - the swap done at the end was a mistake (referring to modification 03.09.2023)
// for hermitian jump it was not a problem. For non hermitian one yes.

vector<MyBondGateDiss>
gates_dissipative_impurity(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt)
{

	int N = length(sites)/2;

	vector<MyBondGateDiss> gates;

	vector<ITensor> Id;

	for(int q=N ; q<=N+1; q++)
	{
		Id.push_back(     op(sites,"Id",q) );
	}
	// the bond gate acts on site [N,N+1]
	// lj1 acts on bra
	// lj2 acts on ket
	ITensor lj1 = Lj[0];
	ITensor lj2 = Lj[1];

	ITensor lj1d = conj(lj1);
	ITensor lj2d = conj(lj2);
	
	lj1.mapPrime(0,2); // I have to do L^* L
	lj2d.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)

	ITensor ljdlj1 = lj1d * lj1; 
	ITensor ljdlj2 = lj2d * lj2; 
	
	// acting on bra (L^dag L)^T = (L^T L*)
	ljdlj1.mapPrime(2,1);
	// acting on ket (L^dag L)
	ljdlj2.mapPrime(2,1);
	
	// I am here using convention that LdL_I = (L^\dag L) \otimes I is: first operator act on ket, the second on bra
	// Here the first half describe the bra , and the second half ket. This is why I have Id[0] * ljdlj2
	ITensor LdL_I = Id[0] * ljdlj2;
	ITensor I_LdL = ljdlj1 * Id[1];

	lj1 = Lj[0];
	lj2 = Lj[1];
	lj1d = conj(lj1);


	ITensor D = gamma * (lj1d * lj2 - 0.5 * LdL_I - 0.5 * I_LdL);

	vector<int> jnket = {N,N+1};
	vector<int> jnbra = {N,N+1};

	MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,D);

	gates.push_back(g);

	return gates;
}


// Dissipative impurity acting on the unfolded density matrix in an impurity problem.
// the first half sites represent the bra and evolve via -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place
// Here we apply the jump/nonunitary part on the bond in between.

// Differences with gates_dissipative_impurity: above we used a first order Kraus approximation of the 
// dissipative part. Here, I use the class BondGate of ITensor which is able to exponentiate
// also non hermitian things since it uses a high grade Pade approximation

vector<BondGate>
gates_dissipative_impurity_high_pade(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt)
{

	int N = length(sites)/2;

	vector<BondGate> gates;

	vector<ITensor> Id;

	for(int q=N ; q<=N+1; q++)
	{
		Id.push_back(     op(sites,"Id",q) );
	}
	// the bond gate acts on site [N,N+1]
	// lj1 acts on bra
	// lj2 acts on ket
	ITensor lj1 = Lj[0];
	ITensor lj2 = Lj[1];

	ITensor lj1d = conj(lj1);
	ITensor lj2d = conj(lj2);
	
	lj1.mapPrime(0,2); // I have to do L^* L
	lj2d.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)

	ITensor ljdlj1 = lj1d * lj1; 
	ITensor ljdlj2 = lj2d * lj2; 
	
	// acting on bra (L^dag L)^T = (L^T L*)
	ljdlj1.mapPrime(2,1);
	// acting on ket (L^dag L)
	ljdlj2.mapPrime(2,1);
	
	// I am here using convention that LdL_I = (L^\dag L) \otimes I is: first operator act on ket, the second on bra
	// Here the first half describe the bra , and the second half ket. This is why I have Id[0] * ljdlj2
	ITensor LdL_I = Id[0] * ljdlj2;
	ITensor I_LdL = ljdlj1 * Id[1];


	lj1 = Lj[0];
	lj2 = Lj[1];
	lj1d = conj(lj1);

	ITensor D = gamma * (lj1d * lj2 - 0.5 * LdL_I - 0.5 * I_LdL);

	vector<int> jnket = {N,N+1};
	vector<int> jnbra = {N,N+1};

	BondGate g = BondGate(sites,N,N+1,BondGate::tImag,-1*dt,D); 
	gates.push_back(g);

	return gates;
}


// Gates of impurity problem (bath treated exactly in energy basis)
// this results in a hihgly non local Hamiltonian. Specifically
// H = 
// we do not implement the swap gates as gates, but directly in the TEBD algorithm.

// We move from the center of the chain (where the bond connecting bra and ket is located)
// moving towards the outer part

// we feed in input J_eps, hup, hdn in the right order (acting on sites from 1 to N) in the TRUE system
// wee perform first evolution of the bra [1,N]
// then, we perform the evolution of the ket [N+1,2*N]
// we start from the center, so that the first interaction is nearest-neighbor and then
// we should apply swapgates

vector<MyBondGate>
doubling_space_gates(const vector<BondGate> gates_single, const SiteSet sites_single,  const SiteSet sites_doubled)
{

    vector<MyBondGate> gates_doubled;
	int N =  length(sites_single);
	cerr << N << "\n";
	cerr << "Entering ket\n";
	// ket dynamics
    for( BondGate g : gates_single)
    {
		cerr << "Here.\n";
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
		cerr << i1 << "\n";
		cerr << i2 << "\n";
        // position along the doubled space
        int inew_1 = N + 1 + i1;
        int inew_2 = N + 1 + i2;
        
        ITensor gate 	= g.gate();
        Index si1 		= sites_single(i1);
        Index si2 		= sites_single(i2);
        Index sinew1 	= sites_doubled(inew_1);
        Index sinew2 	= sites_doubled(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));
        
        MyBondGate gnew = MyBondGate(sites_doubled,{inew_1,inew_2},0,gate);
        gnew.modify_gate(gate);
        gates_doubled.push_back(gnew);
		
    }

	cerr << "Entering bra\n";

    // bra dynamics
    for( BondGate g : gates_single)
    {
        // physical space position
        int i1 = g.i1();
        int i2 = g.i2();
		cerr << i1 << "\n";

        // position along the doubled space
        int inew_1 = N + 2 - i1;
        int inew_2 = N + 2 - i2;
        
        ITensor gate 	= g.gate();
        Index si1 		= sites_single(i1);
        Index si2 		= sites_single(i2);
        Index sinew1 	= sites_doubled(inew_1);
		Index sinew2 	= sites_doubled(inew_2); 

        // changed sites
        gate *= delta(si1,sinew1);
        gate *= delta(si2,sinew2);
        gate *= delta(prime(si1),prime(sinew1));
        gate *= delta(prime(si2),prime(sinew2));

        MyBondGate gnew = MyBondGate(sites_doubled,{inew_1,inew_2},0,gate);
        gnew.modify_gate(dag(gate));
        gates_doubled.push_back(gnew);
    }

	return gates_doubled;

}
