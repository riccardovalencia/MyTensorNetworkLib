#include "TEBD.h"
#include "observables.h"
#include "MyClasses.h"
#include "/home/ricval/Documenti/MyLibrary/TDVP/tdvp.h"
#include "/home/ricval/Documenti/MyLibrary/TDVP/basisextension.h"
#include "state_manipulation.h"
#include <itensor/all.h>
#include <iostream>
#include <math.h>       
#include <complex>
using namespace std;
using namespace itensor;




// gates Hamiltonian H = -Jxx \sum_j X_j X_{j+1} - Jyy \sum_j Y_j Y_{j+1}  - Jzz \sum_j Z_j Z_{j+1} - \sum_j \vec{h} \vec{S}_j
vector<MyBondGate>
gates_spin_model(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt)
{

	int N = length(sites);

    vector<MyBondGate> gates;

	cerr << "vector J = (J_xx, J_yy , J_zz)\n";
	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];
	cerr << Jxx << "\n" << Jyy << "\n" << Jzz << "\n";

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];
	// cerr << hx << "\n" << hy << "\n" << hz << "\n";
	// exit(0);
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

		H_SS  = Jxx * X[0] * X[1];
		H_SS += Jyy * Y[0] * Y[1];
		H_SS += Jzz * Z[0] * Z[1];
		// H_SS = Jzz * Z[0] * Z[1];
		
		
		ITensor H = H_S + H_SS;

		vector<int> jn = {j,j+1};
		// BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}
	
	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());
	for(MyBondGate gate : gates_) gates.push_back(gate);

	return gates;
}



// gates Hamiltonian H = -Jxx \sum_j X_j X_{j+1} - Jyy \sum_j Y_j Y_{j+1}  - Jzz \sum_j Z_j Z_{j+1} - \sum_j \vec{h} \vec{S}_j
// Same as the one retuning <MyBondGate>: overload of the function. Depending on whether I need flexibility and control 
// on its action, I might use one or the other
vector<BondGate>
gates_spin_model_bondgate(const SiteSet sites , const vector<double> J, const vector<double> h, const double dt)
{

	int N = length(sites);

    vector<BondGate> gates;

	cerr << "vector J = (J_xx, J_yy , J_zz)\n";
	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];
	cerr << Jxx << "\n" << Jyy << "\n" << Jzz << "\n";

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];
	// cerr << hx << "\n" << hy << "\n" << hz << "\n";
	// exit(0);
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

		H_SS  = Jxx * X[0] * X[1];
		H_SS += Jyy * Y[0] * Y[1];
		H_SS += Jzz * Z[0] * Z[1];
		// H_SS = Jzz * Z[0] * Z[1];
		
		
		ITensor H = H_S + H_SS;

		vector<int> jn = {j,j+1};
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
		// MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}
	
	vector<BondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());
	for(BondGate gate : gates_) gates.push_back(gate);

	return gates;
}

// effective Hamiltonian of a local dissipative process
// The coherent dynamics is given by the generic short range spin model

vector<BondGate>
gates_spin_eff_model_bondgate(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const vector<int> Lj_sites, const vector<double> gamma, const double dt)
{

	int N = length(sites);

    vector<BondGate> gates;

	// cerr << "vector J = (J_xx, J_yy , J_zz)\n";
	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];
	// cerr << Jxx << "\n" << Jyy << "\n" << Jzz << "\n";

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];
	// cerr << hx << "\n" << hy << "\n" << hz << "\n";
	// exit(0);
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

		for(int k=0 ; k<= Lj_sites.size(); k++)
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
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
		// MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}
	
	vector<BondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());
	for(BondGate gate : gates_) gates.push_back(gate);

	return gates;
}

// // WORKING - THE PREVIOUS VERSION (LIKELY) HAS ISSUES WITH NON HERMITIAN JUMPS

// vector<MyBondGateDiss>
// gates_dissipative_impurity_v2(const SiteSet sites , const vector<ITensor> Lj, const double gamma, const double dt)
// {

// 	int N = length(sites)/2;

// 	vector<MyBondGateDiss> gates;

// 	vector<ITensor> Id;

// 	for(int q=N ; q<=N+1; q++)
// 	{
// 		Id.push_back(     op(sites,"Id",q) );
// 	}
// 	// the bond gate acts on site [N,N+1]
// 	// lj1 acts on bra
// 	// lj2 acts on ket
// 	ITensor lj1 = Lj[0];
// 	ITensor lj2 = Lj[1]; 

// 	ITensor lj1d = conj(lj1);
// 	ITensor lj2d = conj(lj2);
	
// 	lj1.mapPrime(0,2); // I have to do L^T L^*
// 	lj2d.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)

// 	ITensor ljdlj1 = lj1 * lj1d; 
// 	ITensor ljdlj2 = lj2d * lj2; 
	
// 	// acting on bra (L^dag L)^T = (L^T L*)
// 	ljdlj1.mapPrime(2,1);
// 	// acting on ket (L^dag L)
// 	ljdlj2.mapPrime(2,1);
// 	// ljdlj2.swapPrime(1,0);

// 	// I am here using convention that LdL_I = (L^\dag L) \otimes I is: first operator act on ket, the second on bra
// 	// Here the first half describe the bra , and the second half ket. This is why I have Id[0] * ljdlj2
// 	ITensor LdL_I = Id[0] * ljdlj2;
// 	ITensor I_LdL = ljdlj1 * Id[1];

// 	lj1 = Lj[0];
// 	lj2 = Lj[1];
// 	lj1d = conj(lj1);

// 	// test 03.09.2023 (it is not relevant for hermitian jumps)
// 	// lj1d = swapPrime(lj1d,0,1);

// 	ITensor D = gamma * (lj1d * lj2 - 0.5 * LdL_I - 0.5 * I_LdL);

// 	vector<int> jnket = {N,N+1};
// 	vector<int> jnbra = {N,N+1};

// 	MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,D);

// 	gates.push_back(g);

// 	return gates;
// }



// Dissipative impurity acting on the unfolded density matrix in an impurity problem.
// the first half sites represent the bra and evolve via -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place
// Here we apply the jump/nonunitary part on the bond in between

// There was an error - the swap done at the end was a mistake (referring to modification 03.09.2023 )
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

	// test 03.09.2023 (it is not relevant for symmetric jumps) - WAS A MISTAKE!
	// lj1d = swapPrime(lj1d,0,1);


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
// Here we apply the jump/nonunitary part on the bond in between
// DIFFERENENCES WITH ABOVE: above we used a first order Kraus approximation of the 
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

	// cerr << "L* : \n";
	// PrintData(lj1d);
	// cerr << "L : \n";
	// PrintData(lj2);
	// cerr << "L^d L : \n";
	// PrintData(ljdlj2);
	// cerr << "L^T L* : \n";
	// PrintData(ljdlj1);
	// exit(0);

	ITensor D = gamma * (lj1d * lj2 - 0.5 * LdL_I - 0.5 * I_LdL);

	vector<int> jnket = {N,N+1};
	vector<int> jnbra = {N,N+1};

	BondGate g = BondGate(sites,N,N+1,BondGate::tImag,-1*dt,D); 
	gates.push_back(g);

	return gates;
}

// Simulation of an impurity problem located on the first physical site.
// Can be applied to the unfolded density matrix: 
// the first half sites represent the bra and evolve via  -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place
// This part takes care of the coherent dynamics - we have hopping of free spinful fermions

vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{ 

	// careful with local fields if not homogeneous - I have to flip the array. Remember that [1-N] corresponds to bra 
	vector<double> Jket = J;
	vector<double> Jbra = J;

	vector<double> hupket = hup;
	vector<double> hupbra = hup;

	vector<double> hdnket = hdn;
	vector<double> hdnbra = hdn;

	reverse(Jbra.begin()  , Jbra.end());
	reverse(hupbra.begin(), hupbra.end());
	reverse(hdnbra.begin(), hdnbra.end());

	// hopping in the doubled space
	vector<double> J2;
	J2.insert( J2.end(), Jbra.begin(), Jbra.end() );
	// acting on the bond connecting bra - ket (so it should be 0)
	J2.push_back(0.);
	J2.insert( J2.end(), Jket.begin(), Jket.end() );

	// local fields in the doubled space

	vector<double> hup2;
	vector<double> hdn2;

	hup2.insert( hup2.end(), hupbra.begin(), hupbra.end() );
	hup2.insert( hup2.end(), hupket.begin(), hupket.end() );
	hdn2.insert( hdn2.end(), hdnbra.begin(), hdnbra.end() );
	hdn2.insert( hdn2.end(), hdnket.begin(), hdnket.end() );

	// for(double j : hup2) cerr << j << " ";
	// cerr << "\n";
	// for(double j : hdn2) cerr << j << " ";
	// exit(0);



	// N is always even and corresponds to the physical size
	int N2 = length(sites);
	int N = N2/2;

    vector<BondGate> gates;
	vector<BondGate> gates_ket;
	vector<BondGate> gates_bra;
	

	for(int j=1 ; j <= N2-1 ; j+=1)
	{

		// act on bra j \in [1,N] or ket space j \in [N+1,2*N]
		if(j!=N)
		{
			cerr << "Site : " << j << "\n";

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

				// ITensor o = prime(op(sites,"Adagup",q)) * op(sites,"F",q);
				// o.mapPrime(2,1);
				// AdagupFi.push_back(o);

				// o = prime(op(sites,"Aup",q)) * op(sites,"F",q) ;
				// o.mapPrime(2,1);
				// AupFi.push_back(o);

				// o = prime(op(sites,"F",q)) * op(sites,"Adn",q);
				// o.mapPrime(2,1);
				// FiAdn.push_back(o);

				// o = prime(op(sites,"F",q)) * op(sites,"Adagdn",q);
				// o.mapPrime(2,1);
				// FiAdagdn.push_back(o);

			
				AdagupFi.push_back(op(sites,"Adagup*F",q));
				AupFi.push_back(op(sites,"Aup*F",q));
				FiAdn.push_back(op(sites,"F*Adn",q));
				FiAdagdn.push_back(op(sites,"F*Adagdn",q));


			}

			ITensor H, H_S , H_SS;

			// build single site
		
			if(j == 1 || j == (N+1) ) H_S = (hup2[j-1] * Nup[0] + hdn2[j-1] * Ndn[0]) * Id[1] ;
			else                      H_S = (hup2[j-1] * Nup[0] + hdn2[j-1] * Ndn[0]) * Id[1] / 2.;

			if(j == N-1 || j == 2*N-1) H_S += Id[0] * (hup2[j] * Nup[1] + hdn2[j] * Ndn[1]) ;
			else                       H_S += Id[0] * (hup2[j] * Nup[1] + hdn2[j] * Ndn[1]) / 2.      ;

			cerr << "Site : " << j     << " -> h : " << hup2[j-1] << "\n";
			cerr << "Site : " << j + 1 << " -> h : " << hup2[j] << "\n";

			// two sites hopping
			
			if(j>N)
			{
				H_SS   = J2[j-1] * (  AdagupFi[0] * Aup[1]   - AupFi[0] * Adagup[1] );
				H_SS  += J2[j-1] * (  Adagdn[0]   * FiAdn[1] - Adn[0] * FiAdagdn[1] );
			}
			else
			{
				// the ordering of the physical indices is the opposite with respect to the indices
				// used here (the bra is inverted)
				H_SS   = J2[j-1] * (  Aup[0] * AdagupFi[1] - Adagup[0] * AupFi[1] );
				H_SS  += J2[j-1] * (  FiAdn[0] * Adagdn[1] - FiAdagdn[0]  * Adn[1] );
			}
			
			
			cerr << "Site : " << j << " -> " <<  J2[j-1] << endl;
			if(j<N){
				// action on bra (it is exp(+i H^T t))
				// attempt - the sign is due to the convention of the sites 
				H = - H_S - H_SS ;
				H = swapPrime(H,0,1); // NOT SURE
				BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
				gates_bra.push_back(g);
			}
			// action on ket
			else
			{
				H = H_S + H_SS ;
				BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
				gates_ket.push_back(g);
			}
	
		}
		
	}
	

	vector<BondGate> gates_ket_reversed = gates_ket;
	vector<BondGate> gates_bra_reversed = gates_bra;
	
	reverse(gates_ket_reversed.begin(), gates_ket_reversed.end());
	reverse(gates_bra_reversed.begin(), gates_bra_reversed.end());

	for(BondGate g : gates_bra)  		 gates.push_back(g);
	for(BondGate g : gates_ket) 		 gates.push_back(g);
	for(BondGate g : gates_ket_reversed) gates.push_back(g);
	for(BondGate g : gates_bra_reversed) gates.push_back(g);

	return gates;
}


// Simulation of free spinful fermions 

vector<BondGate>
gates_free_spinful_fermions(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{ 


	int N = length(sites);

    vector<BondGate> gates;
	
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
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
		gates.push_back(g);
			
	}


	vector<BondGate> gates_reversed = gates;
	
	reverse(gates_reversed.begin(), gates_reversed.end());

	for(BondGate g : gates_reversed) gates.push_back(g);

	return gates;
}




// Simulation of an impurity problem located on the first physical site.
// Can be applied to the unfolded density matrix: 
// the first half sites represent the bra and evolve via -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place


vector<BondGate>
gates_coherent_part_spin_dissipative_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> h, const vector<ITensor> Lj, const double gamma, const double dt)
{

	// N is always even and corresponds to the physical size
	int N2 = length(sites);
	int N = N2/2;


    vector<BondGate> gates;

	// cerr << "vector J = (J_xx, J_yy , J_zz)\n";
	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];
	// cerr << Jxx << "\n" << Jyy << "\n" << Jzz << "\n";

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];
	// cerr << hx << "\n" << hy << "\n" << hz << "\n";
	// exit(0);
	for(int j=1 ; j <= N2-1 ; j+=1)
	{
		cerr << "Site : " << j << "\n";

		vector<ITensor> X;
		vector<ITensor> Y;
		vector<ITensor> Z;
		vector<ITensor> Id;

		for(int q=j ; q<=j+1; q++)
		{
			Id.push_back(     op(sites,"Id",q) );
			X.push_back(  2 * op(sites,"Sx",q) );
			Y.push_back(  2 * op(sites,"Sy",q) );
			Z.push_back(  2 * op(sites,"Sz",q) );
		}

		ITensor H, H_S , H_SS;

		// act on bra j \in [1,N] or ket space j \in [N+1,2*N]
		if(j!=N)
		{
			if(j % N ==1) H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] ;
			else          H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] / 2.;

			if(j == N-1 || j == 2*N-1) H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) ;
			else                       H_S  += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) / 2.      ;

			H_SS  = Jxx * X[0] * X[1] + Jyy * Y[0] * Y[1] + Jzz * Z[0] * Z[1];

			// action on bra (it is exp(+i H^T t))
			if(j<N){
				H = -H_S - H_SS;
				H = swapPrime(H,0,1);
			}
			// action on ket
			else    H =  H_S + H_SS;
		}

		// action on the bond connecting bra-ket
		// non hermitian dynamics
		// else
		// {
		// 	// lj1 acts on bra
		// 	// lj2 acts on ket
		// 	ITensor lj1 = Lj[0];
		// 	ITensor lj2 = Lj[1];

		// 	ITensor lj1d = conj(lj1);
		// 	ITensor lj2d = conj(lj2);
			
		// 	lj1.mapPrime(0,2); // I have to do L^* L
		// 	lj2d.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)

		// 	ITensor ljdlj1 = lj1d * lj1; 
		// 	ITensor ljdlj2 = lj2d * lj2; 
			
		// 	// acting on bra
		// 	ljdlj1.mapPrime(2,1);
		// 	// acting on ket
		// 	ljdlj2.mapPrime(2,1);
			
		// 	// non hermitian part
		// 	// non hermitian part
		// 	// H_S -= 0.5 * gamma * Cplx_i * ljdlj1 * Id[1]; 
		// 	// H_S -= 0.5 * gamma * Cplx_i * Id[0] * ljdlj2;

		// 	H_S = - 0.5 * gamma * ( ljdlj1 * Id[1]  + Id[0] * ljdlj2); 


		// 	lj1 = Lj[0];
		// 	lj2 = Lj[1];
		// 	lj1d = conj(lj1);
		// 	H_SS = gamma * lj1d * lj2; 
		// 	H = H_S + H_SS;	
		
		// }

		if(j !=N)
		{
			BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
			gates.push_back(g);
		}
		
		
	}

	for(int j=N2-1 ; j >= 1; j-=1)
	{
		cerr << "Site : " << j << "\n";

		vector<ITensor> X;
		vector<ITensor> Y;
		vector<ITensor> Z;
		vector<ITensor> Id;

		for(int q=j ; q<=j+1; q++)
		{
			Id.push_back(     op(sites,"Id",q) );
			X.push_back(  2 * op(sites,"Sx",q) );
			Y.push_back(  2 * op(sites,"Sy",q) );
			Z.push_back(  2 * op(sites,"Sz",q) );
		}

		ITensor H, H_S , H_SS;

		// act on bra j \in [1,N] or ket space j \in [N+1,2*N]
		if(j!=N)
		{
			if(j ==1 || j==N+1) H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] ;
			else       H_S = (hx * X[0] + hy * Y[0] + hz * Z[0]) * Id[1] / 2.;

			if(j == N-1 || j == 2*N-1) H_S += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) ;
			else           H_S   += Id[0] * (hx * X[1] + hy * Y[1] + hz * Z[1]) / 2.      ;

			H_SS  = Jxx * X[0] * X[1] + Jyy * Y[0] * Y[1] + Jzz * Z[0] * Z[1];
	

			// action on bra
			if(j<N){
				H = -H_S - H_SS;
				H = swapPrime(H,0,1);
			}
			// action on ket
			else    H =  H_S + H_SS;
		}

		// action on the bond connecting bra-ket
		// non hermitian dynamics
		// else
		// {
		// 	// lj1 acts on bra
		// 	// lj2 acts on ket
		// 	ITensor lj1 = Lj[0];
		// 	ITensor lj2 = Lj[1];

		// 	ITensor lj1d = conj(lj1);
		// 	ITensor lj2d = conj(lj2);
			
		// 	lj1.mapPrime(0,2); // I have to do L^* L
		// 	lj2d.mapPrime(0,2); // I have to do L^dag L, which is (L^*)^T L (this is why I make the 'row' index the 'column' one)

		// 	ITensor ljdlj1 = lj1d * lj1; 
		// 	ITensor ljdlj2 = lj2d * lj2; 
			
		// 	// acting on bra
		// 	ljdlj1.mapPrime(2,1);
			
		// 	// acting on ket
		// 	ljdlj2.mapPrime(2,1);
			
		// 	// non hermitian part
		// 	// H_S -= 0.5 * gamma * Cplx_i * ljdlj1 * Id[1]; 
		// 	// H_S -= 0.5 * gamma * Cplx_i * Id[0] * ljdlj2;

		// 	H_S = - 0.5 * gamma * ( ljdlj1 * Id[1]  + Id[0] * ljdlj2); 

		// 	lj1 = Lj[0];
		// 	lj2 = Lj[1];
		// 	lj1d = conj(lj1);
		// 	H_SS = gamma * lj1d * lj2; 
		// 	H = H_S + H_SS;	
		
		// }

		
		if(j !=N)
		{
			BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,H); 
			gates.push_back(g);
		}
		
	}


	
	// vector<BondGate> gates_ = gates;
	// reverse(gates_.begin(), gates_.end());
	// for(BondGate gate : gates_) gates.push_back(gate);

	return gates;
}



// Simulation of an impurity problem located on the first physical site.
// Can be applied to the unfolded density matrix: 
// the first half sites represent the bra and evolve via -H
// the second half sites represent the ket and evolve via +H
// the bond in between site N and N+1 is where jump/nonunitary dynamics take place


vector<MyBondGate>
gates_coherent_part_spin_dissipative_NNN_interactions_impurity_model(const SiteSet sites , const vector<double> J, const vector<double> J_NNN, const vector<double> h, const double dt)
{

	// N is always even and corresponds to the physical size
	int N2 = length(sites);
	int N = N2/2;


    vector<MyBondGate> gates, gates_bra, gates_ket;

	double Jxx = J[0];
	double Jyy = J[1];
	double Jzz = J[2];

	double hx  = h[0];
	double hy  = h[1];
	double hz  = h[2];


	vector<double> hx_tot;
	vector<double> Jzz_tot;
	for(int j : range1(N))   hx_tot.push_back(0);
	for(int j : range1(N-1)) Jzz_tot.push_back(0);

	if(N<3)
	{
		cerr << "Less than 3-sites. Not possible to implement 3-site TEBD.\n";
		exit(0);
	}

	else if(N==3)
	{
		for(int j = 1; j < N2 ; j+=3)
		{
			ITensor H, H_S , H_SS, H_SSS;

			int js = j;
			int jf = js+2;
				
			vector<double> hx_vec = { hx, hx, hx};
			vector<double> hy_vec = { hy, hy, hy};
			vector<double> hz_vec = { hz, hz, hz};
			
			// nearest-neighbor interactions - when acting in the bulk
			vector<double> Jxx_vec = {Jxx,Jxx};
			vector<double> Jyy_vec = {Jyy,Jyy};
			vector<double> Jzz_vec = {Jzz,Jzz};

			vector<ITensor> Xj;
			vector<ITensor> Yj;
			vector<ITensor> Zj;
			vector<ITensor> Ij;

			for(int q=js ; q<= jf; q++)
			{
				Ij.push_back(      op(sites,"Id",q) );
				Xj.push_back(  2 * op(sites,"Sx",q) );
				Yj.push_back(  2 * op(sites,"Sy",q) );
				Zj.push_back(  2 * op(sites,"Sz",q) );
			}

			H_S   = hx_vec[0] * Xj[0] * Ij[1] * Ij[2];
			H_S  += hx_vec[1] * Ij[0] * Xj[1] * Ij[2];
			H_S  += hx_vec[2] * Ij[0] * Ij[1] * Xj[2];

			H_S  += hy_vec[0] * Yj[0] * Ij[1] * Ij[2];
			H_S  += hy_vec[1] * Ij[0] * Yj[1] * Ij[2];
			H_S  += hy_vec[2] * Ij[0] * Ij[1] * Yj[2];

			H_S  += hz_vec[0] * Zj[0] * Ij[1] * Ij[2];
			H_S  += hz_vec[1] * Ij[0] * Zj[1] * Ij[2];
			H_S  += hz_vec[2] * Ij[0] * Ij[1] * Zj[2];

			H_SS   = Jxx_vec[0] * Xj[0] * Xj[1] * Ij[2] ;
			H_SS  += Jxx_vec[1] * Ij[0] * Xj[1] * Xj[2] ;

			H_SS  += Jyy_vec[0] * Yj[0] * Yj[1] * Ij[2] ;
			H_SS  += Jyy_vec[1] * Ij[0] * Yj[1] * Yj[2] ;

			H_SS  += Jzz_vec[0] * Zj[0] * Zj[1] * Ij[2] ;
			H_SS  += Jzz_vec[1] * Ij[0] * Zj[1] * Zj[2] ;

			H_SSS  = J_NNN[0] * Xj[0] * Ij[1] * Xj[2] ;
			H_SSS += J_NNN[1] * Yj[0] * Ij[1] * Yj[2] ;
			H_SSS += J_NNN[2] * Zj[0] * Ij[1] * Zj[2] ;
			
	
			vector<int> jn = {js,js+1,js+2};
			if(j<N){
				// action on bra (it is exp(+i H^T t))

				H = -H_S - H_SS - H_SSS;
				H = swapPrime(H,0,1);

				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_bra.push_back(g);
			}
			// action on ket
			else
			{
				H = H_S + H_SS + H_SSS;
				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_ket.push_back(g);

			}
		}
	}

	else
	{
		for(int layer = 1 ; layer <= 3 ; layer ++)
		{
			for(int j = layer ; j <= N-2 ; j+=3)
			{
			
				cerr << "Layer : " << layer << " site : " << j << "\n";

				ITensor H, H_S , H_SS, H_SSS;

				int js = j;
				int jf = js + 2;
				
				vector<double> hx_vec,   hy_vec,  hz_vec; 
				// nearest-neighbor interactions - when acting in the bulk
				vector<double> Jxx_vec = {Jxx/2.,Jxx/2.};
				vector<double> Jyy_vec = {Jyy/2.,Jyy/2.};
				vector<double> Jzz_vec = {Jzz/2.,Jzz/2.};

				vector<ITensor> Xj;
				vector<ITensor> Yj;
				vector<ITensor> Zj;
				vector<ITensor> Ij;

				for(int q=js ; q<= jf; q++)
				{
					Ij.push_back(     op(sites,"Id",q) );
					Xj.push_back(  2 * op(sites,"Sx",q) );
					Yj.push_back(  2 * op(sites,"Sy",q) );
					Zj.push_back(  2 * op(sites,"Sz",q) );
				}

				
				if( js % N == 1 )
				{
					hx_vec = {hx, hx/2. , hx/3.};
					hy_vec = {hy, hy/2. , hy/3.};
					hz_vec = {hz, hz/2. , hz/3.};
				}
				else if( js % N == 2 )
				{
					hx_vec = {hx/2. , hx/3. , hx/3.};
					hy_vec = {hy/2. , hy/3. , hy/3.};
					hz_vec = {hz/2. , hz/3. , hz/3.};
				}
				else if(jf % (N-1) == 0)
				{
					hx_vec = {hx/3. , hx/3. , hx/2.};
					hy_vec = {hy/3. , hy/3. , hy/2.};
					hz_vec = {hz/3. , hz/3. , hz/2.};
				}
				else if(jf % N == 0)
				{
					hx_vec = {hx/3. , hx/2. , hx};
					hy_vec = {hy/3. , hy/2. , hy};
					hz_vec = {hz/3. , hz/2. , hz};
				}
				else
				{
					hx_vec = {hx/3. , hx/3. , hx/3.};
					hy_vec = {hy/3. , hy/3. , hy/3.};
					hz_vec = {hz/3. , hz/3. , hz/3.};
				}

				if( js%N == 1) // First physical site (it technically correspond to the last site of the bra, but since couplings are homogeneous it does not matter)
				{
					Jxx_vec[0] = Jxx;
					Jyy_vec[0] = Jyy;
					Jzz_vec[0] = Jzz;
				}

				if( jf%N == 0) // Last physical site
				{
					Jxx_vec[1] = Jxx;
					Jyy_vec[1] = Jyy;
					Jzz_vec[1] = Jzz;
				}

				// for(double jj : Jxx_vec) cerr << jj << " ";
				// cerr << "\n";
				// for(double jj : Jyy_vec) cerr << jj << " ";
				// cerr << "\n";
				// for(double jj : Jzz_vec) cerr << jj << " ";
				// cerr << "\n";
				for(double hj : hx_vec) cerr << hj << " ";
				cerr << "\n";


				Jzz_tot[j-1] += Jzz_vec[0];
				Jzz_tot[j]   += Jzz_vec[1];

				// hx_tot[j-1]   += hx_vec[0];
				// hx_tot[j] += hx_vec[1];
				// hx_tot[j+1] += hx_vec[2];
				
				// for(double hj : hy_vec) cerr << hj << " ";
				// cerr << "\n";
				// for(double hj : hz_vec) cerr << hj << " ";
				// cerr << "\n";
				
				
				

				H_S   = hx_vec[0] * Xj[0] * Ij[1] * Ij[2];
				H_S  += hx_vec[1] * Ij[0] * Xj[1] * Ij[2];
				H_S  += hx_vec[2] * Ij[0] * Ij[1] * Xj[2];

				H_S  += hy_vec[0] * Yj[0] * Ij[1] * Ij[2];
				H_S  += hy_vec[1] * Ij[0] * Yj[1] * Ij[2];
				H_S  += hy_vec[2] * Ij[0] * Ij[1] * Yj[2];

				H_S  += hz_vec[0] * Zj[0] * Ij[1] * Ij[2];
				H_S  += hz_vec[1] * Ij[0] * Zj[1] * Ij[2];
				H_S  += hz_vec[2] * Ij[0] * Ij[1] * Zj[2];

				H_SS   = Jxx_vec[0] * Xj[0] * Xj[1] * Ij[2] ;
				H_SS  += Jxx_vec[1] * Ij[0] * Xj[1] * Xj[2] ;

				H_SS  += Jyy_vec[0] * Yj[0] * Yj[1] * Ij[2] ;
				H_SS  += Jyy_vec[1] * Ij[0] * Yj[1] * Yj[2] ;

				H_SS  += Jzz_vec[0] * Zj[0] * Zj[1] * Ij[2] ;
				H_SS  += Jzz_vec[1] * Ij[0] * Zj[1] * Zj[2] ;

				H_SSS  = J_NNN[0] * Xj[0] * Ij[1] * Xj[2] ;
				H_SSS += J_NNN[1] * Yj[0] * Ij[1] * Yj[2] ;
				H_SSS += J_NNN[2] * Zj[0] * Ij[1] * Zj[2] ;
				
		
				vector<int> jn = {js,js+1,js+2};
				if(j<N){
					// action on bra (it is exp(+i H^T t))

					H = -H_S - H_SS - H_SSS;
					H = swapPrime(H,0,1);

					MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
					gates_bra.push_back(g);
				}
				// action on ket
				else
				{
					H = H_S + H_SS + H_SSS;
					MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
					gates_ket.push_back(g);

				}
			}

			// acting on ket
			for(int j = N + layer ; j <= N2-2 ; j+=3)
			{
			
				cerr << "Layer : " << layer << " site : " << j << "\n";

				ITensor H, H_S , H_SS, H_SSS;

				int js = j;
				int jf = js + 2;
				
				vector<double> hx_vec,   hy_vec,  hz_vec; 
				// nearest-neighbor interactions - when acting in the bulk
				vector<double> Jxx_vec = {Jxx/2.,Jxx/2.};
				vector<double> Jyy_vec = {Jyy/2.,Jyy/2.};
				vector<double> Jzz_vec = {Jzz/2.,Jzz/2.};

				vector<ITensor> Xj;
				vector<ITensor> Yj;
				vector<ITensor> Zj;
				vector<ITensor> Ij;

				for(int q=js ; q<= jf; q++)
				{
					Ij.push_back(     op(sites,"Id",q) );
					Xj.push_back(  2 * op(sites,"Sx",q) );
					Yj.push_back(  2 * op(sites,"Sy",q) );
					Zj.push_back(  2 * op(sites,"Sz",q) );
				}

				

				if( js == N + 1 )
				{
					hx_vec = {hx, hx/2. , hx/3.};
					hy_vec = {hy, hy/2. , hy/3.};
					hz_vec = {hz, hz/2. , hz/3.};
				}
				else if( js == N+2 )
				{
					hx_vec = {hx/2. , hx/3. , hx/3.};
					hy_vec = {hy/2. , hy/3. , hy/3.};
					hz_vec = {hz/2. , hz/3. , hz/3.};
				}
				else if(jf == N2 - 1)
				{
					hx_vec = {hx/3. , hx/3. , hx/2.};
					hy_vec = {hy/3. , hy/3. , hy/2.};
					hz_vec = {hz/3. , hz/3. , hz/2.};
				}
				else if(jf == N2)
				{
					hx_vec = {hx/3. , hx/2. , hx};
					hy_vec = {hy/3. , hy/2. , hy};
					hz_vec = {hz/3. , hz/2. , hz};
				}
				else
				{
					hx_vec = {hx/3. , hx/3. , hx/3.};
					hy_vec = {hy/3. , hy/3. , hy/3.};
					hz_vec = {hz/3. , hz/3. , hz/3.};
				}

				if( js == N+1) // First physical site (it technically correspond to the last site of the bra, but since couplings are homogeneous it does not matter)
				{
					Jxx_vec[0] = Jxx;
					Jyy_vec[0] = Jyy;
					Jzz_vec[0] = Jzz;
				}

				if( jf == N2) // Last physical site
				{
					Jxx_vec[1] = Jxx;
					Jyy_vec[1] = Jyy;
					Jzz_vec[1] = Jzz;
				}

				// for(double jj : Jxx_vec) cerr << jj << " ";
				// cerr << "\n";
				// for(double jj : Jyy_vec) cerr << jj << " ";
				// cerr << "\n";
				// for(double jj : Jzz_vec) cerr << jj << " ";
				// cerr << "\n";
				for(double hj : hx_vec) cerr << hj << " ";
				cerr << "\n";
				// for(double hj : hy_vec) cerr << hj << " ";
				// cerr << "\n";
				// for(double hj : hz_vec) cerr << hj << " ";
				// cerr << "\n";
				
				hx_tot[j%N-1]   += hx_vec[0];
				hx_tot[j%N] += hx_vec[1];
				hx_tot[j%N+1] += hx_vec[2];

				// Jzz_tot[j%N-1] += Jzz_vec[0];
				// Jzz_tot[j%N]   += Jzz_vec[1];
						
				
				

				H_S   = hx_vec[0] * Xj[0] * Ij[1] * Ij[2];
				H_S  += hx_vec[1] * Ij[0] * Xj[1] * Ij[2];
				H_S  += hx_vec[2] * Ij[0] * Ij[1] * Xj[2];

				H_S  += hy_vec[0] * Yj[0] * Ij[1] * Ij[2];
				H_S  += hy_vec[1] * Ij[0] * Yj[1] * Ij[2];
				H_S  += hy_vec[2] * Ij[0] * Ij[1] * Yj[2];

				H_S  += hz_vec[0] * Zj[0] * Ij[1] * Ij[2];
				H_S  += hz_vec[1] * Ij[0] * Zj[1] * Ij[2];
				H_S  += hz_vec[2] * Ij[0] * Ij[1] * Zj[2];

				H_SS   = Jxx_vec[0] * Xj[0] * Xj[1] * Ij[2] ;
				H_SS  += Jxx_vec[1] * Ij[0] * Xj[1] * Xj[2] ;

				H_SS  += Jyy_vec[0] * Yj[0] * Yj[1] * Ij[2] ;
				H_SS  += Jyy_vec[1] * Ij[0] * Yj[1] * Yj[2] ;

				H_SS  += Jzz_vec[0] * Zj[0] * Zj[1] * Ij[2] ;
				H_SS  += Jzz_vec[1] * Ij[0] * Zj[1] * Zj[2] ;

				H_SSS  = J_NNN[0] * Xj[0] * Ij[1] * Xj[2] ;
				H_SSS += J_NNN[1] * Yj[0] * Ij[1] * Yj[2] ;
				H_SSS += J_NNN[2] * Zj[0] * Ij[1] * Zj[2] ;
				
		
				vector<int> jn = {js,js+1,js+2};
				if(j<N)
				{
					// action on bra (it is exp(+i H^T t))
					H = -H_S - H_SS - H_SSS;
					H = swapPrime(H,0,1);

					MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
					gates_bra.push_back(g);
				}
				// action on ket
				else
				{
					H = H_S + H_SS + H_SSS;
					MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
					gates_ket.push_back(g);
				}
			}
				
		}
	}
	// exit(0);
	/*
	// second layer
	for(int j=2 ; j <= N2-1 ; j+=3)
	{
		if(j!=N)
		{
			ITensor H, H_S , H_SS, H_SSS;

			int js = j;
			int jf = js + 2;
			
			vector<double> hx_vec,   hy_vec,  hz_vec; 
			// nearest-neighbor interactions - when acting in the bulk
			vector<double> Jxx_vec = {Jxx/2.,Jxx/2.};
			vector<double> Jyy_vec = {Jyy/2.,Jyy/2.};
			vector<double> Jzz_vec = {Jzz/2.,Jzz/2.};

			vector<ITensor> X;
			vector<ITensor> Y;
			vector<ITensor> Z;
			vector<ITensor> Id;

			for(int q=j ; q<=j+2; q++)
			{
				Id.push_back(     op(sites,"Id",q) );
				X.push_back(  2 * op(sites,"Sx",q) );
				Y.push_back(  2 * op(sites,"Sy",q) );
				Z.push_back(  2 * op(sites,"Sz",q) );
			}

			if( js % N == 1 )
			{
				hx_vec = {hx, hx/2. , hx/3.};
				hy_vec = {hy, hy/2. , hy/3.};
				hz_vec = {hz, hz/2. , hz/3.};
			}
			else if(  js % N == 2 )
			{
				hx_vec = {hx/2. , hx/3. , hx/3.};
				hy_vec = {hy/2. , hy/3. , hy/3.};
				hz_vec = {hz/2. , hz/3. , hz/3.};
			}
			else if(jf % (N-1) == 0)
			{
				hx_vec = {hx/3. , hx/3. , hx/2.};
				hy_vec = {hy/3. , hy/3. , hy/2.};
				hz_vec = {hz/3. , hz/3. , hz/2.};
			}
			else if(jf % N == 0)
			{
				hx_vec = {hx/3. , hx/2. , hx};
				hy_vec = {hy/3. , hy/2. , hy};
				hz_vec = {hz/3. , hz/2. , hz};
			}
			else
			{
				hx_vec = {hx/3. , hx/3. , hx/3.};
				hy_vec = {hy/3. , hy/3. , hy/3.};
				hz_vec = {hz/3. , hz/3. , hz/3.};
			}

			if( js%N == 1) // First physical site (it technically correspond to the last site of the bra, but since couplings are homogeneous it does not matter)
			{
				Jxx_vec[0] = Jxx;
				Jyy_vec[0] = Jyy;
				Jzz_vec[0] = Jzz;
			}

			if( jf%N == 0) // Last physical site
			{
				Jxx_vec[1] = Jxx;
				Jyy_vec[1] = Jyy;
				Jzz_vec[1] = Jzz;
			}

			H_S   = hx_vec[0] * Xj[0] * Ij[1] * Ij[2];
			H_S  += hx_vec[1] * Ij[0] * Xj[1] * Ij[2];
			H_S  += hx_vec[2] * Ij[0] * Ij[1] * Xj[2];

			H_S  += hy_vec[0] * Yj[0] * Ij[1] * Ij[2];
			H_S  += hy_vec[1] * Ij[0] * Yj[1] * Ij[2];
			H_S  += hy_vec[2] * Ij[0] * Ij[1] * Yj[2];

			H_S  += hz_vec[0] * Zj[0] * Ij[1] * Ij[2];
			H_S  += hz_vec[1] * Ij[0] * Zj[1] * Ij[2];
			H_S  += hz_vec[2] * Ij[0] * Ij[1] * Zj[2];

			H_SS   = Jxx_vec[0] * Xj[0] * Xj[1] * Ij[2] ;
			H_SS  += Jxx_vec[1] * Ij[0] * Xj[1] * Xj[2] ;

			H_SS  += Jyy_vec[0] * Yj[0] * Yj[1] * Ij[2] ;
			H_SS  += Jyy_vec[1] * Ij[0] * Yj[1] * Yj[2] ;

			H_SS  += Jzz_vec[0] * Zj[0] * Zj[1] * Ij[2] ;
			H_SS  += Jzz_vec[1] * Ij[0] * Zj[1] * Zj[2] ;

			H_SSS  = J_NNN[0] * Xj[0] * Ij[1] * Xj[2] ;
			H_SSS += J_NNN[1] * Yj[0] * Ij[1] * Yj[2] ;
			H_SSS += J_NNN[2] * Zj[0] * Ij[1] * Zj[2] ;
			
	
			vector<int> jn = {j,j+1,j+2};
			if(j<N){
				// action on bra (it is exp(+i H^T t))

				H = -H_S - H_SS - H_SSS;
				H = swapPrime(H,0,1);

				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_bra.push_back(g);
			}
			// action on ket
			else
			{
				H = H_S + H_SS + H_SSS;
				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_ket.push_back(g);

			}
		}
		
		
	}

	// third layer
	for(int j=3 ; j <= N2-1 ; j+=3)
	{
		if(j!=N)
		{
			ITensor H, H_S , H_SS, H_SSS;

			int js = j;
			int jf = js + 2;
			
			vector<double> hx_vec,   hy_vec,  hz_vec; 
			// nearest-neighbor interactions - when acting in the bulk
			vector<double> Jxx_vec = {Jxx/2.,Jxx/2.};
			vector<double> Jyy_vec = {Jyy/2.,Jyy/2.};
			vector<double> Jzz_vec = {Jzz/2.,Jzz/2.};

			vector<ITensor> X;
			vector<ITensor> Y;
			vector<ITensor> Z;
			vector<ITensor> Id;

			for(int q=j ; q<=j+2; q++)
			{
				Id.push_back(     op(sites,"Id",q) );
				X.push_back(  2 * op(sites,"Sx",q) );
				Y.push_back(  2 * op(sites,"Sy",q) );
				Z.push_back(  2 * op(sites,"Sz",q) );
			}

			if( js % N == 1 )
			{
				hx_vec = {hx, hx/2. , hx/3.};
				hy_vec = {hy, hy/2. , hy/3.};
				hz_vec = {hz, hz/2. , hz/3.};
			}
			else if(  js % N == 2 )
			{
				hx_vec = {hx/2. , hx/3. , hx/3.};
				hy_vec = {hy/2. , hy/3. , hy/3.};
				hz_vec = {hz/2. , hz/3. , hz/3.};
			}
			else if(jf % (N-1) == 0)
			{
				hx_vec = {hx/3. , hx/3. , hx/2.};
				hy_vec = {hy/3. , hy/3. , hy/2.};
				hz_vec = {hz/3. , hz/3. , hz/2.};
			}
			else if(jf % N == 0)
			{
				hx_vec = {hx/3. , hx/2. , hx};
				hy_vec = {hy/3. , hy/2. , hy};
				hz_vec = {hz/3. , hz/2. , hz};
			}
			else
			{
				hx_vec = {hx/3. , hx/3. , hx/3.};
				hy_vec = {hy/3. , hy/3. , hy/3.};
				hz_vec = {hz/3. , hz/3. , hz/3.};
			}

			if( js%N == 1) // First physical site (it technically correspond to the last site of the bra, but since couplings are homogeneous it does not matter)
			{
				Jxx_vec[0] = Jxx;
				Jyy_vec[0] = Jyy;
				Jzz_vec[0] = Jzz;
			}

			if( jf%N == 0) // Last physical site
			{
				Jxx_vec[1] = Jxx;
				Jyy_vec[1] = Jyy;
				Jzz_vec[1] = Jzz;
			}

			H_S   = hx_vec[0] * Xj[0] * Ij[1] * Ij[2];
			H_S  += hx_vec[1] * Ij[0] * Xj[1] * Ij[2];
			H_S  += hx_vec[2] * Ij[0] * Ij[1] * Xj[2];

			H_S  += hy_vec[0] * Yj[0] * Ij[1] * Ij[2];
			H_S  += hy_vec[1] * Ij[0] * Yj[1] * Ij[2];
			H_S  += hy_vec[2] * Ij[0] * Ij[1] * Yj[2];

			H_S  += hz_vec[0] * Zj[0] * Ij[1] * Ij[2];
			H_S  += hz_vec[1] * Ij[0] * Zj[1] * Ij[2];
			H_S  += hz_vec[2] * Ij[0] * Ij[1] * Zj[2];

			H_SS   = Jxx_vec[0] * Xj[0] * Xj[1] * Ij[2] ;
			H_SS  += Jxx_vec[1] * Ij[0] * Xj[1] * Xj[2] ;

			H_SS  += Jyy_vec[0] * Yj[0] * Yj[1] * Ij[2] ;
			H_SS  += Jyy_vec[1] * Ij[0] * Yj[1] * Yj[2] ;

			H_SS  += Jzz_vec[0] * Zj[0] * Zj[1] * Ij[2] ;
			H_SS  += Jzz_vec[1] * Ij[0] * Zj[1] * Zj[2] ;

			H_SSS  = J_NNN[0] * Xj[0] * Ij[1] * Xj[2] ;
			H_SSS += J_NNN[1] * Yj[0] * Ij[1] * Yj[2] ;
			H_SSS += J_NNN[2] * Zj[0] * Ij[1] * Zj[2] ;
			
	
			vector<int> jn = {j,j+1,j+2};
			if(j<N){
				// action on bra (it is exp(+i H^T t))

				H = -H_S - H_SS - H_SSS;
				H = swapPrime(H,0,1);

				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_bra.push_back(g);
			}
			// action on ket
			else
			{
				H = H_S + H_SS + H_SSS;
				MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
				gates_ket.push_back(g);

			}
		}
		
		
	}

	*/

	// total magnetic field on each site

	for(double q : hx_tot) cerr << q << " ";
	cerr << "\n";
	for(double q : Jzz_tot) cerr << q << " ";
	cerr << "\n";
	// exit(0);

	vector<MyBondGate> gates_ket_reversed = gates_ket;
	vector<MyBondGate> gates_bra_reversed = gates_bra;
	
	reverse(gates_ket_reversed.begin(), gates_ket_reversed.end());
	reverse(gates_bra_reversed.begin(), gates_bra_reversed.end());

	for(MyBondGate g : gates_bra)  		   gates.push_back(g);
	for(MyBondGate g : gates_ket) 		   gates.push_back(g);
	for(MyBondGate g : gates_ket_reversed) gates.push_back(g);
	for(MyBondGate g : gates_bra_reversed) gates.push_back(g);

	return gates;
}





// return either short-range gates or long-range one.
// we do not implement the swap gates as gates, but directly in the TEBD algorithm.

vector<BondGate>
gates_tavis_cummings(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon)
{
    vector<BondGate> gates;

	int N = length(sites);
	
	// -----------------------------------------------------------------------------------
	// gates for H = \omega_0 n_a + \sum_j h \sigma_j^z
	// single-site gate (it implies no need of SVD when applied)
	// in order to use the BondGate class I write as two-gate Hamiltonian
	// When applied, I can exploit the fact there is no need of SVDs

	if(matter_or_photon == "short-range")
	{
	cerr << "Build short-range interactions" << endl;

	for(int j=1 ; j<N;j++)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);

		ITensor hj;
		if(j==1) hj = local_fields[0] * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1] * I_j * S_j1;
		else hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);
		
		ITensor hj;
		if(j==1) hj = local_fields[0] * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1] * I_j * S_j1;
		else hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	return gates;
	}
	// -----------------------------------------------------------------------------------

	// gate H = g \sum_j s_j^+ a + h.c.
	// I implement swap gate as gates using the built-in gate  BondGate(sites,7,8); (see https://itensor.org/docs.cgi?vers=cppv3&page=classes/bondgate)
	// NO: IT DOES NOT WORK WITH DIFFERENT PHYSICAL DIMENSIONS AS IT IS THE CASE OF SPIN-BOSON MODEL
	// Alternative: I implement the swap gate using SVDs -> drawback: I have to explicitly write it in the main
	// photon operators

	else if(matter_or_photon == "long-range")
	{

	Index sph   = sites(1);
	Index sph_p = prime(sites(1));

	ITensor A  = ITensor(sph,sph_p);
	ITensor Ad = ITensor(sph,sph_p);
	
	for(int d=1; d < dim(sph) ; d++)
	{
		A.set( sph(d)  ,sph_p(d+1),sqrt(d));
		Ad.set(sph(d+1),sph_p(d)  ,sqrt(d));
	}
	cerr << "Build photon-matter terms" << endl;

	for(int j=1 ; j<N;j++)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);

			
		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}

		
		ITensor hj = g * (A * Sp + Ad * Sm);
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);
		// gate to swap afterwards 

		// BondGate swap = BondGate(sites,1,j+1);
		// cerr << "swap done" << endl;

		// gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);

			
		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}

		
		ITensor hj = g * (A * Sp + Ad * Sm);
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);
		// gate to swap afterwards 
		// BondGate swap = BondGate(sites,j,j+1);
		// gates.push_back(g);
    }
	return gates;
	}

	else{
		cerr << "string argument inserted is neither 'short-range' nor 'long-range" << endl;
		cerr << "Return empy gates" << endl;
		return gates;
	}
}

// more general light-matter interaction. I could select either
// a) dicke type interactions 
// b) jaynes cummings type interactions


// return either short-range gates or long-range one.
// we do not implement the swap gates as gates, but directly in the TEBD algorithm.

vector<BondGate>
gates_photon_matter(const SiteSet sites , const double omega0 , const double h , const double g, const double dt, string matter_or_photon, string type_of_coupling)
{
    vector<BondGate> gates;

	int N = length(sites);
	
	// -----------------------------------------------------------------------------------
	// gates for H = \omega_0 n_a + \sum_j h \sigma_j^z
	// single-site gate (it implies no need of SVD when applied)
	// in order to use the BondGate class I write as two-gate Hamiltonian
	// When applied, I can exploit the fact there is no need of SVDs

	if(matter_or_photon == "short-range")
	{
	cerr << "Build short-range interactions" << endl;

	for(int j=1 ; j<N;j++)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);

		ITensor hj;
		if(j==1) hj = local_fields[0] * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1] * I_j * S_j1;
		else hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
        
        Index sj = sites(j);
        Index sjp = prime(sites(j));

		Index sj1 = sites(j+1);
        Index sj1p = prime(sites(j+1));
        
        ITensor S_j  = ITensor(sj ,sjp );
		ITensor I_j  = ITensor(sj ,sjp );
		
        ITensor S_j1 = ITensor(sj1,sj1p);
        ITensor I_j1 = ITensor(sj1,sj1p);

		vector<double> local_fields;
	
		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
			local_fields.push_back(omega0);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
			S_j.set(sj(1),sjp(1),1.);
			S_j.set(sj(2),sjp(2),-1.);	
			local_fields.push_back(h);
		}

		if(hasTags(sj1,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj1) ; d++) S_j1.set(sj1(d),sj1p(d),d-1.);
			local_fields.push_back(omega0);

        }
				
		if(hasTags(sj1,"Site,S=1/2"))
		{
			S_j1.set(sj1(1),sj1p(1),1.);
			S_j1.set(sj1(2),sj1p(2),-1.);	
			local_fields.push_back(h);

		}

		for(int d=1; d <= dim(sj)  ; d++) I_j.set(sj(d),sjp(d),1.);
		for(int d=1; d <= dim(sj1) ; d++) I_j1.set(sj1(d),sj1p(d),1.);
		
		ITensor hj;
		if(j==1) hj = local_fields[0] * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		else if(j==N-1) hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1] * I_j * S_j1;
		else hj = local_fields[0]/2. * S_j * I_j1 + local_fields[1]/2. * I_j * S_j1;
		BondGate g = BondGate(sites,j,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);

    }

	return gates;
	}
	// -----------------------------------------------------------------------------------

	// gate H = g \sum_j s_j^+ a + h.c.
	// I implement swap gate as gates using the built-in gate  BondGate(sites,7,8); (see https://itensor.org/docs.cgi?vers=cppv3&page=classes/bondgate)
	// NO: IT DOES NOT WORK WITH DIFFERENT PHYSICAL DIMENSIONS AS IT IS THE CASE OF SPIN-BOSON MODEL
	// Alternative: I implement the swap gate using SVDs -> drawback: I have to explicitly write it in the main
	// photon operators

	else if(matter_or_photon == "long-range")
	{

	Index sph   = sites(1);
	Index sph_p = prime(sites(1));

	ITensor A  = ITensor(sph,sph_p);
	ITensor Ad = ITensor(sph,sph_p);
	
	for(int d=1; d < dim(sph) ; d++)
	{
		A.set( sph(d)  ,sph_p(d+1),sqrt(d));
		Ad.set(sph(d+1),sph_p(d)  ,sqrt(d));
	}
	cerr << "Build photon-matter terms" << endl;

	for(int j=1 ; j<N;j++)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);

			
		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}

		ITensor hj;
		if(type_of_coupling == "dicke")       hj = g * (A + Ad) * (Sp + Sm);
		else if (type_of_coupling == "tavis") hj = g * (A * Sp + Ad * Sm);
		else
		{
			cerr << "Selected type_of_coupling : " << type_of_coupling << " not supported.\n";
			exit(-1); 
		}
		BondGate g = BondGate(sites,1,j+1,BondGate::tReal,dt/2.,hj); 
		gates.push_back(g);
		// gate to swap afterwards 

		// BondGate swap = BondGate(sites,1,j+1);
		// cerr << "swap done" << endl;

		// gates.push_back(g);

    }

	for(int j=N-1 ; j>=1;j--)
    {
		Index sj1   = sites(j+1);
        Index sj1p  = prime(sites(j+1));
        ITensor Sp = ITensor(sj1,sj1p);
        ITensor Sm = ITensor(sj1,sj1p);

			
		if(hasTags(sj1,"Site,S=1/2"))
		{
			Sm.set(sj1(1),sj1p(2),1.);
			Sp.set(sj1(2),sj1p(1),1.);	
		}

		
		ITensor hj;
		if(type_of_coupling == "dicke")       hj = g * (A + Ad) * (Sp + Sm);
		else if (type_of_coupling == "tavis") hj = g * (A * Sp + Ad * Sm);
		else
		{
			cerr << "Selected type_of_coupling : " << type_of_coupling << " not supported.\n";
			exit(-1); 
		}
		BondGate g = BondGate(sites,1,j+1,BondGa	te::tReal,dt/2.,hj); 
		gates.push_back(g);
		// gate to swap afterwards 
		// BondGate swap = BondGate(sites,j,j+1);
		// gates.push_back(g);
    }
	return gates;
	}

	else{
		cerr << "string argument inserted is neither 'short-range' nor 'long-range" << endl;
		cerr << "Return empy gates" << endl;
		return gates;
	}
}


// PXP - returns MPO(s) to perform time evolution of time-step dt.
// 1. We split H = H_1 + H_2 + H_3, so that [H_i,H_j] \neq 0 while the elements within each H_i commute.
// 2. We prepare the gates for H_j, and then we put them inside a time-evolving operator U_j (of time step dt/2) via SVDs. Namely: we construct the gates and then the resulting MPO
// 3. Either we return the vector [U_1,U_2,U_3,U_3,U_2,U_1]. Or we multiply the MPOs in order to have a single one.

vector<MyBondGate>
gates_pxp(const SiteSet sites , const double omega, const double dt)
{

	int N = length(sites);

	
	vector<MyBondGate> gates;

	// first layer (acts on sites [1,2,3] , [4,5,6] , ... )
	for(int j=1 ; j <= N-2 ; j+=3)
	{
		ITensor P1 = (op(sites,"Id",j) + 2*op(sites,"Sz",j))/2;
		ITensor X2 = 2*op(sites,"Sx",j+1);
		ITensor P3 = (op(sites,"Id",j+2) + 2*op(sites,"Sz",j+2))/2;
		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,omega*P1*X2*P3);
		gates.push_back(g);
	}

	// second layer (acts on sites [2,3,4] , [5,6,7] , ... )
	for(int j=2 ; j <= N-2 ; j+=3)
	{
		ITensor P1 = (op(sites,"Id",j) + 2*op(sites,"Sz",j))/2;
		ITensor X2 = 2*op(sites,"Sx",j+1);
		ITensor P3 = (op(sites,"Id",j+2) + 2*op(sites,"Sz",j+2))/2;
		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,omega*P1*X2*P3);
		gates.push_back(g);
	}

	// third layer (acts on sites [3,4,5] , [6,7,8] , ... )
	for(int j=3 ; j <= N-2 ; j+=3)
	{
		ITensor P1 = (op(sites,"Id",j) + 2*op(sites,"Sz",j))/2;
		ITensor X2 = 2*op(sites,"Sx",j+1);
		ITensor P3 = (op(sites,"Id",j+2) + 2*op(sites,"Sz",j+2))/2;
		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,omega*P1*X2*P3);
		gates.push_back(g);
	}


	// if(open_system == true) cerr << "TO DO : Include non-hermitian part" << endl;

	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}



// Rydberg Hamiltonian - we keep up to nearest neighbor interactions

vector<MyBondGate>
gates_rydberg_up_to_VNN(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{

	int N = length(sites);

	vector<MyBondGate> gates;


	for(int j=1 ; j <= N-1 ; j+=1)
	{
		
		vector<ITensor> Nj;
		vector<ITensor> Ij;
		vector<ITensor> Xj;
		for(int q=j ; q<=j+1; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
			Xj.push_back(   op(sites,"Sx",q) );
		}

		// for(ITensor A : Nj) PrintData(A);

		
		double V = Vj[j-1];
		double Omega1 = Omegaj[j-1];
		double Omega2 = Omegaj[j];
		double Delta1 = Deltaj[j-1];
		double Delta2 = Deltaj[j];

		if(j<N-1)
		{
			Omega2 /= 2.;
		 	Delta2 /= 2.;
		}

		if(j>1)
		{
			Omega1 /= 2.;
			Delta1 /= 2.;
		}



		ITensor H_om, H_N, H_NN; 
		H_om = Omega1 * Xj[0] * Ij[1] + Omega2 * Ij[0] * Xj[1];
		H_N  = Delta1 * Nj[0] * Ij[1] + Delta2 * Ij[0] * Nj[1];
		H_NN = V * Nj[0] * Nj[1];

		ITensor H = H_NN + H_N + H_om;
	


		vector<int> jn = {j,j+1};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}
	

	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}

// Rydberg Hamiltonian - we keep up to next-nearest neighbor interactions
// 1. We split H = H_1 + H_2 + H_3, so that [H_i,H_j] \neq 0 while the elements within each H_i commute.
// 2. We prepare the gates for H_j, and then we put them inside a time-evolving operator U_j (of time step dt/2) via SVDs. Namely: we construct the gates and then the resulting MPO
// 3. Either we return the vector [U_1,U_2,U_3,U_3,U_2,U_1]. Or we multiply the MPOs in order to have a single one.

// PLUS: it does not split the single-site terms separately. You earn ~30% in computation time

vector<MyBondGate>
gates_rydberg_up_to_VNNN(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{

	int N = length(sites);

	vector<MyBondGate> gates;
	vector<double> omega;
	vector<double> delta;

	// first layer (acts on sites [1,2,3] , [4,5,6] , ... )
	for(int j=1 ; j <= N-2 ; j+=3)
	{
		int js = j;
		int jf = js + 2;
		
		vector<ITensor> Nj;
		vector<ITensor> Ij;
		vector<ITensor> Xj;

		if(js==1)
			{
			omega = {Omegaj[j-1] , Omegaj[j]/2. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1] , Deltaj[j]/2. , Deltaj[j+1]/3.};
		}
		else if(js==2)
		{
			omega = {Omegaj[j-1]/2. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/2. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}
		else if(jf==N-1)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/2.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/2.};
		}
		else if(jf==N)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/2. , Omegaj[j+1]};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/2. , Deltaj[j+1]};
		}
		else
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}

		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
			Xj.push_back(   op(sites,"Sx",q)) ; 
		}

		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);

		if(j<N-2) V23 /= 2.;		
		if(j>1)   V12 /= 2.;

		ITensor H1, H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;


		H1  =  omega[0] * Xj[0] * Ij[1] * Ij[2];
		H1  += omega[1] * Ij[0] * Xj[1] * Ij[2];
		H1  += omega[2] * Ij[0] * Ij[1] * Xj[2];

		H1  += delta[0] * Nj[0] * Ij[1] * Ij[2];
		H1  += delta[1] * Ij[0] * Nj[1] * Ij[2];
		H1  += delta[2] * Ij[0] * Ij[1] * Nj[2];


		ITensor H = H_NN + H1;

		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}

	// second layer (acts on sites [2,3,4] , [5,6,7] , ... )
	for(int j=2 ; j <= N-2 ; j+=3)
	{
		int js = j;
		int jf = js + 2;
		
		vector<ITensor> Nj;
		vector<ITensor> Ij;
		vector<ITensor> Xj;

		if(js==1)
			{
			omega = {Omegaj[j-1] , Omegaj[j]/2. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1] , Deltaj[j]/2. , Deltaj[j+1]/3.};
		}
		else if(js==2)
		{
			omega = {Omegaj[j-1]/2. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/2. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}
		else if(jf==N-1)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/2.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/2.};
		}
		else if(jf==N)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/2. , Omegaj[j+1]};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/2. , Deltaj[j+1]};
		}
		else
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}

		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
			Xj.push_back(   op(sites,"Sx",q)) ; 

		}

		
		
		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);

		if(j < N-2) V23 /= 2.;
		if(j > 1)   V12 /= 2.;

		ITensor H1, H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;


		H1  =  omega[0] * Xj[0] * Ij[1] * Ij[2];
		H1  += omega[1] * Ij[0] * Xj[1] * Ij[2];
		H1  += omega[2] * Ij[0] * Ij[1] * Xj[2];

		H1  += delta[0] * Nj[0] * Ij[1] * Ij[2];
		H1  += delta[1] * Ij[0] * Nj[1] * Ij[2];
		H1  += delta[2] * Ij[0] * Ij[1] * Nj[2];

		ITensor H = H_NN + H1;

		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}

	// third layer (acts on sites [3,4,5] , [6,7,8] , ... )
	for(int j=3 ; j <= N-2 ; j+=3)
	{
		int js = j;
		int jf = js + 2;
		
		vector<ITensor> Nj;
		vector<ITensor> Ij;
		vector<ITensor> Xj;

		if(js==1)
			{
			omega = {Omegaj[j-1] , Omegaj[j]/2. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1] , Deltaj[j]/2. , Deltaj[j+1]/3.};
		}

		else if(js==2)
		{
			omega = {Omegaj[j-1]/2. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/2. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}
		else if(jf==N-1)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/2.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/2.};
		}
		else if(jf==N)
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/2. , Omegaj[j+1]};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/2. , Deltaj[j+1]};
		}
		else
		{
			omega = {Omegaj[j-1]/3. , Omegaj[j]/3. , Omegaj[j+1]/3.};
			delta = {Deltaj[j-1]/3. , Deltaj[j]/3. , Deltaj[j+1]/3.};
		}


		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
			Xj.push_back(   op(sites,"Sx",q)) ; 
		}
		
		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);

		if(j<N-2) V23 /= 2.;
		if(j>1)   V12 /= 2.;

		ITensor H1, H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;


		H1  =  omega[0] * Xj[0] * Ij[1] * Ij[2];
		H1  += omega[1] * Ij[0] * Xj[1] * Ij[2];
		H1  += omega[2] * Ij[0] * Ij[1] * Xj[2];

		H1  += delta[0] * Nj[0] * Ij[1] * Ij[2];
		H1  += delta[1] * Ij[0] * Nj[1] * Ij[2];
		H1  += delta[2] * Ij[0] * Ij[1] * Nj[2];

		ITensor H = H_NN + H1;

		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}


	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}


// Rydberg Hamiltonian - we keep up to next-nearest neighbor interactions
// 1. We split H = H_1 + H_2 + H_3, so that [H_i,H_j] \neq 0 while the elements within each H_i commute.
// 2. We prepare the gates for H_j, and then we put them inside a time-evolving operator U_j (of time step dt/2) via SVDs. Namely: we construct the gates and then the resulting MPO
// 3. Either we return the vector [U_1,U_2,U_3,U_3,U_2,U_1]. Or we multiply the MPOs in order to have a single one.

// deprecated in favour of gates_rydberg_up_to_VNNN: in the new version we have single-site terms applied together with the 3-site one.

vector<MyBondGate>
gates_rydberg_up_to_VNNN_deprecated(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{

	int N = length(sites);

	vector<MyBondGate> gates;

	// on site terms
	for(int j=1 ; j<= N ; j++)
	{
		ITensor Nj =   (op(sites,"Id",j)   - 2*op(sites,"Sz",j))  /2. ;
		ITensor Xj =    op(sites,"Sx",j);

		ITensor H = Omegaj[j-1] * Xj + Deltaj[j-1] * Nj;
		vector<int> jn = {j};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
		gates.push_back(g);
	}

	// first layer (acts on sites [1,2,3] , [4,5,6] , ... )
	for(int j=1 ; j <= N-2 ; j+=3)
	{
		
		vector<ITensor> Nj;
		vector<ITensor> Ij;

		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
		}

		
		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);
		// cerr << V13 << "\n";

		if(j<N-2) V23 /= 2.;		
		if(j>1)   V12 /= 2.;

		ITensor H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;


		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H_NN);
		gates.push_back(g);
	}
	
	// exit(0);

	// second layer (acts on sites [2,3,4] , [5,6,7] , ... )
	for(int j=2 ; j <= N-2 ; j+=3)
	{
		vector<ITensor> Nj;
		vector<ITensor> Ij;

		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );
		}

		
		
		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);
		// V13 = 0.;

		if(j < N-2) V23 /= 2.;
		if(j > 1)   V12 /= 2.;

		ITensor H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;

		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H_NN);
		gates.push_back(g);
	}

	// third layer (acts on sites [3,4,5] , [6,7,8] , ... )
	for(int j=3 ; j <= N-2 ; j+=3)
	{
		vector<ITensor> Nj;
		vector<ITensor> Ij;

		for(int q=j ; q<=j+2; q++)
		{
			Nj.push_back(  (op(sites,"Id",q)   - 2*op(sites,"Sz",q))  /2. );
			Ij.push_back(   op(sites,"Id",q) );

		}
		
		double V12 = Vj[j-1];
		double V23 = Vj[j];

		double r1 = pow(1/V12, 1./6);
		double r2 = pow(1/V23, 1./6);
		double V13 = pow(1/(r1+r2),6.);
		// V13 = 0.;

		if(j<N-2) V23 /= 2.;
		if(j>1)   V12 /= 2.;

		ITensor H_NN ;
		H_NN  = V12 * Nj[0] * Nj[1] * Ij[2] ;
		H_NN += V23 * Ij[0] * Nj[1] * Nj[2] ;
		H_NN += V13 * Nj[0] * Ij[1] * Nj[2] ;

		vector<int> jn = {j,j+1,j+2};
		MyBondGate g = MyBondGate(sites,jn,dt/2.,H_NN);
		gates.push_back(g);
	}

	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}




vector<MyBondGate>
gates_spin_local_field(const SiteSet sites , vector<double> omegaj, const double dt)
{

	int N = length(sites);
	vector<MyBondGate> gates;

	// H = \sum_j omveja
	for(int j=1 ; j <= N; j++)
	{
		ITensor Sx = 2*op(sites,"Sx",j);
		ITensor Sy = 2*op(sites,"Sy",j);
		ITensor Sz = 2*op(sites,"Sz",j);

		ITensor hj = omegaj[0] * Sx + omegaj[1] * Sy + omegaj[2] * Sz;

		vector<int> jn = {j};

		MyBondGate g = MyBondGate(sites,jn,dt/2.,hj);
		gates.push_back(g);
	}

	
	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}



// using this, we have the more standard usage of indices, but we have to use the ITensor convention
// for S^- , which corresponds to the standard S^+ convention. 
// CHECK: 03.09.23 IT COULD BE THAT I AM MESSING UP INDICES TECHNICALLY, SINCE I AM ASSOCIATING THE INDEX 0 TO THE KET,
// AND PRIME TO THE BRA. BUT FROM ITENSOR DEFAULT CONVENTION IT COULD BE THAT THEY ARE SWAPPED. THIS IS WHY IT LOOKS LIKE
// MY CONVENTION OF S^- IS THE OPPOSITE OF THE ONE OF ITENSOR (IN REALITY THEY ARE NOT DIFFERENT). I SHOULD CHECK THIS
// REWRITING A PIECE OF CODE CONCERNING THIS AND TESTING WITH A NON-HERMITIAN JUMP.
vector<MyBondGateDiss>
gates_local_lindbland(const SiteSet sites , vector<ITensor> Lj, vector<int> lj_sites, vector<double> gammaj , const double dt)
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

		// cerr << inds(Dj) << "\n";
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
gates_local_lindbland(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
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

		// cerr << inds(Dj) << "\n";
		vector<int> jnket = {j};
		vector<int> jnbra = {j};

		MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,Dj);
	
		gates.push_back(g);

	}

	
	// vector<MyBondGate> gates_ = gates;
	// reverse(gates_.begin(), gates_.end());

	// for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}

// keeping as backup - 4.05.23
// vector<MyBondGateDiss>
// gates_local_lindbland(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
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


// 	for(int j=1 ; j <= N; j++)
// 	{
		
// 		Index sj  = sites(j);
// 		Index sj1 = prime(sj);
// 		Index sj2 = prime(sj,2);
// 		Index sj3 = prime(sj,3);

// 		ITensor Idket = ITensor(sj,sj2);
// 		ITensor Idbra = ITensor(sj1,sj3);

// 		for(int q=1 ; q<=dim(sj) ; q++)
// 		{
// 			Idket.set(sj(q),sj2(q),1.);
// 			Idbra.set(sj1(q),sj3(q),1.);
// 		}


// 		ITensor lj = Lj[j-1];
// 		ITensor ljd = conj(lj);

// 		// non-hermitian Hamiltonian

// 		// v2: I think correct version

// 		ljd.mapPrime(0,2);
// 		ITensor ljdlj_I = ljd * lj * Idbra; // acts on ket (sj,sj') -> (sj'',sj''') as desired
	
// 		ljd.mapPrime(2,0);
// 		ljd.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
// 		ITensor I_ljdlj = ljd * lj; // acts on bra  - (3,0)
// 		// I_ljdlj.mapPrime(0,1);      // (3,0) -> (3,1)
// 		I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
// 		I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
// 		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
// 		ljd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
		
// 		// end v2

// 		// start v1: I think it does not do correctly the transposition
// 		/*

// 		ljd.mapPrime(0,2);
// 		ITensor ljdlj_I = ljd * lj * Idbra; // acts on ket (sj,sj') -> (sj'',sj''') as desired
	
// 		ljd.mapPrime(2,0);
// 		ljd.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
// 		ITensor I_ljdlj = ljd * lj; // acts on bra  - (3,0)
// 		I_ljdlj.mapPrime(0,1);      // (3,0) -> (3,1)
// 		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
// 		ljd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
// `		
// 		// end v1
// 		*/ 

// 		// jumps 
// 		lj.mapPrime(1,2);
// 		ljd.mapPrime(1,3);
// 		ljd.mapPrime(0,1);
// 		ITensor lj_ljd = lj * ljd;

	
// 		ITensor Dj = gammaj[j-1] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

// 		// cerr << inds(Dj) << "\n";
// 		vector<int> jnket = {j};
// 		vector<int> jnbra = {j};

// 		MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt,Dj);
	
// 		gates.push_back(g);

// 	}

	
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
gates_nearest_neighbour_local_lindbland(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt)
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
			// ITensor Dj = gamma * ( - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

			// cerr << inds(Dj) << "\n";
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

			// second attempt
			// ljd_.mapPrime(1,3);
			// ljd_.mapPrime(0,1); 
			// lid_.mapPrime(1,3);
			// lid_.mapPrime(0,1); 
			// li_.mapPrime(0,3);
			// lj_.mapPrime(0,3);
			// end second attempt
		
			ITensor I_LdL = ljd_ * li_ + lid_ * lj_; // acts on bra  - (3,0)
			// cerr << I_LdL << endl;
			I_LdL = swapPrime(I_LdL,1,3); // transposition
			// cerr << I_LdL << endl;
			// exit(0);
			// I_LdL.mapPrime(0,1);      // (3,0) -> (3,1)
			// I_LdL.mapPrime(3,1);      // (3,0) -> (1,0)
			// I_LdL.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
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
			// ITensor Dij = gamma * (- 0.5 * LdL_I - 0.5 * I_LdL);
			

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


// vector<MyBondGateDiss>
// gates_nearest_neighbour_local_lindbland(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt)
// {

// 	int N = length(sites);
// 	vector<MyBondGateDiss> gates;
// 	// check size of the two containers

// 	for(MyTrainITensor T : TTrain)
// 	{
		
// 		int i = T.i();
// 		int j = T.j();

// 		ITensor li = T.Ti();
// 		ITensor lj = T.Tj();
// 		ITensor lid = dag(li);
// 		ITensor ljd = dag(lj);

// 		double gamma = T.gamma();

// 		// site index

// 		cerr << "Sites : " << i << " " << j << "\n";

// 		if( abs(i-j)>1 )
// 		{
// 			cerr << "non-local dissipation still not implemented! Returning empty set of gates.\n";
// 			return gates;
// 		}

// 		if(abs(i-j)==0)
// 		{
// 			Index sj  = sites(j);
// 			Index sj1 = prime(sj);
// 			Index sj2 = prime(sj,2);
// 			Index sj3 = prime(sj,3);

// 			ITensor Idket = ITensor(sj,sj2);
// 			ITensor Idbra = ITensor(sj1,sj3);

// 			for(int q=1 ; q<=dim(sj) ; q++)
// 			{
// 				Idket.set(sj(q),sj2(q),1.);
// 				Idbra.set(sj1(q),sj3(q),1.);
// 			}

// 			// non-hermitian Hamiltonian

// 			ljd.mapPrime(0,2);
// 			ITensor ljdlj_I = ljd * lj * Idbra; // acts on ket (sj,sj') -> (sj'',sj''') as desired
		
// 			ljd.mapPrime(2,0);
// 			ljd.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
// 			ITensor I_ljdlj = ljd * lj; // acts on bra  - (3,0)
// 			I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
// 			I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
// 			I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
// 			ljd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
			
// 			// jumps 
// 			lj.mapPrime(1,2);
// 			ljd.mapPrime(1,3);
// 			ljd.mapPrime(0,1);
// 			ITensor lj_ljd = lj * ljd;

// 			ITensor Dj = gamma * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

// 			// cerr << inds(Dj) << "\n";
// 			vector<int> jnket = {j};
// 			vector<int> jnbra = {j};

// 			MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dj);
		
// 			gates.push_back(g);

// 		}

// 		if(abs(i-j)==1)
// 		{
// 			Index si  = sites(i);
// 			Index si1 = prime(si);
// 			Index si2 = prime(si,2);
// 			Index si3 = prime(si,3);

// 			Index sj  = sites(j);
// 			Index sj1 = prime(sj);
// 			Index sj2 = prime(sj,2);
// 			Index sj3 = prime(sj,3);

// 			ITensor Idket = ITensor(si,sj,si2,sj2);
// 			ITensor Idbra = ITensor(si1,sj1,si3,sj3);

// 			ITensor Id_iket_jbra = ITensor(si,si2,sj1,sj3);
// 			ITensor Id_ibra_jket = ITensor(sj,sj2,si1,si3);
						
// 			for(int q=1 ; q<=dim(sj) ; q++)
// 			{
// 				Idket.set( si(q) , si2(q) , sj(q)  , sj2(q) , 1.);
// 				Idbra.set(si1(q) , si3(q) , sj1(q) , sj3(q) , 1.);
// 				Id_iket_jbra.set(si(q),si2(q),sj1(q),sj3(q) , 1.);
// 				Id_ibra_jket.set(sj(q),sj2(q),si1(q),si3(q) , 1.);
// 			}

// 			// effective Hamiltonian part - jumps act either on the ket or bra, but not both

// 			// acts on ket

// 			ITensor lid_ = lid;
// 			ITensor ljd_ = ljd;
// 			ITensor li_  = li;
// 			ITensor lj_  = lj;


// 			lid_.mapPrime(0,2);
// 			lid_.mapPrime(1,0);
// 			ljd_.mapPrime(0,2);
// 			ljd_.mapPrime(1,0);
			
// 			ITensor ljd_li_I = ljd_ * li_ ;
// 			ITensor lid_lj_I = lid_ * lj_ ;
			
// 			ljd_li_I.mapPrime(1,2);
// 			lid_lj_I.mapPrime(1,2);

// 			ITensor LdL_I = (ljd_li_I + lid_lj_I ) * Idbra ; // acts on ket (sj,sj') -> (sj'',sj''') as desired

		
// 			// acts on bra

// 			lid_ = lid;
// 			ljd_ = ljd;
// 			li_  = li;
// 			lj_  = lj;

// 			// first attempt
// 			ljd_.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
// 			lid_.mapPrime(0,3);
// 			li_.mapPrime(1,3);
// 			lj_.mapPrime(1,3);
// 			li_.mapPrime(0,1);
// 			lj_.mapPrime(0,1);

// 			// second attempt
// 			// ljd_.mapPrime(1,3);
// 			// ljd_.mapPrime(0,1); 
// 			// lid_.mapPrime(1,3);
// 			// lid_.mapPrime(0,1); 
// 			// li_.mapPrime(0,3);
// 			// lj_.mapPrime(0,3);
// 			// end second attempt
		
// 			ITensor I_LdL = ljd_ * li_ + lid_ * lj_; // acts on bra  - (3,0)
// 			// I_LdL.mapPrime(0,1);      // (3,0) -> (3,1)
// 			// I_LdL.mapPrime(3,1);      // (3,0) -> (1,0)
// 			// I_LdL.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
// 			I_LdL *= Idket;			// acts on bra from (0,1) to (2,3) as desired
			
// 			// jumps

// 			lid_ = lid;
// 			ljd_ = ljd;
// 			li_  = li;
// 			lj_  = lj;

// 			lj_.mapPrime(1,2);
// 			ljd_.mapPrime(1,3);
// 			ljd_.mapPrime(0,1);

// 			li_.mapPrime(1,2);
// 			lid_.mapPrime(1,3);
// 			lid_.mapPrime(0,1);

// 			ITensor L_Ld = li_ * ljd_ * Id_ibra_jket + lj_ * lid_ * Id_iket_jbra;
 
// 			ITensor Dij = gamma * (L_Ld - 0.5 * LdL_I - 0.5 * I_LdL);

// 			vector<int> jnket = {i,j};
// 			vector<int> jnbra = {i,j};

// 			MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dij);
		
// 			gates.push_back(g);

// 		}

		

// 	}

	
// 	vector<MyBondGateDiss> gates_ = gates;
// 	reverse(gates_.begin(), gates_.end());

// 	for(MyBondGateDiss gate : gates_) gates.push_back(gate);
	
// 	return gates;
// }




// Dissipative gate : we have a list of tensors which act 
// We have a Lidbland of the form L_{i,j} \rho L_{i,j}^\dagger + ...

// NOT USEFUL AT THE MOMENT - NOT TESTED (SHOULD WORK)
vector<MyBondGateDiss>
gates_local_nsites_lindbland(const SiteSet sites , vector<ITensor> Lij_list, vector<vector<int> > Lj_sites, vector<double> gammaj , const double dt)
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

		// vector<ITensor> ljn = Lj_list[k];
		// ITensor lij = ljn[0];
		// ITensor lj = ljn[1];
		// ITensor lid = dag(li);
		// ITensor ljd = dag(lj);

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
		// I_ljdlj.mapPrime(0,1);      // (3,0) -> (3,1)
		I_ljdlj.mapPrime(3,1);      // (3,0) -> (1,0)
		I_ljdlj.mapPrime(0,3);      // (1,0) -> (1,3) I have performed transposition
		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3)
		lijd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
		
		// start v1: I think it does not do correctly the transposition
		/*

		lijd.mapPrime(0,2);
		ITensor ljdlj_I = lijd * lij * Idbra; // acts on ket
	
		lijd.mapPrime(2,0);
		lijd.mapPrime(0,3); // (0,1) -> (3,1) (prime order)
		ITensor I_ljdlj = lijd * lij; // acts on bra  - (3,0)
		I_ljdlj.mapPrime(0,1);      // (3,0) -> (3,1)
		I_ljdlj *= Idket;			// acts on bra from (0,1) to (2,3) as desired
		lijd.mapPrime(3,0); // (0,1) -> (3,1) (prime order)
`		
		// end v1
		*/ 

		// jumps 
		lij.mapPrime(1,2);
		lijd.mapPrime(1,3);
		lijd.mapPrime(0,1);
		ITensor lj_ljd = lij * lijd;

	
		ITensor Dj = gammaj[k] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

		// cerr << inds(Dj) << "\n";
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


// vector<MyBondGate>
// MPO_global_dissipation_plane_waves(const double gamma, const double k0, const douple sign)
// {
// 	complex<double> ci;
// 	complex<double> i(0,sign*1);
// 	auto ampo = AutoMPO(sites);
// 	for(int j = 1; j < N; ++j)
// 		{
// 		ampo += exp(i * k0 * (j-1)),"S-",j;
// 		}

// 	//Convert the AutoMPO object to an MPO
// 	MPO L = toMPO(ampo);

// }




// vector<MyBondGate>
// MPO_global_dissipation_plane_waves(const double gamma, const double k0, const douple sign)
// {
// 	complex<double> ci;
// 	complex<double> i(0,sign*1);
// 	auto ampo = AutoMPO(sites);
// 	for(int j = 1; j < N; ++j)
// 		{
// 		ampo += exp(i * k0 * (j-1)),"S-",j;
// 		}

// 	//Convert the AutoMPO object to an MPO
// 	MPO L = toMPO(ampo);

// }

MPS
TEBD_lindbland_time_evolve(MPS psi_t, vector<BondGate> gates , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");
    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

    	gateTEvol( gates , dt , dt , psi_t , TEBD_args); 

        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                

                int j = jket[0];

                cerr <<  j << " "; 

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }

            gateTEvol( gates , dt , dt , psi_t , TEBD_args); 
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}



MPS
TDVP_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;

	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 10;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",false,
                                      "Silent",true});
            }
        
        // TDVP sweep

        double energy = tdvp(psi_t,H,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

    
        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                

                int j = jket[0];

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",false,
                                      "Silent",true});
            }
            // TDVP sweep
			energy = tdvp(psi_t,H,-1_i*dt_step,sweeps,{"Truncate",true,
										"DoNormalize",false,
										"Silent",true,
										"NumCenter",2,
										"ErrGoal",1E-10});
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}

MPS
TDVP_split_lindbland_time_evolve(MPS psi_t, MPO Hbra , MPO Hket, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;

	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 10;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        // if(k < 3)
        //     {
        //     // Global subspace expansion
        //     std::vector<Real> epsilonK = {1E-12, 1E-12};
        //     addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
        //                               "Method","DensityMatrix",
        //                               "KrylovOrd",3,
        //                               "DoNormalize",false,
        //                               "Silent",true});
        //     }
        
        // TDVP sweep

        double energy = tdvp(psi_t,Hbra,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

	
		energy = tdvp(psi_t,Hket,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

    
        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                

                int j = jket[0];

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			// if(k < 3)
            // {
            // // Global subspace expansion
            // std::vector<Real> epsilonK = {1E-12, 1E-12};
            // addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
            //                           "Method","DensityMatrix",
            //                           "KrylovOrd",3,
            //                           "DoNormalize",false,
            //                           "Silent",true});
            // }
            // TDVP sweep
			energy = tdvp(psi_t,Hbra,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});

	
			energy = tdvp(psi_t,Hket,-1_i*dt_step,sweeps,{"Truncate",true,
                                        "DoNormalize",false,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}



MPS
MPO_lindbland_time_evolve(MPS psi_t, MPO H , vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");

	double dt_step = dt;
	if(dissipative) dt_step = dt/2.;
	
	MPO Ht =  -1_i * dt_step * H; 

   
    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;

        MPS dpsi = applyMPO(Ht,psi_t,{"Method=","DensityMatrix","MaxDim=",MaxDim,"Cutoff=",cut_off});
        psi_t = sum(psi_t, dpsi.noPrime());
	
        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                

                int j = jket[0];

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			dpsi = applyMPO(Ht,psi_t,{"Method=","DensityMatrix","MaxDim=",MaxDim,"Cutoff=",cut_off});
        	psi_t = sum(psi_t, dpsi.noPrime());
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
}


// Time evolution using swap gates to handle long range interactions

MPS
TEBD_long_range_int_lindbland_time_evolve(MPS psi_t, vector<BondGate> gates_H, vector<MyBondGateDiss> gates_D , Args TEBD_args, bool dissipative , double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int total_steps = int(T/dt);
	int MaxDim     = TEBD_args.getInt("MaxDim");
	double cut_off = TEBD_args.getReal("Cutoff");
	// cerr << "Total steps : " << total_steps << "\n";
	int N2 = length(psi_t);
	int N = int(N2/2);

	// cerr << "Start evolution.\n";

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = t_start + (k+1)*dt;
		cerr << "Time : " << t << "\n";
		// long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
        for (BondGate g : gates_H)
        {   
            int j1 = g.i1();
			int j2 = g.i2();
			int j ;

			// cerr << "Applying swap gate between " << j1 << " " << j2 << "\n";

			// the swap gates move around the impurity site, to make it near 
			// at then end, it will be back to the original position
			// it acts on bra
			if(j1 <= N && j2 <= N)
			{
				j = min(j1,j2);
				ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
				auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
				
				psi_t.set(j,U);
				psi_t.set(j+1,S*V);
				swap_gate(&psi_t,j,j+1,cut_off,MaxDim);

			}

			// it acts on ket

			else
			{
				j = max(j1,j2);
				ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
				auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
				psi_t.set(j-1,U);
				psi_t.set(j,S*V);
				swap_gate(&psi_t,j-1,j,cut_off,MaxDim);
			}

            // // change orthogonality center to minimize error
            // psi.position(j);
            // psi.normalize();
            // apply gate
            
        }

    
        if(dissipative)
        {
            for (MyBondGateDiss gate : gates_D)
            {
                vector<int> jket = gate.jnket(); // sites where it acts on ket
                ITensor g        = gate.gate();
                

                int j = jket[0];

                ITensor AA = psi_t(j) * psi_t(j+1);
                ITensor dpsi =  g * AA;
                dpsi.mapPrime(1,0);

                AA = AA + dpsi;

                auto [U,S,V] = svd(AA,inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
                psi_t.set(j,U);
                psi_t.set(j+1,S*V);

            }


			// long range interaction gates -> need to swap gates (see https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043255)
			for (BondGate g : gates_H)
			{   

				int j1 = g.i1();
				int j2 = g.i2();
				int j ;
				// cerr << "Applying swap gate between " << j1 << " " << j2 << "\n";

				// the swap gates move around the impurity site, to make it near 
				// at then end, it will be back to the original position
				// it acts on bra
				if(j1 <= N && j2 <= N)
				{
					j = min(j1,j2);
					ITensor AA = psi_t(j)*psi_t(j+1)*g.gate();
					auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
					
					psi_t.set(j,U);
					psi_t.set(j+1,S*V);
					swap_gate(&psi_t,j,j+1,cut_off,MaxDim);

				}

				// it acts on ket

				else
				{
					j = max(j1,j2);
					ITensor AA = psi_t(j-1)*psi_t(j)*g.gate();
					auto [U,S,V] = svd(noPrime(AA),inds(psi_t(j-1)),{"Cutoff=",cut_off,"MaxDim=",MaxDim});
					psi_t.set(j-1,U);
					psi_t.set(j,S*V);
					swap_gate(&psi_t,j-1,j,cut_off,MaxDim);
				}
			}
            
        }


        if(normalize)
        {
            double norm = compute_norm_purifed_impurity(&psi_t);
            psi_t /= norm;
        }
        


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";
			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
				return psi_t;
			}
            writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       
    

	return psi_t;
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

vector<BondGate>
gates_coherent_unfolded_kondo_impurity_model_energy_basis(const SiteSet sites , const vector<double> J, const vector<double> hup, const vector<double> hdn, const double dt)
{
	// N is always even and corresponds to the physical size
	int N2 = length(sites);
	int N = N2/2;


    vector<BondGate> gates;
	vector<BondGate> gates_ket;

	// collection of gates

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

	
	for(int q=1 ; q<=N2; q++)
	{
		Adagup.push_back(op(sites,"Adagup",q) );
		Aup.push_back(   op(sites,"Aup",q) );
		Adagdn.push_back(op(sites,"Adagdn",q) );
		Adn.push_back(   op(sites,"Adn",q) );
		Nup.push_back(   op(sites,"Nup",q) );
		Ndn.push_back(   op(sites,"Ndn",q) );
		Id.push_back(op(sites,"Id",q));
		AdagupFi.push_back(op(sites,"Adagup*F",q));
		AupFi.push_back(op(sites,"Aup*F",q));
		FiAdn.push_back(op(sites,"F*Adn",q));
		FiAdagdn.push_back(op(sites,"F*Adagdn",q));

	}
	
	cerr << "Couplings.\n\n";

	for(double z : J) cerr << z << " ";
	cerr << "\n\n";

	// Action of the Hamiltonian on the bra
	// we start from the center of the chain

	cerr << "Entering bra cycle.\n\n";
	for(int j=N-1 ; j>=1 ; j--)
	{

		ITensor H , H_S , H_SS;

		// build single site - the impurity is set only at the first bond
		if(j==N-1)
		{
			H_S = (hup[N-j-1] * Nup[N-1] + hdn[N-j-1] * Ndn[N-1]) * Id[N-2] ;
			H_S += Id[N-1] * (hup[N-j] * Nup[j-1] + hdn[N-j] * Ndn[j-1]);
		}
		else  H_S = Id[N-1] * (hup[N-j] * Nup[j-1] + hdn[N-j] * Ndn[j-1]);
		
		

		// cerr << "Site : " << j     << " -> h : " << hup[j-1] << "\n";
		// cerr << "Site : " << j + 1 << " -> h : " << hup[j] << "\n";

		// two sites hopping
		
		H_SS   = J[N-j-1] * (  Aup[j-1] * AdagupFi[N-1] - Adagup[j-1]   * AupFi[N-1] );
		H_SS  += J[N-j-1] * (  FiAdn[j-1] * Adagdn[N-1] - FiAdagdn[j-1] * Adn[N-1]   );
		
		cerr << J[N-j-1] << " ";

		
		H = - H_S - H_SS ;
		H = swapPrime(H,0,1); 
		BondGate g = BondGate(sites,j,N,BondGate::tReal,dt/2.,H); 
		gates.push_back(g);
			
	
    }
	cerr << "\n\n";

	vector<BondGate> gates_reversed = gates;
	reverse(gates_reversed.begin(), gates_reversed.end());
	for(BondGate g : gates_reversed)  gates.push_back(g);


	cerr << "Entering ket cycle.\n\n";

	for(int j=N+2 ; j<=N2 ; j++)
	{

		ITensor H , H_S , H_SS;

		// build single site - the impurity is set only at the first bond
		if(j==N+2)
		{
			H_S =  (hup[abs(N-j)-1] * Nup[N] + hdn[abs(N-j)-1] * Ndn[N]) * Id[N+1] ;
			H_S += Id[N] * (hup[abs(N-j)] * Nup[j-1] + hdn[abs(N-j)] * Ndn[j-1]);
		}
		else H_S = Id[N] * (hup[abs(N-j)] * Nup[j-1] + hdn[abs(N-j)] * Ndn[j-1]);
		

		// cerr << "Site : " << j     << " -> h : " << hup[j-1] << "\n";
		// cerr << "Site : " << j + 1 << " -> h : " << hup[j] << "\n";

		// two sites hopping
		
		H_SS   = J[abs(N-j)-2] * (  AdagupFi[N] * Aup[j-1]   - AupFi[N] * Adagup[j-1] );
		H_SS  += J[abs(N-j)-2] * (  Adagdn[N]   * FiAdn[j-1] - Adn[N] * FiAdagdn[j-1] );
		
		cerr << J[abs(N-j)-2] << " ";

		H = H_S + H_SS ;
		BondGate g = BondGate(sites,N+1,j,BondGate::tReal,dt/2.,H); 
		gates_ket.push_back(g);
			
	
    }

	gates_reversed = gates_ket;
	reverse(gates_reversed.begin(), gates_reversed.end());
	for(BondGate g : gates_reversed)  gates_ket.push_back(g);
	for(BondGate g : gates_ket) gates.push_back(g);

	// cerr << "\n\n";

	// cerr << "Returning gates\n";
	
	// for(BondGate g : gates) cerr << g.i1() << " ";
	// cerr << "\n\n";

	// for(BondGate g : gates) cerr << g.i2() << " ";

	// exit(0);

	return gates;

}



MPS
TDVP_time_evolve(MPS psi_t, SiteSet sites, MPO H , Args TDVP_args, double dt , double T , int steps_save_state, bool normalize, string file_root, double t_start)
{
	int N = length(sites);
	
 
    
    string file_obs    = tinyformat::format("%s.txt",file_root);
    ofstream save_file( file_obs) ;
    save_file << "# t . norm \n";

    file_obs    = tinyformat::format("%s_nj_up.txt",file_root);
    ofstream save_file_nj_up( file_obs) ;

    file_obs    = tinyformat::format("%s_nj_dn.txt",file_root);
    ofstream save_file_nj_dn( file_obs) ;

    int total_steps = int(T/dt);
	int MaxDim = TDVP_args.getInt("MaxDim");
	double cut_off = TDVP_args.getReal("Cutoff");
	auto sweeps = Sweeps(1);
    sweeps.maxdim() = MaxDim;
    sweeps.cutoff() = cut_off;
    sweeps.niter() = 20;

    for(int k=0 ; k<total_steps ; k++)
    {
        double t = (k+1)*dt;

        if(k < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi_t,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",normalize,
                                      "Silent",true});
            }
        
        // TDVP sweep

        double energy = tdvp(psi_t,H,-1_i*dt,sweeps,{"Truncate",true,
                                        "DoNormalize",normalize,
                                        "Silent",true,
                                        "NumCenter",2,
                                        "ErrGoal",1E-10});


        double normalization = norm(psi_t);
        
        save_file << t << " " << normalization << "\n";
        save_file_nj_up << t ;
		save_file_nj_dn << t ;
		for(int j : range1(N))
		{
			psi_t.position(j);
			ITensor ket = psi_t(j);
			ITensor bra = dag(prime(ket,"Site"));
			ITensor Njup = op(sites,"Nup",j);
			ITensor Njdn = op(sites,"Ndn",j);

			save_file_nj_up << " " << real(eltC(ket * Njup * bra))/normalization ;
        	save_file_nj_dn << " " << real(eltC(ket * Njdn * bra))/normalization ; 
		}
		
		save_file_nj_up << "\n";
		save_file_nj_dn << "\n";
        
        save_file_nj_dn.flush();
        save_file_nj_up.flush();


        if ( (k+1) % steps_save_state == 0)
        {
			cerr << "Saving state time : " << t << "\n";

			if( maxLinkDim(psi_t) > MaxDim)
			{
				cerr << "Reached max bond dimension. Abort.\n";
                exit(0);
			}
            // writeToFile(tinyformat::format("%s_psi_t%.3f",file_root,t),psi_t); 
        }
    }       


	return psi_t;
}