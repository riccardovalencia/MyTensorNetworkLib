#include "TEBD_testing.h"
#include <itensor/all.h>
#include <iostream>
#include <math.h>       
#include <complex>
using namespace std;
using namespace itensor;

// Rydberg Hamiltonian - we keep up to nearest neighbor interactions

vector<MyBondGate>
gates_rydberg_up_to_VNN_testing(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
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

vector<MyBondGate>
gates_rydberg_up_to_VNNN_testing(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{

	int N = length(sites);

	vector<MyBondGate> gates;

	// on site terms
	// for(int j=1 ; j<= N ; j++)
	// {
	// 	ITensor Nj =   (op(sites,"Id",j)   - 2*op(sites,"Sz",j))  /2. ;
	// 	ITensor Xj =    op(sites,"Sx",j);

	// 	ITensor H = Omegaj[j-1] * Xj + Deltaj[j-1] * Nj;
	// 	vector<int> jn = {j};
	// 	MyBondGate g = MyBondGate(sites,jn,dt/2.,H);
	// 	gates.push_back(g);
	// }

	vector<double> omega;
	vector<double> delta;

	vector<double> Omega_check;
	for(int j=1; j<=N ;j++) Omega_check.push_back(0.);

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

		Omega_check[j-1] += delta[0];
		Omega_check[j]   += delta[1];
		Omega_check[j+1] += delta[2];
		

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
	// exit(0);

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

		Omega_check[j-1] += delta[0];
		Omega_check[j]   += delta[1];
		Omega_check[j+1] += delta[2];

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
		// V13 = 0.;

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

		Omega_check[j-1] += delta[0];
		Omega_check[j]   += delta[1];
		Omega_check[j+1] += delta[2];
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
		// V13 = 0.;

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

	// for(double om : Omega_check) cerr << om << "\n";
	// exit(0);

	vector<MyBondGate> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGate gate : gates_) gates.push_back(gate);
	
	return gates;
}




// For lidbland we use a convention of indices so that the lowering operator S^- is given by [[0,0],[1,0]]
// instead of [[0,1],[0,0]] given by ITensor. Despite being in some sense more natural, it leads to messy 
// handling of the indices. For this reason, it is DEPRECATED to use it. Based on my tests, it works correctly.

vector<MyBondGateDiss>
gates_local_lindbland_testing(const SiteSet sites , vector<ITensor> Lj, vector<double> gammaj , const double dt)
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
		ljd_.mapPrime(1,2);  
		ITensor ljdlj_I = ljd_ * lj_ ;
		ljdlj_I.mapPrime(1,0);
		ljdlj_I *= Idbra;   

		// bra
		ITensor I_ljdlj = ljd_ * lj_; // it has index (2,1)
		// perform transpose of I_ljdlj
		I_ljdlj.mapPrime(1,3);
		I_ljdlj.mapPrime(2,1);
		I_ljdlj *= Idket;	

		// reset operators

		ljd_ = ljd;
		lj_  = lj;
		
		// jumps
		// before 
		lj_.mapPrime(0,2);
		lj_.mapPrime(1,0);
		// now - test
		// lj_.mapPrime(1,2);
		// 

		ljd_.mapPrime(0,3);
		// // test
		// ljd_.mapPrime(1,3);
		// ljd_.mapPrime(0,1);
		
		

		ITensor lj_ljd = lj_ * ljd_;
		

		// ITensor Dj = gammaj[j-1] * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);
		ITensor Dj = gammaj[j-1] * lj_ljd ;

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

// Non-diagonal local Lindland 
// We consider Lindbland of the form: Li \rho L_{i+1}^\dagger + 0.5 * {Li L_{i+1}, \rho}
// It appears as a 4-sites gate if we apply the jumps and the non-hermitian part at the same time.
// A possibility is to split differently: the non hermitiain part in the hermitian one
// In this way I have at most 2-sites gates. The number of gates that have to be applied is the same.  <- POSSIBLE EFFICIENCY GAIN(?)
// Drawback: we would have long-range interactions in the final case study both in the Hamiltonian part 
// and jump part -> MULTIPLE LOOPS NECESSARY                                                           -> HUGE INEFFICENCY FROM LOOPING

// MyTrainITensor is a personalized class containing ITensors which have to act either on 


// It is a 4-sites object
// We apply Li \rho L_j^\dagger + L_j \rho L_i^\dagger - 1/2( {L_i^\dagger L_j , \rho} + {L_j^\dagger L_i,\rho} ) (OR SIMILAR)
vector<MyBondGateDiss>
gates_nearest_neighbour_local_lindbland_testing(const SiteSet sites , vector<MyTrainITensor> TTrain, const double dt)
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

			ITensor lj_ = lj;
			ITensor ljd_ = ljd;

			// non-hermitian Hamiltonian
			// ket 
			ljd_.mapPrime(1,2);  
			ITensor ljdlj_I = ljd_ * lj_ ;
			ljdlj_I.mapPrime(1,0);
			ljdlj_I *= Idbra;   

			// bra
			ITensor I_ljdlj = ljd_ * lj_; // it has index (2,1)
			// perform transpose of I_ljdlj
			I_ljdlj.mapPrime(1,3);
			I_ljdlj.mapPrime(2,1);
			I_ljdlj *= Idket;	

			// reset operators

			ljd_ = ljd;
			lj_  = lj;
			
			// jumps 
			lj_.mapPrime(0,2);
			lj_.mapPrime(1,0);
			ljd_.mapPrime(0,3);

			ITensor lj_ljd = lj_ * ljd_;

			ITensor Dj = gamma * (lj_ljd - 0.5 * ljdlj_I - 0.5 * I_ljdlj);

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

			// useful ITensors

			ITensor lid_ = lid;
			ITensor ljd_ = ljd;
			ITensor li_  = li;
			ITensor lj_  = lj;

			// effective Hamiltonian part - jumps act either on the ket or bra, but not both

			// non-hermitian Hamiltonian
			  
// // old
// ITensor lid_ = lid;
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


// // old

			// acts on ket

			lid_.mapPrime(1,2);
			ljd_.mapPrime(1,2); 
			// li_.mapPrime(1,2);
			// lj_.mapPrime(1,2);
			li_.mapPrime(0,2);
			lj_.mapPrime(0,2);
			li_.mapPrime(1,0);
			lj_.mapPrime(1,0);
			
			
			ITensor LdL_I = ljd_ * li_ + lid_ * lj_;
			// LdL_I.mapPrime(1,0);
			LdL_I *= Idbra ; 

		
			// acts on bra
			// lid_ = lid;
			// ljd_ = ljd;
			// li_  = li;
			// lj_  = lj;
			// lid_.mapPrime(1,3);
			// lid_.mapPrime(0,1);

			// ljd_.mapPrime(1,3); 
			// ljd_.mapPrime(0,1); 

			// li_.mapPrime(1,3);
			// li_.mapPrime(0,1);

			// lj_.mapPrime(1,3);
			// lj_.mapPrime(0,1);

			// bra
			ITensor I_LdL = ljd_ * li_ + lid_ * lj_; // it has index (2,1)
			// perform transpose of I_ljdlj
			I_LdL.mapPrime(0,3);
			I_LdL.mapPrime(2,1);
			I_LdL *= Idket;	

			// jumps 
			
			lid_ = lid;
			ljd_ = ljd;
			li_  = li;
			lj_  = lj;

			lj_.mapPrime(0,2);
			lj_.mapPrime(1,0);
			ljd_.mapPrime(0,3);

			li_.mapPrime(0,2);
			li_.mapPrime(1,0);
			lid_.mapPrime(0,3);


			// PrintData(li_);
			// PrintData(lj_);

			ITensor L_Ld = li_ * ljd_ * Id_ibra_jket + lj_ * lid_ * Id_iket_jbra;
 
			// cerr << inds(L_Ld) << endl;
			// cerr << "\n\n";
			// cerr << inds(LdL_I) << endl;
			// cerr << "\n\n";

			// cerr << inds(I_LdL) << endl;
			// cerr << equals(inds(L_Ld),inds(LdL_I)) << endl;
			// cerr << equals(inds(LdL_I),inds(I_LdL)) << endl;

			ITensor Dij = gamma * (L_Ld - 0.5 * LdL_I - 0.5 * I_LdL);
			// ITensor Dij = gamma * I_LdL ;

			vector<int> jnket = {i,j};
			vector<int> jnbra = {i,j};

			// PrintData(Dij);
			// cerr << gamma << endl;
			// cerr << dt << endl;
			// exit(0);

			MyBondGateDiss g = MyBondGateDiss(sites,jnket,jnbra,dt/2.,Dij);
		
			gates.push_back(g);

		}

		

	}

	
	vector<MyBondGateDiss> gates_ = gates;
	reverse(gates_.begin(), gates_.end());

	for(MyBondGateDiss gate : gates_) gates.push_back(gate);
	
	return gates;
}
