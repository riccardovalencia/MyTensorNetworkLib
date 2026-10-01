/**
 * @file rydberg.cc
 * @brief Implementation of rydberg.h (the functions are documented in the header).
 */
#include "rydberg.h"
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


MPO
mpo_pxp(const SiteSet s, const double omega)
{
	int N = length(s);
    MPO H = MPO(s);

    // list of necessary operators (nb I use convention S are spin-matrices)
    // row and column index
    Index i_idx = Index(2);
    Index j_idx = Index(2);    

    // projectors (1+2*Sz)/2
    ITensor P = ITensor(i_idx,j_idx);
    P.set(i_idx(1),j_idx(1),1.);

    // 2*Sx  
    ITensor X = ITensor(i_idx,j_idx);
    X.set(i_idx(1),j_idx(2),1.);
    X.set(i_idx(2),j_idx(1),1.);

    // Identity
    ITensor I = ITensor(i_idx,j_idx);
    I.set(i_idx(1),j_idx(1),1.);
    I.set(i_idx(2),j_idx(2),1.);


    // build MPO
    // link indeces
    Index r_idx, l_idx;

    for(int j=1; j <= N; j++)
    {
        Index sj = s(j);
        Index sjp = prime(s(j));
        ITensor Hj;

        if(j==1)
        {
            r_idx = Index(2,tinyformat::format("l=%d,Link",j));
            Hj = ITensor(sj,sjp,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set I to the first bond-index
                    Hj.set(sj(l),sjp(k),r_idx(1),omega*elt(I,i_idx(l),j_idx(k)));
                    // set P to the second bond-index
                    Hj.set(sj(l),sjp(k),r_idx(2),omega*elt(P,i_idx(l),j_idx(k)));
                }
            }

        }

        if(j==2)
        {
            l_idx = r_idx;
            r_idx = Index(3,tinyformat::format("l=%d,Link",j));
            Hj = ITensor(sj,sjp,l_idx,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (1,1)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(1),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (1,2)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(2),elt(P,i_idx(l),j_idx(k)));

                    // set X to element (2,3)
                    Hj.set(sj(l),sjp(k),l_idx(2),r_idx(3),elt(X,i_idx(l),j_idx(k)));
                }
            }
        }

        if(j==3)
        {
            l_idx = r_idx;
            r_idx = Index(4,tinyformat::format("l=%d,Link",j));

            Hj = ITensor(sj,sjp,l_idx,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (1,1)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(1),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (1,2) and (3,4)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(2),elt(P,i_idx(l),j_idx(k)));
                    Hj.set(sj(l),sjp(k),l_idx(3),r_idx(4),elt(P,i_idx(l),j_idx(k)));

                    // set X to element (2,3)
                    Hj.set(sj(l),sjp(k),l_idx(2),r_idx(3),elt(X,i_idx(l),j_idx(k)));
                }
            }
        }

        if(j >=4 && j <= N-3)
        {
            l_idx = r_idx;
            r_idx = Index(4,tinyformat::format("l=%d,Link",j));

            Hj = ITensor(sj,sjp,l_idx,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (1,1) and (4,4)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(1),elt(I,i_idx(l),j_idx(k)));
                    Hj.set(sj(l),sjp(k),l_idx(4),r_idx(4),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (1,2) and (3,4)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(2),elt(P,i_idx(l),j_idx(k)));
                    Hj.set(sj(l),sjp(k),l_idx(3),r_idx(4),elt(P,i_idx(l),j_idx(k)));

                    // set X to element (2,3)
                    Hj.set(sj(l),sjp(k),l_idx(2),r_idx(3),elt(X,i_idx(l),j_idx(k)));
                }
            }
        }

        if(j==N-2)
        {
            l_idx = r_idx;
            r_idx = Index(4,tinyformat::format("l=%d,Link",j));

            Hj = ITensor(sj,sjp,l_idx,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (4,4)
                    Hj.set(sj(l),sjp(k),l_idx(4),r_idx(4),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (1,2) and (3,4)
                    Hj.set(sj(l),sjp(k),l_idx(1),r_idx(2),elt(P,i_idx(l),j_idx(k)));
                    Hj.set(sj(l),sjp(k),l_idx(3),r_idx(4),elt(P,i_idx(l),j_idx(k)));

                    // set X to element (2,3)
                    Hj.set(sj(l),sjp(k),l_idx(2),r_idx(3),elt(X,i_idx(l),j_idx(k)));
                }
            }
        }

        if(j==N-1)
        {
            l_idx = r_idx;
            r_idx = Index(4,tinyformat::format("l=%d,Link",j));

            Hj = ITensor(sj,sjp,l_idx,r_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (4,4)
                    Hj.set(sj(l),sjp(k),l_idx(4),r_idx(4),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (3,4)
                    Hj.set(sj(l),sjp(k),l_idx(3),r_idx(4),elt(P,i_idx(l),j_idx(k)));

                    // set X to element (2,3)
                    Hj.set(sj(l),sjp(k),l_idx(2),r_idx(3),elt(X,i_idx(l),j_idx(k)));
                }
            }
        }

        if(j==N)
        {
            l_idx = r_idx;

            Hj = ITensor(sj,sjp,l_idx);

            for(int l=1; l<= dim(sj); l++)
            {
                for(int k=1; k<= dim(sjp); k++)
                {
                    // set Id to element (4,)
                    Hj.set(sj(l),sjp(k),l_idx(4),elt(I,i_idx(l),j_idx(k)));

                    // set P in elemennt (3,)
                    Hj.set(sj(l),sjp(k),l_idx(3),elt(P,i_idx(l),j_idx(k)));
                }
            }
        }


        H.ref(j) = Hj;
        
    }


    return H;

}


// Gates of the PXP Hamiltonian
// H = omega \sum_j P_j X_{j-1} P_{j+1}
// where P_j = (1+Z_j)/2

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


// given spatial configurations, it returs the potentials 1/rj^alpha
 
vector<double>
compute_potential(const vector< vector<double> > rj , const double alpha)
{

    vector<double> V;
    int N = rj.size();

    for(int j = 0 ; j < N - 1 ; j++)
    {
        double drx = rj[j][0] - rj[j+1][0];
        double dry = rj[j][1] - rj[j+1][1];
        double drz = rj[j][2] - rj[j+1][2];

        double d = sqrt(drx*drx + dry*dry + drz*drz );

        V.push_back(pow(1/d,alpha) );        

    }

    return V ; 
}
