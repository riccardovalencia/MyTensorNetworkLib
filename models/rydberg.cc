/**
 * @file rydberg.cc
 * @brief Implementation of rydberg.h (interfaces documented in the header, logic commented here).
 */
#include "rydberg.h"
#include "../mps/gates.h"
#include <itensor/all.h>
#include <functional>
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


// MPO written by hand (bond dimension 4): the link value records how much of a term P X P has been
// placed, 1 = nothing (identity), 2 = P, 3 = P X, 4 = complete (identity afterwards). The first and
// last sites and their neighbours use the allowed subsets of these transitions; omega sits on site 1.
MPO
make_pxp_mpo(const SiteSet s, const double omega)
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


// three layers of non-overlapping three-site gates: [1,2,3] [4,5,6] ..., [2,3,4] ..., [3,4,5] ...,
// each with the term make_term(j) on (j, j+1, j+2), followed by the reversed sequence

static vector<TebdGate>
make_three_site_layers(const int N, const double dt, const function<ITensor(int)>& make_term)
{
    vector<TebdGate> gates;
    for(int layer = 1 ; layer <= 3 ; layer++)
        for(int j = layer ; j <= N-2 ; j += 3)
            gates.push_back(TebdGate({j,j+1,j+2}, dt/2., make_term(j)));
    return make_symmetric_sweep(gates);
}


// Gates of the PXP Hamiltonian H = omega sum_j P_j X_{j+1} P_{j+2}, P_j = (1+Z_j)/2: one
// three-site term per triple, in the three layers of make_three_site_layers

vector<TebdGate>
make_pxp_gates(const SiteSet sites , const double omega, const double dt)
{
    function<ITensor(int)> pxp = [&](int j)
    {
        ITensor P1 = (op(sites,"Id",j) + 2*op(sites,"Sz",j))/2;
        ITensor X2 = 2*op(sites,"Sx",j+1);
        ITensor P3 = (op(sites,"Id",j+2) + 2*op(sites,"Sz",j+2))/2;
        return omega*P1*X2*P3;
    };
    return make_three_site_layers(length(sites), dt, pxp);
}


// Rydberg Hamiltonian with nearest-neighbour interactions: one two-site gate per bond with
// V n_j n_{j+1} and the on-site terms (Omega S^x, Delta n) divided among the bonds sharing them

vector<TebdGate>
make_rydberg_gates_nn(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{
    int N = length(sites);
    vector<TebdGate> gates;
    for(int j = 1 ; j <= N-1 ; j++)
    {
        vector<ITensor> Nj, Ij, Xj;
        for(int q = j ; q <= j+1 ; q++)
        {
            Nj.push_back( (op(sites,"Id",q) - 2*op(sites,"Sz",q)) / 2. );
            Ij.push_back(  op(sites,"Id",q) );
            Xj.push_back(  op(sites,"Sx",q) );
        }

        // on-site terms shared with the neighbouring bonds
        double count1 = count_gates_containing(j,   1, 2, N);
        double count2 = count_gates_containing(j+1, 1, 2, N);

        ITensor H_om = Omegaj[j-1] / count1 * Xj[0] * Ij[1] + Omegaj[j] / count2 * Ij[0] * Xj[1];
        ITensor H_N  = Deltaj[j-1] / count1 * Nj[0] * Ij[1] + Deltaj[j] / count2 * Ij[0] * Nj[1];
        ITensor H_NN = Vj[j-1] * Nj[0] * Nj[1];

        gates.push_back(TebdGate({j,j+1}, dt/2., H_NN + H_N + H_om));
    }
    return make_symmetric_sweep(gates);
}


// Rydberg Hamiltonian - we keep up to next-nearest neighbor interactions
// 1. We split H = H_1 + H_2 + H_3, so that [H_i,H_j] \neq 0 while the elements within each H_i commute.
// 2. Each H_i is made of three-site gates of time step dt/2; the sequence is [U_1,U_2,U_3,U_3,U_2,U_1].
// The single-site terms are not split off in separate gates: ~30% faster.

vector<TebdGate>
make_rydberg_gates_nnn(const SiteSet sites , const vector<double> Deltaj, const vector<double> Omegaj, const vector<double> Vj, const double dt)
{
    int N = length(sites);
    function<ITensor(int)> rydberg_term = [&](int j)
    {
        vector<ITensor> Nj, Ij, Xj;
        for(int q = j ; q <= j+2 ; q++)
        {
            Nj.push_back( (op(sites,"Id",q) - 2*op(sites,"Sz",q)) / 2. );
            Ij.push_back(  op(sites,"Id",q) );
            Xj.push_back(  op(sites,"Sx",q) );
        }

        // the next-nearest-neighbour interaction follows from the distances r = V^(-1/6)
        double r1  = pow(1/Vj[j-1], 1./6);
        double r2  = pow(1/Vj[j],   1./6);
        double V13 = pow(1/(r1+r2), 6.);

        // nearest-neighbour interactions and on-site terms are shared among the gates containing them
        double V12 = Vj[j-1] / count_gates_containing(j,   2, 3, N);
        double V23 = Vj[j]   / count_gates_containing(j+1, 2, 3, N);
        vector<double> omega, delta;
        for(int a = 0 ; a < 3 ; a++)
        {
            double count = count_gates_containing(j+a, 1, 3, N);
            omega.push_back(Omegaj[j-1+a] / count);
            delta.push_back(Deltaj[j-1+a] / count);
        }

        ITensor H_NN  = V12 * Nj[0] * Nj[1] * Ij[2];
        H_NN         += V23 * Ij[0] * Nj[1] * Nj[2];
        H_NN         += V13 * Nj[0] * Ij[1] * Nj[2];

        ITensor H1  = omega[0] * Xj[0] * Ij[1] * Ij[2];
        H1         += omega[1] * Ij[0] * Xj[1] * Ij[2];
        H1         += omega[2] * Ij[0] * Ij[1] * Xj[2];
        H1         += delta[0] * Nj[0] * Ij[1] * Ij[2];
        H1         += delta[1] * Ij[0] * Nj[1] * Ij[2];
        H1         += delta[2] * Ij[0] * Ij[1] * Nj[2];

        return H_NN + H1;
    };
    return make_three_site_layers(N, dt, rydberg_term);
}


// positions along x from the cyclic spacings, then the gaussian displacements (only if some sigma > 0,
// so that the random generator is not used otherwise); one unit normal distribution per direction,
// scaled by sigma, so that a zero sigma is allowed
vector<vector<double> >
make_chain_positions(const int N, const vector<double> spacings, const vector<double> sigma, const unsigned seed)
{
    vector<vector<double> > rj;
    double x = 0.;
    for(int j = 0 ; j < N ; j++)
    {
        rj.push_back({x, 0., 0.});
        x += spacings[j % spacings.size()];
    }
    if(sigma[0] > 0 || sigma[1] > 0 || sigma[2] > 0)
    {
        default_random_engine generator;
        generator.seed(seed);
        normal_distribution<double> noise_x(0, 1), noise_y(0, 1), noise_z(0, 1);
        for(vector<double>& r : rj)
        {
            r[0] += sigma[0] * noise_x(generator);
            r[1] += sigma[1] * noise_y(generator);
            r[2] += sigma[2] * noise_z(generator);
        }
    }
    return rj;
}


// V_j = 1/|r_j - r_{j+1}|^alpha from the Euclidean distance of consecutive atoms
 
vector<double>
compute_power_law_couplings(const vector< vector<double> > rj , const double alpha)
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
