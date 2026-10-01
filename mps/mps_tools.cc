/**
 * @file mps_tools.cc
 * @brief Implementation of mps_tools.h (the functions are documented in the header).
 */
#include "mps_tools.h"
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


complex<double>
compute_two_point( MPS *psi, const SiteSet sites, ITensor op_i, ITensor op_j, int i, int j)
{
	if(j<i)
	{
		int k = i;
		ITensor op_k = op_i;
		i = j;
		j = k;
		op_i = op_j;
		op_j = op_k;
	}

	//'gauge' the MPS to site i
	//any 'position' between i and j, inclusive, would work here
	(*psi).position(i); 

	//Create the bra/dual version of the MPS psi
	auto psidag = dag(*psi);

	//Prime the link indices to make them distinct from
	//the original ket links
	psidag.prime("Link");

	//index linking i-1 to i:
	auto li_1 = leftLinkIndex(*psi,i);

	auto C = prime((*psi)(i),li_1)*op_i;
	C *= prime(psidag(i),"Site");
	for(int k = i+1; k < j; ++k)
		{
		C *= (*psi)(k);
		C *= psidag(k);
		}
	//index linking j to j+1:
	auto lj = rightLinkIndex((*psi),j);

	C *= prime((*psi)(j),lj)*op_j;
	C *= prime(psidag(j),"Site");

	complex<double> result = eltC(C); //or eltC(C) if expecting complex	
	return result;
}


// Insert a state within another state, such that you have a state |state_to_insert> that you want to put in another state |psi_t0> from site start to start+L

void
insert_state(MPS* psi, MPS psi_seed, const int start, bool inverted,bool dagger)
{

    IndexSet phys_idx = siteInds(*psi);
    IndexSet phys_idx_seed = siteInds(psi_seed);

    IndexSet link_idx = linkInds(*psi);
    IndexSet link_idx_seed = linkInds(psi_seed);


    int N      = length(phys_idx);
    int N_seed = length(phys_idx_seed);

    Index lj ;
    Index rj ;
    Index sj ;
    ITensor Tj ;

    if(inverted)
    {
        MPS psi_seed_inv = psi_seed;
        int k = 1;
    
        for(int j= N_seed ; j >= 1; j--)
        {
            if(dagger) psi_seed_inv.set(k,dag(psi_seed(j)));
            else psi_seed_inv.set(k,psi_seed(j)); 
            k += 1;
            
        }

        psi_seed = psi_seed_inv;
        phys_idx_seed = siteInds(psi_seed);
        link_idx_seed = linkInds(psi_seed);
    }

    // add tags to link indices, in order to avoid issues if we insert multiple copyes of the same state
    IndexSet new_links = IndexSet(N_seed-1);
    for(int j : range1(N_seed-1))
	{
		int j_psi = start + j - 1;
        Index new_idx = addTags(link_idx_seed(j),"seed="+str(j_psi));
        new_links[j-1] = new_idx;
    }

    psi_seed.replaceLinkInds(new_links);
    link_idx_seed = linkInds(psi_seed);


	for(int j : range1(N_seed))
	{
		int j_psi = start + j - 1;


        if(j==1 && j_psi != 1)
        {

            // I have to add a right index
            rj = link_idx(j_psi-1);
            lj = link_idx_seed(j);
            sj = phys_idx_seed(j);
            Tj = ITensor(sj,rj,lj);

            for(int l=1; l<= dim(lj) ; l++)
			{
            for(int r=1; r<= dim(rj); r++)
            {
            for(int d=1; d<= dim(sj); d++)
            {
                if(r==1) Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, lj=l) );
                else Tj.set(lj=l,sj=d,rj=r , 0 );
            }
            }
			}

            (*psi).set(j_psi,Tj);

        }

        else if (j==N_seed && j_psi != N)
        {
            // I have to add a left index
            rj = link_idx_seed(j-1);
            lj = link_idx(j_psi);
            sj = phys_idx_seed(j);
            Tj = ITensor(sj,rj,lj);

            for(int l=1; l<= dim(lj) ; l++)
			{
            for(int r=1; r<= dim(rj); r++)
            {
            for(int d=1; d<= dim(sj); d++)
            {
                if(l==1) Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, rj=r) );
                else Tj.set(lj=l,sj=d,rj=r , 0 );		
            }
            }
			}

            (*psi).set(j_psi,Tj);
        }

        else if(j_psi <= N)
        {
            (*psi).set(j_psi,psi_seed(j));
        }

        else
        {
            cerr << "There is no room for accomodating the state." << endl;
            break;
        }

    }
    (*psi).replaceSiteInds(phys_idx);
}


// // Insert a state within another state, such that you have a state |state_to_insert> that you want to put in another state |psi_t0> from site start to start+L
// As we are dealing with sitesets with QN quantities, we have to define the link indices with a flux direction (In or Out)
// -> psi(j) ->

void
insert_qn_state(MPS* psi, MPS psi_seed, const int start, bool inverted,bool dagger)
{

    IndexSet phys_idx = siteInds(*psi);
    IndexSet phys_idx_seed = siteInds(psi_seed);

    IndexSet link_idx = linkInds(*psi);
    IndexSet link_idx_seed = linkInds(psi_seed);

    int N      = length(phys_idx);
    int N_seed = length(phys_idx_seed);

    Index lj ;
    Index rj ;
    Index sj ;
    ITensor Tj ;

    if(inverted)
    {
        MPS psi_seed_inv = psi_seed;
        int k = 1;
    
        for(int j= N_seed ; j >= 1; j--)
        {
            if(dagger) psi_seed_inv.set(k,dag(psi_seed(j)));
            else psi_seed_inv.set(k,psi_seed(j)); 
            k += 1;
            
        }

        psi_seed = psi_seed_inv;
        psi_seed.replaceSiteInds(phys_idx_seed);
        psi_seed.replaceLinkInds(link_idx_seed);
        link_idx_seed = linkInds(psi_seed);
    }

    // add tags to link indices, in order to avoid issues if we insert multiple copyes of the same state
    IndexSet new_links = IndexSet(N_seed-1);
    for(int j : range1(N_seed-1))
	{
		int j_psi = start + j - 1;
        Index new_idx = addTags(link_idx_seed(j),"seed="+str(j_psi));
        new_links[j-1] = new_idx;
    }

    psi_seed.replaceLinkInds(new_links);
    link_idx_seed = linkInds(psi_seed);


	for(int j : range1(N_seed))
	{
		int j_psi = start + j - 1;


        if(j==1 && j_psi != 1)
        {

            // I have to add a right index
            rj = link_idx(j_psi-1);
            lj = link_idx_seed(j);
            sj = phys_idx_seed(j);
            Tj = ITensor(sj,rj,lj);

            for(int l=1; l<= dim(lj) ; l++)
			{
            for(int r=1; r<= dim(rj); r++)
            {
            for(int d=1; d<= dim(sj); d++)
            {
                if(r==1) Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, lj=l) );
                else Tj.set(lj=l,sj=d,rj=r , 0 );
            }
            }
			}

            (*psi).set(j_psi,Tj);

        }

        else if (j==N_seed && j_psi != N)
        {
            // I have to add a left index
            rj = link_idx_seed(j-1);
            lj = link_idx(j_psi);
            sj = phys_idx_seed(j);
            Tj = ITensor(sj,rj,lj);

            cerr << Tj << endl;
            exit(0);
            for(int l=1; l<= dim(lj) ; l++)
			{
            for(int r=1; r<= dim(rj); r++)
            {
            for(int d=1; d<= dim(sj); d++)
            {
                if(l==1) Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, rj=r) );
                else Tj.set(lj=l,sj=d,rj=r , 0 );		
            }
            }
			}

            (*psi).set(j_psi,Tj);
        }

        else if(j_psi <= N)
        {
            (*psi).set(j_psi,psi_seed(j));
            cerr << "I have inserted the state on site : " << j_psi << "\n";
        }

        else
        {
            cerr << "There is no room for accomodating the state." << endl;
            break;
        }

    }
    (*psi).replaceSiteInds(phys_idx);
}


// ----------------------------------------------------------
// Insert a state within another state, such that you have a state |state_to_insert> that you want to put in another state |psi_t0> from site start to start+L

void
insert_state(MPS* psi_t0, MPS state_to_insert, const SiteSet sites, const SiteSet sites_state_to_insert, const int start, const int L, const int N)
{
    vector<ITensor> copy_of_state_to_insert; 


    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
	Index physical_state_to_insert;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

	int L_effective = 0;

	for(int j=1 ; j<=L ; j++)
	{
		L_effective += 1;
		int current_position_psi_t0 = start + j -1;
		if(current_position_psi_t0<N)
		{
			if(j==1)     leftindexj  = leftLinkIndex(  *psi_t0, current_position_psi_t0);
			else         leftindexj  = leftLinkIndex(  state_to_insert, j);
			if(j==L)     rightindexj = rightLinkIndex( *psi_t0, current_position_psi_t0);
			else         rightindexj = rightLinkIndex( state_to_insert, j);
			physical    = sites(current_position_psi_t0);
			physical_state_to_insert    = sites_state_to_insert(j);
			psi_tocopy_j = state_to_insert(j);

			Tj = ITensor(leftindexj, physical, rightindexj);
			
			for(int l=1; l<= dim(leftindexj) ; l++)
			{
				for(int r=1; r<= dim(rightindexj); r++)
				{
					for(int d=1; d<= dim(physical); d++)
					{
						if(j!=1 && j!=L) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, leftindexj=l,physical_state_to_insert=d, rightindexj=r) );
						else if(j==1) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_state_to_insert=d, rightindexj=r) );			
						else if(j==L) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_state_to_insert=d, leftindexj=l) );
					}
				}
			}

			copy_of_state_to_insert.push_back(Tj);
		}
		
		else if(current_position_psi_t0==N)
		{
			cerr << "Siamo a : " << current_position_psi_t0 << endl;
			leftindexj  = leftLinkIndex( state_to_insert, j);
			rightindexj  = rightLinkIndex( state_to_insert, j);
			physical = sites(current_position_psi_t0);
			physical_state_to_insert = sites_state_to_insert(j);
			psi_tocopy_j = state_to_insert(j);
			cerr << psi_tocopy_j << endl;
			Tj = ITensor(leftindexj, physical);
			
			for(int l=1; l<= dim(leftindexj) ; l++)
			{
					for(int d=1; d<= dim(physical); d++)
					{
						Tj.set(leftindexj=l,physical=d , eltC(psi_tocopy_j, leftindexj=l,physical_state_to_insert=d, rightindexj=1));					
					}
			}

			copy_of_state_to_insert.push_back(Tj);
		}

		else break;
	}
	cerr << "L effective : " << L_effective << endl;	
			
    for(int j=1 ; j<=L_effective ; j++) (*psi_t0).set(j+start-1, copy_of_state_to_insert[j-1]);	

}


// move index from position j1 to position j2
// implementing virtual swap (not applying a gate)
// see http://itensor.org/support/2330/non-consecutive-swap-gates

void
swap_gate( MPS *psi, int j1, int j2, double cut_off, int maxDim)
{

    SiteSet sj = siteInds(*psi);
    ITensor T, U, V, D;
    Index sj_p;
    Index lj;


    if(j1 > j2)
    {
        int k = j1;
        j1 = j2;
        j2 = k; 
    }


    if( j1 < 1 || j2 > length(sj))
    {
        cerr << "Index out of physical bound" << endl;
        exit(0);
    }

    for(int j=j1 ; j<j2; j++)
    {
        T = (*psi)(j) * (*psi)(j+1);
        sj_p = sj(j+1);
        if (j > 1){
            lj = commonIndex((*psi)(j-1),(*psi)(j));
        }
        else{
            // NOT SURE. THINK ABOUT IT
            lj = commonIndex((*psi)(j),(*psi)(j+1));
            }

        U = ITensor(sj_p,lj);
        svd(T,U,D,V,{"Cutoff=",cut_off,"MaxDim=",maxDim});
        (*psi).set(j,U);
        (*psi).set(j+1,D*V);

    }


}


// ----------------------------------------------------------
// Given a pure state psi, presented as an MPS, it return its density matrix representation |psi> <psi| as an MPO

MPO 
from_mps_to_mpdo(MPS psi )
{
    // ket and bra (bra is primed)
    MPS ket = psi;
    MPS bra = dag(prime(psi));

    // In order to merge the two MPS into an MPO (outer product), I use a similar procedure used in
    // nmultMPO -> the idea is to perform a transformation so that we introduce a new link index.

    // IndexSet 
    IndexSet sA  = siteInds(psi);
    IndexSet sB  = siteInds(bra);

    int N = length(sA);

    // MPO hosting the final MatrixProductDensityOperator
    MPO rho = MPO(sA);
    if(N==1)
    {
        rho.ref(1) = psi(1) * bra(1);
    }

    else
    {    
        IndexSet lA = linkInds(ket);
        IndexSet lB = linkInds(bra);

        // length

        ITensor clust, nfork;

        rho.ref(1) = ITensor(sA(1),sB(1),lA(1));

        for(int i : range1(N))
        {
            if(i==1) clust = psi(i) * bra(i);
            else     clust = nfork * psi(i) * bra(i);
            if(i==N-1) break;   

            // for site i=1 -> we have Tensor 
            //      |
            //      o =
            //      |
            // We want a single link index. To do so, we cut along the two link indices
            //      |
            //      o -  -o<
            //      | 
            // we need a 'fork' Tensor: -o< tensor

            nfork = ITensor(lA(i),lB(i),linkIndex(rho,i));

            // it perform a denmatDecomp (similar to SVD): will cut along the two link indeces
            denmatDecomp(clust,rho.ref(i),nfork,Fromleft,{{"MaxDim",500,"Cutoff",1E-16},"Tags=",tags(linkIndex(rho,i))});

            Index mid = commonIndex(rho(i),nfork);
            mid.dag(); // why the dag? I am taking from the nmultMPO ITensor code
            rho.ref(i+1) = ITensor(mid,sA(i+1),sB(i+1),rightLinkIndex(rho,i+1));
        }

        nfork = clust * psi(N) * bra(N);
        rho.svdBond(N-1,nfork,Fromright, {"MaxDim",500,"Cutoff",1E-16});
        rho.orthogonalize();
    }
    return rho;


}


// ----------------------------------------------------------
// Given a pure state psi, presented as an MPS, it return its density matrix representation |psi> <psi| as an MPO
// Similar to above, but it uses a variation for fusing the link indices in a single one

MPO 
from_mps_to_mpdo_v2(MPS psi )
{
    // ket and bra (bra is primed)
    MPS ket = psi;
    MPS bra = dag(prime(psi));

    // In order to merge the two MPS into an MPO (outer product), I use a similar procedure used in
    // nmultMPO -> the idea is to perform a transformation so that we introduce a new link index.

    // IndexSet 
    IndexSet sA  = siteInds(psi);
    IndexSet sB  = siteInds(bra);

    int N = length(sA);
    // MPO hosting the final MatrixProductDensityOperator
    MPO rho = MPO(sA);
    if(N==1)
    {
        rho.ref(1) = psi(1) * bra(1);
    }

    else
    {    
        IndexSet lA = linkInds(ket);
        IndexSet lB = linkInds(bra);
        IndexSet lrho = linkInds(rho);

        ITensor clust, nfork; // helper ITensor


        for(int i = 1 ; i <= N ; i++)
        {
            clust = psi(i) * bra(i);
            if(i==1)
            {
                auto [C,c] = combiner(lA(i),lB(i));
                clust = clust * C;
                clust *= delta(c,lrho(i));
            }

            else if(i>1 && i < N)
            {
                auto [C,c] = combiner(lA(i-1),lB(i-1));
                clust = clust * C;
                auto [C2,c2] = combiner(lA(i),lB(i));
                clust = clust * C2;
                   
                clust *= delta(c,lrho(i-1));
                clust *= delta(c2,lrho(i));                
            }

            else
            {
                auto [C,c] = combiner(lA(i-1),lB(i-1));
                clust = clust * C;
                clust *= delta(c,lrho(i-1));
            }
            rho.ref(i) = clust;
            cerr << rho(i) << endl;
        }
    }
    return rho;


}


// ----------------------------------------------------------
// compute the reduced density matrix bewteen sites i and j
ITensor
extract_reduced_density_matrix(MPS *psi, int i, int j)
{
    int L = length(*psi);

    if( i < 1 || i > L || j < 1 || j > L){
        cerr << "Invalid set on which compute the reduced density matrix" << endl;
        exit(-1);
    }

    if( j < i){
        int tmp = j;
        j = i;
        i = tmp;
    }

    //'gauge' the MPS to site j
    //any 'position' between j and i, inclusive, would work here
    ITensor rho ;

    if( i != j)
    {
        (*psi).position(i); 
        MPS psidag = dag((*psi));
        psidag.prime("Link");


        Index li_1 = leftLinkIndex((*psi),i);

        rho = prime((*psi)(i),li_1)*prime(psidag(i),"Site");
        for(int k = i+1; k < j; ++k)
            {
            rho *= ((*psi)(k)) * prime(psidag(k),"Site");
            }
        //index linking i to i+1:
        Index lj = rightLinkIndex((*psi),j);
        rho *= prime((*psi)(j),lj);
        rho *= prime(psidag(j),"Site");
    }

    else
    {
        (*psi).position(i); 
        MPS psidag = dag((*psi));
        rho = (*psi)(i)*prime(psidag(i),"Site");     
    }
    return rho;

}
