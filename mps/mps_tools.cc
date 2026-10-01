/**
 * @file mps_tools.cc
 * @brief Implementation of mps_tools.h (interfaces documented in the header, logic commented here).
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


// Orthogonality center on i (after ordering i < j), then contract bra and ket from i to j with
// op_i and op_j inserted: the open link indices at the two ends are shared by bra and ket.
complex<double>
measure_two_point_function( MPS *psi, const SiteSet sites, ITensor op_i, ITensor op_j, int i, int j)
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

	return eltC(C);
}


// dense identity: element (q, q) = 1 for every basis state q

ITensor
make_identity_operator( const Index& in, const Index& out )
{
    ITensor I = ITensor(dag(in), out);
    for(int q = 1 ; q <= dim(in) ; q++) I.set(in(q), out(q), 1.);
    return I;
}


// with the orthogonality center on the site, <psi|O|psi> only involves the tensor of that site

Cplx
measure_local_operator( MPS* psi, const ITensor& O, const int site )
{
    (*psi).position(site);
    ITensor ket = (*psi)(site);
    ITensor bra = dag(prime(ket,"Site"));
    return eltC(bra * O * ket);
}


// Optionally reverse (and conjugate) psi_seed, retag its link indices (unique per insertion), then
// copy its tensors on start, start+1, ...: the first and last tensors get the external link of *psi
// as an extra index, filled only at value 1 (the seed is a product with the rest of *psi across it).
// The site indices of *psi are restored at the end.
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

        else throw ITError("insert_state: psi_seed does not fit in psi from site start");

    }
    (*psi).replaceSiteInds(phys_idx);
}


// ----------------------------------------------------------
// Copy the tensors of state_to_insert element by element onto the site indices of *psi_t0 (sites
// start..start+L-1): the link indices at the edges of the block are those of *psi_t0, the inner
// ones those of state_to_insert; sites beyond N are dropped.
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
			leftindexj  = leftLinkIndex( state_to_insert, j);
			rightindexj  = rightLinkIndex( state_to_insert, j);
			physical = sites(current_position_psi_t0);
			physical_state_to_insert = sites_state_to_insert(j);
			psi_tocopy_j = state_to_insert(j);
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
			
    for(int j=1 ; j<=L_effective ; j++) (*psi_t0).set(j+start-1, copy_of_state_to_insert[j-1]);	

}


// Move the site at j1 to j2 by exchanging neighbours j, j+1 for j = j1..j2-1: contract the two
// tensors and split them again with an SVD that puts the site index of j+1 on the left tensor
// (virtual swap, no gate; see http://itensor.org/support/2330/non-consecutive-swap-gates).
void
swap_sites( MPS *psi, int j1, int j2, double cut_off, int maxDim)
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


    if( j1 < 1 || j2 > length(sj)) throw ITError("swap_sites: sites out of the chain");

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
// Outer product |psi><psi| site by site, as in ITensor's nmultMPO: the pair of links (ket, bra)
// of each bond is fused into a single MPO link with a density-matrix decomposition (truncation
// 1E-16, at most 500 states), sweeping from left to right.
MPO 
make_density_matrix_mpo(MPS psi )
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
// Orthogonality center on i, then contract ket and primed bra on sites i..j; the open link
// indices at the two ends are shared by bra and ket, so the result is Tr_{rest} |psi><psi|.
ITensor
compute_reduced_density_matrix(MPS *psi, int i, int j)
{
    int L = length(*psi);

    if( i < 1 || i > L || j < 1 || j > L) throw ITError("compute_reduced_density_matrix: sites out of the chain");

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


// ----------------------------------------------------------
// new tensor on the site with the same indices (site, left and right links), non-zero only for
// link values 1: element d is amplitudes[d-1]

void
set_site_tensor( MPS* psi, const SiteSet& sites, int site, const vector<Cplx>& amplitudes )
{
    int N = length(*psi);
    Index s = sites(site);

    IndexSet is = {s};
    if(site > 1) is = IndexSet(is, leftLinkIndex(*psi, site));
    if(site < N) is = IndexSet(is, rightLinkIndex(*psi, site));
    ITensor wf = ITensor(is);

    bool real = true;
    for(Cplx a : amplitudes) if(a.imag() != 0.) real = false;

    for(int d = 1 ; d <= dim(s) ; d++)
    {
        // all link indices take value 1
        vector<IndexVal> iv = {s(d)};
        for(Index l : is) if(l != s) iv.push_back(l(1));
        if(real) wf.set(iv, amplitudes[d-1].real());
        else     wf.set(iv, amplitudes[d-1]);
    }
    (*psi).set(site, wf);
}


// ----------------------------------------------------------
// contract every site with the rectangular identity P(s, t) between its index and the target index

MPS
make_resized_state( MPS psi, const SiteSet& target_sites )
{
    for(int j = 1 ; j <= length(psi) ; j++)
    {
        Index s = siteIndex(psi, j);
        Index t = target_sites(j);
        ITensor P = ITensor(s, t);
        for(int d = 1 ; d <= min(dim(s), dim(t)) ; d++) P.set(s(d), t(d), 1.);
        psi.set(j, psi(j) * P);
    }
    return psi;
}
