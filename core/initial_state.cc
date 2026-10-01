#include "initial_state.h"
#include <itensor/all.h>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace itensor;


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
        // test - no need to replace siteinds in the inverted one
        // psi_seed.replaceSiteInds(phys_idx_seed);
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


    // cerr << psi_seed(2) << endl;
    // cerr << link_idx_seed << endl;
    // double a ;
    // cin >> a;

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
                // Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, rj=r) );	
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
insert_QN_state(MPS* psi, MPS psi_seed, const int start, bool inverted,bool dagger)
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


    // cerr << psi_seed(2) << endl;
    // cerr << link_idx_seed << endl;
    // double a ;
    // cin >> a;

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
                // Tj.set(lj=l,sj=d,rj=r , eltC(psi_seed(j), sj=d, rj=r) );	
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

// ----------------------------------------------------------
// Product state of spin-1/2, |psi> = |c_1> |c_2> ... |c_N>, with c_j = config[j-1] in {'0','1'}:
//   basis "z": '0' -> |up_z>, '1' -> |down_z>
//   basis "x": '0' -> |+x> = (|up_z> + |down_z>)/sqrt(2), '1' -> |-x> = (|up_z> - |down_z>)/sqrt(2)
//   basis "y": '0' -> |+y> = (|up_z> + i|down_z>)/sqrt(2), '1' -> |-y> = (|up_z> - i|down_z>)/sqrt(2)
// Examples: "0000" (all up), "1111" (all down), "0011" (domain wall), "0101" (Neel).

MPS
initial_computational_state(const SiteSet sites , const string config , const string basis)
{
    int N = length(sites);

    if((int)config.size() != N)
        throw ITError(tinyformat::format("initial_computational_state: config \"%s\" has %d characters, but there are %d sites",config,config.size(),N));
    if(basis != "z" && basis != "x" && basis != "y")
        throw ITError("initial_computational_state: basis must be \"z\", \"x\" or \"y\", got \"" + basis + "\"");

    // product state with link indices of dimension 1
    MPS psi(sites);

    for(int j : range1(N))
    {
        Index sj = sites(j);
        char c = config[j-1];

        if(!hasTags(sj,"Site,S=1/2"))
            throw ITError(tinyformat::format("initial_computational_state: site %d is not a spin-1/2",j));
        if(c != '0' && c != '1')
            throw ITError(tinyformat::format("initial_computational_state: invalid character '%c' in config (only '0' and '1' allowed)",c));

        // amplitudes on |up_z> and |down_z>
        double sign = (c == '0') ? 1. : -1.;
        Cplx a, b;
        if(basis == "z")
        {
            a = (c == '0') ? 1. : 0.;
            b = (c == '0') ? 0. : 1.;
        }
        else if(basis == "x")
        {
            a = 1./sqrt(2.);
            b = sign/sqrt(2.);
        }
        else
        {
            a = 1./sqrt(2.);
            b = sign*Cplx_i/sqrt(2.);
        }

        IndexSet is = {sj};
        if(j > 1) is = IndexSet(is, leftLinkIndex(psi,j));
        if(j < N) is = IndexSet(is, rightLinkIndex(psi,j));
        ITensor wf = ITensor(is);

        // all link indices take value 1
        vector<IndexVal> up = {sj(1)}, dn = {sj(2)};
        for(Index l : is) if(l != sj) { up.push_back(l(1)); dn.push_back(l(1)); }

        if(basis == "y")
        {
            wf.set(up, a);
            wf.set(dn, b);
        }
        else
        {
            wf.set(up, a.real());
            wf.set(dn, b.real());
        }
        psi.set(j, wf);
    }

    psi.position(1);
    return psi;
}
