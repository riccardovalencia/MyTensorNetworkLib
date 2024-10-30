#include <itensor/all.h>
#include "build_hamiltonian.h"
#include "observables.h"
#include <vector>
#include <fstream>	//output file
#include <iomanip>
#include <iostream>
#include <math.h>
#include <complex>

using namespace std;
using namespace itensor;

tuple<double, double>
scalar_product_different_cutoff( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int cut_off_fock_space1, const int cut_off_fock_space2 , const int n0, const int symmetry_sector, const double s, const double c)
{

    int max_cut_off = max(cut_off_fock_space1, cut_off_fock_space2);

    vector<ITensor> expanded_state; 

    MPS psi_maxcutoff; 
    MPS psi_mincutoff;
    Boson sites_maxcutoff;
    Boson sites_mincutoff;

    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
    Index physical_previous ;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

    if (max_cut_off == cut_off_fock_space1)
    {
        psi_maxcutoff = *(psi1);
        psi_mincutoff = *(psi2);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_maxcutoff = Boson(ind1);
        sites_mincutoff = Boson(ind2);
    }
    else
    {
        psi_maxcutoff = *(psi2);
        psi_mincutoff = *(psi1);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_maxcutoff = Boson(ind2);
        sites_mincutoff = Boson(ind1);
    }

    // site 1

    rightindexj = rightLinkIndex( psi_mincutoff, 1);
    physical = sites_maxcutoff(1);
    physical_previous = sites_mincutoff(1);
    psi_tocopy_j = psi_mincutoff(1);
    
    Tj = ITensor(physical, rightindexj);

    for(int r=1; r<= dim(rightindexj); r++)
    {
        for(int d=1; d<= dim(physical); d++)
        {
            if( d <= dim(physical_previous) ) Tj.set(physical=d,rightindexj=r , eltC(psi_tocopy_j, physical_previous=d, rightindexj=r) );
            else Tj.set(physical=d,rightindexj=r , 0);
        }
    }
    
    expanded_state.push_back(Tj);

    // in the middle of the chain

    for (int j =2 ; j<= size-1 ; j ++)
    {  
        leftindexj  = leftLinkIndex( psi_mincutoff, j);
        rightindexj = rightLinkIndex( psi_mincutoff, j);
        physical = sites_maxcutoff(j);
        physical_previous = sites_mincutoff(j);
        psi_tocopy_j = psi_mincutoff(j);
    
        Tj = ITensor(rightindexj, physical, leftindexj);


        for(int l=1; l<= dim(leftindexj) ; l++)
        {
            for(int r=1; r<= dim(rightindexj); r++)
            {
                for(int d=1; d<= dim(physical); d++)
                {
                    if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d,rightindexj=r , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d, rightindexj=r) );
                    else Tj.set(leftindexj=l,physical=d, rightindexj=r, 0);
                
                }
            }
        }
 
        expanded_state.push_back(Tj);

    }

    // last site

    leftindexj  = leftLinkIndex( psi_mincutoff, size);
    physical = sites_maxcutoff(size);
    physical_previous = sites_mincutoff(size);
    psi_tocopy_j = psi_mincutoff(size);
    
    Tj = ITensor(leftindexj, physical);
    
    for(int l=1; l<= dim(leftindexj) ; l++)
    {
            for(int d=1; d<= dim(physical); d++)
            {
                if( d <= dim(physical_previous) ) Tj.set(leftindexj=l,physical=d , eltC(psi_tocopy_j, leftindexj=l,physical_previous=d));
                else Tj.set(leftindexj=l,physical=d , 0);
            
            }
    }

    expanded_state.push_back(Tj);

    // scalar product home made (NOT USING ITENSOR INNER())

    // ITensor overlap;

    // MPS psi_maxcutoff_dag = dag(psi_maxcutoff);

    // psi_maxcutoff_dag.replaceLinkInds(sim(linkInds(psi_maxcutoff_dag)));

 
    // overlap = expanded_state[0] * psi_maxcutoff_dag(1) ;

    // for( int j = 2 ; j <= size ; j++)
    // {
    //     overlap *= expanded_state[j-1] * psi_maxcutoff_dag(j) ;
    // }

    // complex<double> overlap_amplitude = eltC(overlap);

    // define a new MPS of the same kind of sites_maxcutoff whose elements are set to be equal to the expanded_state
    // in this way I can use built-in functions in ITensor

    MPS psi_expanded = MPS(sites_maxcutoff);

    for(int j=1 ; j<=size ; j++) psi_expanded.set(j, expanded_state[j-1]);

    complex<double> overlap_amplitude = innerC(psi_maxcutoff,psi_expanded);
    double overlap_absolute = overlap_amplitude.real() * overlap_amplitude.real() + overlap_amplitude.imag() * overlap_amplitude.imag();

    // stringstream file_symmetry_name;
	// file_symmetry_name << "/home/ricval/Documenti/Bosonic/bQEM_Cpp_final/mm_symmetry_sector/symmetry_sector" << symmetry_sector << "_maxcutoff50_s" ;
	// file_symmetry_name << fixed << setprecision(2);
	// file_symmetry_name << s << "_c" << c << ".dat";

	// fstream file_symmetry(file_symmetry_name.str(),std::ios_base::in );
	// double symmetry ;
	// for(int j=1 ; j <= max(cut_off_fock_space1,cut_off_fock_space2) ; j++)
	// {
	// 	file_symmetry >> symmetry; 
	// }
	// file_symmetry.close();
    // MPO H_max_cutoff = H_mmGcbQEM( sites_maxcutoff, size , n0, symmetry , s, c);

    // double variance_expanded_space = inner(psi_expanded, H_max_cutoff, H_max_cutoff, psi_expanded) - inner(psi_expanded, H_max_cutoff, psi_expanded) * inner(psi_expanded, H_max_cutoff, psi_expanded); 

    double variance_expanded_space = compute_variance_H_mmGcbQEM(&psi_expanded, sites_maxcutoff , size , max_cut_off, n0, symmetry_sector, s, c);

    // if( overlap_absolute < 1E-10 ) overlap_absolute = 0.;

    return {overlap_absolute,variance_expanded_space};
}


