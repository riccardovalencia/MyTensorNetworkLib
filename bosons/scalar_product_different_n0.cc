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
scalar_product_different_n0( MPS *psi1, MPS *psi2, const SiteSet sites1, const SiteSet sites2, const int size, const int n0_1, const int n0_2 , const int lambda, const int symmetry_sector, const double s, const double c)
{
    // if you want to compute the variance over the Hamiltonian with smaller n_0 over the state with greater n_0
    int max_n0 = max(n0_1, n0_2);
    int min_n0 = min(n0_1, n0_2);
    
    // if you want to compute the variance over the Hamiltonian with greater n_0 over the state with smaller n_0
    // int max_n0 = min(n0_1, n0_2);
    // int min_n0 = max(n0_1, n0_2);

    vector<ITensor> expanded_state; 

    MPS psi_max_n0; 
    MPS psi_min_n0;
    Boson sites_max_n0;
    Boson sites_min_n0;

    Index leftindexj ;
    Index rightindexj ;
    Index physical ;
    Index physical_previous ;
    ITensor psi_tocopy_j ;
    ITensor Tj ;

    if (min_n0 == n0_1)
    {
        psi_min_n0 = *(psi1);
        psi_max_n0 = *(psi2);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_min_n0 = Boson(ind1);
        sites_max_n0 = Boson(ind2);
    }
    else
    {
        psi_min_n0 = *(psi2);
        psi_max_n0 = *(psi1);
        auto ind1 = inds(sites1);
        auto ind2 = inds(sites2);
        sites_min_n0 = Boson(ind2);
        sites_max_n0 = Boson(ind1);
    }
    

    // site 1

    rightindexj = rightLinkIndex( psi_max_n0, 1);
    physical = sites_min_n0(1);
    physical_previous = sites_max_n0(1);
    psi_tocopy_j = psi_max_n0(1);
    
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
        leftindexj  = leftLinkIndex( psi_max_n0, j);
        rightindexj = rightLinkIndex( psi_max_n0, j);
        physical = sites_min_n0(j);
        physical_previous = sites_max_n0(j);
        psi_tocopy_j = psi_max_n0(j);
    
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

    leftindexj  = leftLinkIndex( psi_max_n0, size);
    physical = sites_min_n0(size);
    physical_previous = sites_max_n0(size);
    psi_tocopy_j = psi_max_n0(size);
    
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


    // define a new MPS of the same kind of sites_max_n0 whose elements are set to be equal to the expanded_state
    // in this way I can use built-in functions in ITensor

    MPS psi_constrained = MPS(sites_min_n0);

    for(int j=1 ; j<=size ; j++) psi_constrained.set(j, expanded_state[j-1]);

    complex<double> overlap_amplitude = innerC(psi_min_n0,psi_constrained);
    double overlap_absolute = overlap_amplitude.real() * overlap_amplitude.real() + overlap_amplitude.imag() * overlap_amplitude.imag();

    // stringstream file_symmetry_name;
	// file_symmetry_name << "/home/ricval/Documenti/Bosonic/bQEM_Cpp_final/mm_symmetry_sector/symmetry_sector" << symmetry_sector << "_maxcutoff50_s" ;
	// file_symmetry_name << fixed << setprecision(2);
	// file_symmetry_name << s << "_c" << c << ".dat";

	// fstream file_symmetry(file_symmetry_name.str(),std::ios_base::in );
	// double symmetry ;
	// for(int j=1 ; j <= max(n0_1,cut_off_fock_space2) ; j++)
	// {
	// 	file_symmetry >> symmetry; 
	// }
	// file_symmetry.close();
    // MPO H_max_cutoff = H_mmGcbQEM( sites_max_n0, size , n0, symmetry , s, c);

    // double variance_expanded_space = inner(psi_expanded, H_max_cutoff, H_max_cutoff, psi_expanded) - inner(psi_expanded, H_max_cutoff, psi_expanded) * inner(psi_expanded, H_max_cutoff, psi_expanded); 

    cerr << "Computing variance of n0: " << max_n0 << " over Hamiltonian with n0=" << min_n0 << endl;
    double variance_constraned_space = compute_variance_H_mmGcbQEM(&psi_constrained, sites_min_n0 , size , lambda, min_n0, symmetry_sector, s, c);
    // double variance_constraned_space = compute_variance_H_mmGcbQEM(&psi_constrained, sites_min_n0 , size , lambda, 1, symmetry_sector, s, c);

    if( overlap_absolute < 1E-10 ) overlap_absolute = 0.;

    return {overlap_absolute,variance_constraned_space};
}


