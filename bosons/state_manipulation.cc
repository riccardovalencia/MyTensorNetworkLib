#include "state_manipulation.h"
#include <itensor/all.h>

using namespace std;
using namespace itensor;

// compute the reduced density matrix bewteen sites i and j
ITensor
exctract_reduced_density_matrix(MPS *psi, int i, int j)
{
    // SiteSet sites = siteInds(*psi);
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
