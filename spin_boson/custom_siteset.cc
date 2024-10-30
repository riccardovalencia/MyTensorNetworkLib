#include "TEBD.h"
#include <itensor/all.h>
#include <iostream>

// ----------------------------------------------------------
// Custom SiteSet for handling a system of N sites, where the first site is a (truncated)
// boson with max occupation max_occ, while the other (N-1) are spin-1/2

SiteSet
custom_spin_boson(const int N , const int max_occ)
{
    IndexSet is = IndexSet(N);

    TagSet ts = TagSet("Site,Boson");
    ts.addTags("n=1");
    is[0] = Index(max_occ+1,ts);
    // initialize all the other as spin-1/2
    for(int j=2; j<=N;j++)
    {
        ts = TagSet("Site,S=1/2");
        ts.addTags("n="+str(j));
        is[j-1] = Index(2,ts);
    }
    SiteSet sites = SiteSet(is);

    return sites;
    
}

// ----------------------------------------------------------
// return Siteset with a boson and N-1 spin-1/2 

// It is defined in the doubling space (ket-bra) and follows the ordering: 
// bra (first half of the chain) - ket (second half of the chain)
// The bra is inverted in space with respect to the ket. 

// Example for N = 4
//  |    |    |   |   |    |    |    |
// s3 - s2 - s1 - b - b - s1 - s2 - s2
// | ---- bra ---- |  | ---- ket ---- |

// Useful if you have dissipative channels acting solely on the bosonic DOF.

SiteSet
custom_spin_boson_doubling(const int N , const int max_occ)
{
    // doubling space (and so system size)
    int N2 = 2 * N;

    IndexSet is = IndexSet(N2);
    TagSet ts ; 
    for(int j=1; j<=N-1;j++)
    {
        ts = TagSet("Site,S=1/2");
        ts.addTags("n="+str(j));
        is[j-1] = Index(2,ts);
    }

    ts = TagSet("Site,Boson");
    ts.addTags("n="+str(N));
    is[N-1] = Index(max_occ+1,ts);

    ts = TagSet("Site,Boson");
    ts.addTags("n="+str(N+1));
    is[N] = Index(max_occ+1,ts);
    // initialize all the other as spin-1/2
    for(int j=N+2; j<=N2;j++)
    {
        ts = TagSet("Site,S=1/2");
        ts.addTags("n="+str(j));
        is[j-1] = Index(2,ts);
    }

    SiteSet sites = SiteSet(is);

    return sites;
    
}

