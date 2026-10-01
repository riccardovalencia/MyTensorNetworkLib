/**
 * @file spin_boson.cc
 * @brief Implementation of spin_boson.h (interfaces documented in the header, logic commented here).
 */
#include "spin_boson.h"
#include "../mps/mps_tools.h"
#include <itensor/all.h>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// One index per site: dimension max_occ + 1 with tags "Site,Boson" on site 1, dimension 2 with
// tags "Site,S=1/2" on sites 2..N (tag n=j as in the ITensor site sets).

SiteSet
make_spin_boson_sites(const int N , const int max_occ)
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
// As make_spin_boson_sites on 2N sites, with the bra (mirrored) on 1..N and the ket on N+1..2N,
// so that the two bosons are on the central bond, where the dissipation acts. For N = 4:
//  |    |    |    |   |    |    |    |
// s3 - s2 - s1 - b - b - s1 - s2 - s3
// | ---- bra ---- |  | ---- ket ---- |

SiteSet
make_purified_spin_boson_sites(const int N , const int max_occ)
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


// ----------------------------------------------------------
// Product state: Fock state |n_photon> on boson sites, spin-coherent state
// |theta_j,phi_j> = cos(theta_j/2)|up_z> + e^{i phi_j} sin(theta_j/2)|down_z> on spin site j.
// The tensors are real when phi_j = 0. Each site is set with set_site_tensor.

MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , const vector<double> theta, const vector<double> phi)
{
    MPS psi = MPS(sites);
    for(int j = 1 ; j <= length(sites) ; j++)
    {
        Index sj = sites(j);
        vector<Cplx> amplitudes;
        if(hasTags(sj,"Site,Boson"))
        {
            amplitudes.assign(dim(sj), 0.);
            amplitudes[n_photon] = 1.;
        }
        else if(hasTags(sj,"Site,S=1/2"))
        {
            Cplx phase = (phi[j-1] == 0.) ? Cplx(1.) : exp(Cplx_i*phi[j-1]);
            amplitudes = {cos(theta[j-1]/2.), sin(theta[j-1]/2.) * phase};
        }
        else throw ITError("make_spin_boson_state: sites must be bosons or spins 1/2");
        set_site_tensor(&psi, sites, j, amplitudes);
    }
    return psi;
}


// the same angles on every site
MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , double theta, double phi)
{
    int N = length(sites);
    return make_spin_boson_state(sites, n_photon, vector<double>(N, theta), vector<double>(N, phi));
}
