/**
 * @file spin_boson.cc
 * @brief Implementation of spin_boson.h (the functions are documented in the header).
 */
#include "spin_boson.h"
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


// ----------------------------------------------------------
// Custom SiteSet for handling a system of N sites, where the first site is a (truncated)
// boson with max occupation max_occ, while the other (N-1) are spin-1/2

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


// Amplitude sin(theta/2) e^{i phi} of |down_z> in the spin-coherent state
// |theta,phi> = cos(theta/2)|up_z> + e^{i phi} sin(theta/2)|down_z>.
// The tensor is kept real when phi = 0.
template<typename... IndexVals>
static void
set_down_amplitude(ITensor& wf, double theta, double phi, IndexVals... iv)
{
    if(phi == 0.) wf.set(iv..., sin(theta/2.));
    else          wf.set(iv..., sin(theta/2.)*exp(Cplx_i*phi));
}


// ----------------------------------------------------------
// Given a mixed spin_boson SiteSet, initialize product state |n_photon > \otimes |\theta,\phi>^(N-1)
// where |n_photon> is a Fock state
//       |theta,\phi> is a spin-coherent state of a spin-1/2 pointing on the Bloch sphere 
// 
// TO DO : do a function ITensor spin_coherent_state(IndexSet,theta,phi) which returns a spin coherent state
//       : do a function ITensor fock_state(IndexSet,n) which return |n>

MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , double theta, double phi)
{

    MPS psi = randomMPS(sites);
    int N   = length(sites);

    Index sj ,rj, lj;

    // first site
	sj = sites(1);
	rj = commonIndex(psi(1),psi(2));
	ITensor wf = ITensor(sj,rj);
	
    if(hasTags(sj,"Site,Boson"))
    {
        for( int d=1; d <= n_photon; d++) wf.set(sj(d),rj(1), 0);
        wf.set(sj(n_photon+1),rj(1),1);
        for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),rj(1), 0);
    }
				
    else if(hasTags(sj,"Site,S=1/2"))
    {
        cerr << "Inserting spin coherent state" << endl;
        wf.set(sj(1),rj(1),cos(theta/2.));
        set_down_amplitude(wf, theta, phi, sj(2),rj(1));
    }
    
    else{
        cerr << "SiteSet not recognize : " << sj << endl;
        cerr << "Return a random initial state" << endl;
        return psi; 
    }

	psi.set(1,wf);
	
    cerr << "Inserted spin coherent state" << endl;

    for(int j=2 ; j < N ; j++)
    {
        sj = sites(j);
        lj = commonIndex(psi(j-1),psi(j));
		rj = commonIndex(psi(j),psi(j+1));
		wf = ITensor(sj,lj,rj);

        if(hasTags(sj,"Site,Boson"))
        {
            for( int d=1; d <= n_photon; d++) wf.set(sj(d),rj(1),lj(1), 0);
            wf.set(sj(n_photon+1),rj(1),lj(1),1);
            for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),rj(1),lj(1), 0);
        }
                    
        else if(hasTags(sj,"Site,S=1/2"))
        {
            wf.set(sj(1),lj(1),rj(1),cos(theta/2.));
            set_down_amplitude(wf, theta, phi, sj(2),lj(1),rj(1));
        }

		psi.set(j,wf); 
    }

    sj = sites(N);
    lj = commonIndex(psi(N-1),psi(N));
    wf = ITensor(sj,lj);

    wf.set(sj(1),lj(1),cos(theta/2.));
    set_down_amplitude(wf, theta, phi, sj(2),lj(1));


    if(hasTags(sj,"Site,Boson"))
    {
        for( int d=1; d <= n_photon; d++) wf.set(sj(d),lj(1), 0);
        wf.set(sj(n_photon+1),lj(1),1);
        for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),lj(1), 0);
    }
                    
    else if(hasTags(sj,"Site,S=1/2"))
    {
        wf.set(sj(1),lj(1),cos(theta/2.));
        set_down_amplitude(wf, theta, phi, sj(2),lj(1));
    }


    psi.set(N,wf); 

    return psi;
}


// ----------------------------------------------------------
// Override previous function - it can accept vectors as theta and phi
// Given a mixed spin_boson SiteSet, initialize product state |n_photon > \otimes |\theta,\phi>^(N-1)
// where |n_photon> is a Fock state
//       |theta,\phi> is a spin-coherent state of a spin-1/2 pointing on the Bloch sphere 
// 
// TO DO : do a function ITensor spin_coherent_state(IndexSet,theta,phi) which returns a spin coherent state
//       : do a function ITensor fock_state(IndexSet,n) which return |n>
MPS
make_spin_boson_state(const SiteSet sites , const int n_photon , const vector<double> theta, const vector<double> phi)
{

    MPS psi = randomMPS(sites);
    int N   = length(sites);

    if(N==1)
    {
        Index sj = sites(1);
        ITensor wf = ITensor(sj);
        
        if(hasTags(sj,"Site,Boson"))
        {
            for( int d=1; d <= n_photon; d++) wf.set(sj(d), 0);
            wf.set(sj(n_photon+1),1);
            for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d), 0);
        }
                    
        else if(hasTags(sj,"Site,S=1/2"))
        {
            cerr << "Inserting spin coherent state" << endl;
            wf.set(sj(1),cos(theta[0]/2.));
            set_down_amplitude(wf, theta[0], phi[0], sj(2));
        }
        
        else{
            cerr << "SiteSet not recognize : " << sj << endl;
            cerr << "Return a random initial state" << endl;
            return psi; 
        }

        psi.set(1,wf);
    }

    if(N>1)
    {

        Index sj ,rj, lj;

        // first site
        sj = sites(1);
        rj = commonIndex(psi(1),psi(2));
        ITensor wf = ITensor(sj,rj);
        
        if(hasTags(sj,"Site,Boson"))
        {
            for( int d=1; d <= n_photon; d++) wf.set(sj(d),rj(1), 0);
            wf.set(sj(n_photon+1),rj(1),1);
            for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),rj(1), 0);
        }
                    
        else if(hasTags(sj,"Site,S=1/2"))
        {
            cerr << "Inserting spin coherent state" << endl;
            wf.set(sj(1),rj(1),cos(theta[0]/2.));
            set_down_amplitude(wf, theta[0], phi[0], sj(2),rj(1));
        }
        
        else{
            cerr << "SiteSet not recognize : " << sj << endl;
            cerr << "Return a random initial state" << endl;
            return psi; 
        }

        psi.set(1,wf);
        
        cerr << "Inserted spin coherent state" << endl;

        for(int j=2 ; j < N ; j++)
        {
            sj = sites(j);
            lj = commonIndex(psi(j-1),psi(j));
            rj = commonIndex(psi(j),psi(j+1));
            wf = ITensor(sj,lj,rj);

            if(hasTags(sj,"Site,Boson"))
            {
                for( int d=1; d <= n_photon; d++) wf.set(sj(d),rj(1),lj(1), 0);
                wf.set(sj(n_photon+1),rj(1),lj(1),1);
                for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),rj(1),lj(1), 0);
            }
                        
            else if(hasTags(sj,"Site,S=1/2"))
            {
                wf.set(sj(1),lj(1),rj(1),cos(theta[j-1]/2.));
                set_down_amplitude(wf, theta[j-1], phi[j-1], sj(2),lj(1),rj(1));
            }

            psi.set(j,wf); 
        }

        sj = sites(N);
        lj = commonIndex(psi(N-1),psi(N));
        wf = ITensor(sj,lj);


        if(hasTags(sj,"Site,Boson"))
        {
            for( int d=1; d <= n_photon; d++) wf.set(sj(d),lj(1), 0);
            wf.set(sj(n_photon+1),lj(1),1);
            for( int d=n_photon+2; d <= dim(sj); d++) wf.set(sj(d),lj(1), 0);
        }
                        
        else if(hasTags(sj,"Site,S=1/2"))
        {
            wf.set(sj(1),lj(1),cos(theta[N-1]/2.));
            set_down_amplitude(wf, theta[N-1], phi[N-1], sj(2),lj(1));
        }


        psi.set(N,wf); 
    }

    return psi;
}
