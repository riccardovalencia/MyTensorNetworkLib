/**
 * @file spin_half.cc
 * @brief Implementation of spin_half.h (the functions are documented in the header).
 */
#include "spin_half.h"
#include "../mps/mps_tools.h"
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
// Product state of spin-1/2, |psi> = |c_1> |c_2> ... |c_N>, with c_j = config[j-1] in {'0','1'}:
//   basis "z": '0' -> |up_z>, '1' -> |down_z>
//   basis "x": '0' -> |+x> = (|up_z> + |down_z>)/sqrt(2), '1' -> |-x> = (|up_z> - |down_z>)/sqrt(2)
//   basis "y": '0' -> |+y> = (|up_z> + i|down_z>)/sqrt(2), '1' -> |-y> = (|up_z> - i|down_z>)/sqrt(2)
// Examples: "0000" (all up), "1111" (all down), "0011" (domain wall), "0101" (Neel).

MPS
make_product_state(const SiteSet sites , const string config , const string basis)
{
    int N = length(sites);

    if((int)config.size() != N)
        throw ITError(tinyformat::format("make_product_state: config \"%s\" has %d characters, but there are %d sites",config,config.size(),N));
    if(basis != "z" && basis != "x" && basis != "y")
        throw ITError("make_product_state: basis must be \"z\", \"x\" or \"y\", got \"" + basis + "\"");

    // product state with link indices of dimension 1
    MPS psi(sites);

    for(int j : range1(N))
    {
        Index sj = sites(j);
        char c = config[j-1];

        if(!hasTags(sj,"Site,S=1/2"))
            throw ITError(tinyformat::format("make_product_state: site %d is not a spin-1/2",j));
        if(c != '0' && c != '1')
            throw ITError(tinyformat::format("make_product_state: invalid character '%c' in config (only '0' and '1' allowed)",c));

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

        set_site_tensor(&psi, sites, j, {a, b});
    }

    psi.position(1);
    return psi;
}


// ----------------------------------------------------------
// Given a mized spin-boson or spin-1/2 system it measures:
// - occupation number for bosons
// - magnetization along direction (x,y,z) for spin-1/2 
 
vector<double>
measure_magnetization(MPS* psi, const SiteSet sites , string direction)
{

    int N = length(sites);
    vector<double> mj;

    for(int j=1 ; j<=N ; j++)
    {
        Index sj = sites(j);
        Index sjp = prime(sites(j));

        ITensor S_j  = ITensor(sj ,sjp );


		if(hasTags(sj,"Site,Boson"))
		{
			for(int d=1; d <= dim(sj) ; d++) S_j.set(sj(d),sjp(d),d-1.);
        }
				
		if(hasTags(sj,"Site,S=1/2"))
		{
            if (direction == "x")
            {
                S_j.set(sj(1),sjp(2),1.);
			    S_j.set(sj(2),sjp(1),1.);	
                
            }
            else if(direction == "y")
            {
                S_j.set(sj(1),sjp(2), 1*Cplx_i);   // sigma^y: <down|..|up> = i, <up|..|down> = -i
			    S_j.set(sj(2),sjp(1),-1*Cplx_i);

            }
            else if(direction == "z")
            {
                S_j.set(sj(1),sjp(1),1.);
                S_j.set(sj(2),sjp(2),-1.);
            }
            else
            {
                cerr << "Direction choses is neither 'x' , 'y' or 'z'" << endl;
                return mj;
            }
		}
        (*psi).position(j);
        ITensor ket = (*psi)(j);
		ITensor bra = dag(prime((*psi)(j),"Site"));
		
		complex<double> exp_Sj = eltC(bra * S_j * ket);
		mj.push_back(real(exp_Sj));
        
    }

    return mj;
}


// Compute number of kinks (|\up_z \dw_z>) on a state psi

double 
measure_kink_number( MPS* psi, const SiteSet sites)
{
    int N = length(sites);

    double kink = 0.;

    if(N==1)
    {
        cerr << "Cannot measure number of kinks in a single-site system.\n";
        cerr << "Returnin 0.\n";
        return kink;
    }

    for(int j=1 ; j<N ; j++)
    {
        (*psi).position(j);
        
        ITensor N_1 = (op(sites,"Id",j)   - 2*op(sites,"Sz",j))    /2.;
        ITensor N_2 = (op(sites,"Id",j+1) + 2*op(sites,"Sz",j+1))  /2.;
        
        ITensor ket = (*psi)(j)*(*psi)(j+1);
		ITensor bra = dag(prime((*psi)(j),"Site"))*dag(prime((*psi)(j+1),"Site"));
		
        complex<double> n_j = eltC(bra * N_1 * N_2 * ket);
        kink += n_j.real();
    }

    return kink;
}


// ----------------------------------------------------------
// Compute correlation functions <N_start N_(start+i)> (both connected and disconnected).
// Where N = (1-2*S^z)/2 = |down_z> <down_z|

vector<double>
measure_density_correlations(MPS* psi, const SiteSet sites, const int start, const bool connected)
{
    int N = length(sites);
    vector<double> C;

    // reference site
    ITensor Ns = (op(sites,"Id",start) - 2*op(sites,"Sz",start))  /2.;

    vector<double> nj;
    if(connected)
    {
        vector<double> mz = measure_magnetization( psi,sites,"z");
        for(double m : mz) nj.push_back((1-m)/2.);
    }

    for(int j=1 ; j<= N; j++)
    {
        double Cjs;
        ITensor M;
        
        ITensor N1, N2;
        int jmax = max(j,start);
        int jmin = min(j,start);

        if(jmin == j)
        {
            N1 = (op(sites,"Id",j) - 2*op(sites,"Sz",j))  /2.;
            N2 = Ns;
        }
        else
        {
            N1 = Ns;
            N2 = (op(sites,"Id",j) - 2*op(sites,"Sz",j))  /2.;
        }

        (*psi).position(jmin);
        ITensor ket = (*psi)(jmin);


        if(j==start)
        {
		    ITensor bra = dag(prime((*psi)(j),"Site"));
            Cjs = eltC(ket*N1*bra).real(); //nb N^2 = N (it is a projector)
        }


        else
        {    
            Index ir = commonIndex( (*psi)(jmin) , (*psi)(jmin + 1) ,"Link");
			M = ket * N1 * dag( prime( prime( ket , "Site") , ir ) );

            for(int q = jmin + 1 ; q < jmax ; q++)
            {
                M *= (*psi)(q);
                M *= dag(prime( (*psi)(q) , "Link"));
            }

            Index il = commonIndex( (*psi)( jmax-1 ), (*psi)(jmax), "Link");
            M *= (*psi)(jmax);
            M *= N2;
            M *= dag( prime( prime((*psi)( jmax ), il) , "Site") );
            Cjs = eltC(M).real();
        }


        if(connected)
        {
            Cjs = Cjs - nj[jmin-1] * nj[jmax-1];
        }

        C.push_back(Cjs);

    }

    return C;
}


//----------------------------------------------------------------------

//measure of longitudinal and trasnversal magnetization in each site

void 
print_magnetization( const SpinHalf sites , MPS psi , const int N)
	{
	
	for( int j = 1 ; j <= N ; j++ )
		{
		psi.position(j);
		double Mx1 = 2 * eltC(dag(prime(psi(j),"Site")) * op(sites,"Sx",j) * psi(j)).real();
		double Mz1 = 2 * eltC(dag(prime(psi(j),"Site")) * op(sites,"Sz",j) * psi(j)).real();
		cout << "Sx_" << j << " = " << Mx1 << "\n"
			 << "Sz_" << j << " = " << Mz1 << endl;
		}
	}
