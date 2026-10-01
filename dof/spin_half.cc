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
// Pauli matrix with input index `in` and output index `out` (spin-1/2: 1 = up_z, 2 = down_z)

ITensor
make_pauli_operator(const Index& in, const Index& out, const string& direction)
{
    ITensor S = ITensor(dag(in), out);
    if(direction == "x")
    {
        S.set(in(1),out(2),1.);
        S.set(in(2),out(1),1.);
    }
    else if(direction == "y")
    {
        S.set(in(1),out(2), 1*Cplx_i);   // sigma^y: <down|..|up> = i, <up|..|down> = -i
        S.set(in(2),out(1),-1*Cplx_i);
    }
    else if(direction == "z")
    {
        S.set(in(1),out(1), 1.);
        S.set(in(2),out(2),-1.);
    }
    else throw ITError("make_pauli_operator: direction must be \"x\", \"y\" or \"z\", got \"" + direction + "\"");
    return S;
}


// number operator on a boson site, Pauli matrix on a spin-1/2 site

ITensor
make_magnetization_operator(const Index& s, const string& direction)
{
    if(hasTags(s,"Site,S=1/2")) return make_pauli_operator(s, prime(s), direction);

    ITensor n = ITensor(dag(s), prime(s));
    if(hasTags(s,"Site,Boson")) for(int d = 1 ; d <= dim(s) ; d++) n.set(s(d),prime(s)(d),d-1.);
    return n;
}


// n_j = (1 - Z_j)/2 = |down_z><down_z|

static ITensor
make_density_operator(const SiteSet& sites, const int j)
{
    return (op(sites,"Id",j) - 2*op(sites,"Sz",j)) / 2.;
}


vector<double>
measure_magnetization(MPS* psi, const SiteSet sites , string direction)
{
    vector<double> mj;
    for(int j = 1 ; j <= length(sites) ; j++)
        mj.push_back(real(measure_local_operator(psi, make_magnetization_operator(sites(j), direction), j)));
    return mj;
}


// number of kinks |down_z up_z> on neighbouring sites

double 
measure_kink_number( MPS* psi, const SiteSet sites)
{
    double kink = 0.;
    for(int j = 1 ; j < length(sites) ; j++)
    {
        ITensor up = (op(sites,"Id",j+1) + 2*op(sites,"Sz",j+1)) / 2.;
        kink += real(measure_two_point_function(psi, sites, make_density_operator(sites, j), up, j, j+1));
    }
    return kink;
}


vector<double>
measure_density_correlations(MPS* psi, const SiteSet sites, const int start, const bool connected)
{
    ITensor Ns = make_density_operator(sites, start);

    vector<double> nj;
    if(connected)
        for(double m : measure_magnetization(psi, sites, "z")) nj.push_back((1-m)/2.);

    vector<double> C;
    for(int j = 1 ; j <= length(sites) ; j++)
    {
        // n is a projector: <n_s n_s> = <n_s>
        double Cjs = (j == start) ? real(measure_local_operator(psi, Ns, start))
                                  : real(measure_two_point_function(psi, sites, Ns, make_density_operator(sites, j), start, j));
        if(connected) Cjs -= nj[start-1] * nj[j-1];
        C.push_back(Cjs);
    }
    return C;
}


void 
print_magnetization( const SpinHalf sites , MPS psi , const int N)
{
    for( int j = 1 ; j <= N ; j++ )
    {
        double Mx = real(measure_local_operator(&psi, 2 * op(sites,"Sx",j), j));
        double Mz = real(measure_local_operator(&psi, 2 * op(sites,"Sz",j), j));
        cout << "Sx_" << j << " = " << Mx << "\n"
             << "Sz_" << j << " = " << Mz << endl;
    }
}
