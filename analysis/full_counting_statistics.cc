/**
 * @file full_counting_statistics.cc
 * @brief Implementation of full_counting_statistics.h (interfaces documented in the header, logic commented here).
 */
#include "full_counting_statistics.h"
#include <itensor/all.h>
#include <functional>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;
using namespace itensor;


// accumulate the steps from -pi: the coarse step in the outer quarters, the fine one in the middle
vector<double>
make_theta_grid( int number_points )
{
    vector<double> theta = {-M_PI};
    for(int k = 0 ; k < number_points - 1 ; k++)
    {
        bool outer = (k <= number_points / 4 || k >= (3 * number_points / 4 - 1));
        theta.push_back(theta.back() + (outer ? (M_PI - 1.) / (number_points / 4) : 2. / (number_points / 2.)));
    }
    return theta;
}


// N/2 minus half the block length (rounded up for even block sizes)
int
block_start( int N, int block_size )
{
    int size = block_size - 1;
    return (size % 2 == 0) ? N/2 - size/2 : N/2 - (size + 1)/2;
}


// contraction <psi| prod_{j in A} exp(i theta S^x_j) |psi> with the orthogonality center at the
// first site of the block: the parts of the chain outside the block contract to the identity

vector<complex<double> >
compute_generating_function( MPS* psi, const SpinHalf sites, int block_size, const vector<double>& theta )
{
    int N     = length(*psi);
    int start = block_start(N, block_size);
    int end   = start + block_size - 1;
    (*psi).position(start);

    vector<complex<double> > G;
    for(double th : theta)
    {
        function<ITensor(int)> phase = [&](int j) { return expHermitian(op(sites, "Sx", j), th * 1_i); };
        ITensor contraction;
        if(block_size == 1)
        {
            contraction = (*psi)(start) * phase(start) * dag(prime((*psi)(start), "Site"));
        }
        else
        {
            // first site: the left link is contracted between bra and ket
            Index right = commonIndex((*psi)(start), (*psi)(start + 1), "Link");
            contraction = (*psi)(start) * phase(start) * dag(prime(prime((*psi)(start), "Site"), right));
            for(int j = start + 1 ; j < end ; j++)
            {
                contraction *= (*psi)(j);
                contraction *= phase(j);
                contraction *= dag(prime((*psi)(j)));
            }
            // last site: the right link is contracted between bra and ket
            Index left = commonIndex((*psi)(end), (*psi)(end - 1), "Link");
            contraction *= (*psi)(end);
            contraction *= phase(end);
            contraction *= dag(prime(prime((*psi)(end), left), "Site"));
        }
        G.push_back(eltC(contraction));
    }
    return G;
}


// Tr(rho O): every site of the MPO is contracted with exp(i theta S^x_j) in the block and with the
// identity outside

vector<complex<double> >
compute_generating_function( MPO* rho, const SpinHalf sites, int block_size, const vector<double>& theta )
{
    int N     = length(*rho);
    int start = block_start(N, block_size);
    int end   = start + block_size - 1;

    vector<complex<double> > G;
    for(double th : theta)
    {
        ITensor trace = (*rho)(1);
        if(1 >= start && 1 <= end) trace *= expHermitian(op(sites, "Sx", 1), th * 1_i);
        else                       trace *= op(sites, "Id", 1);
        for(int j = 2 ; j <= N ; j++)
        {
            if(j >= start && j <= end)
            {
                trace *= (*rho)(j);
                trace *= expHermitian(op(sites, "Sx", j), th * 1_i);
            }
            else trace *= (*rho)(j) * op(sites, "Id", j);
        }
        G.push_back(eltC(trace));
    }
    return G;
}


// AutoMPO with S^x on every site of the block
MPO
make_block_sx_mpo( const SpinHalf sites, int start, int block_size )
{
    AutoMPO ampo(sites);
    for(int j = start ; j < start + block_size ; j++) ampo += 1., "Sx", j;
    return toMPO(ampo);
}


// moments <(S^x_A)^k>, k = 1..4, from powers of the MPO of S^x_A, then the standard
// moment-to-cumulant relations
vector<complex<double> >
compute_cumulants( MPS* psi, const SpinHalf sites, int block_size )
{
    int start = block_start(length(*psi), block_size);
    (*psi).position(start);
    MPO Sx = make_block_sx_mpo(sites, start, block_size);

    // powers of the MPO: in ITensor v3, A*B = nmultMPO(A, prime(B)) with primes 2 -> 1
    Args args_mult = {"MaxDim", 500, "Cutoff", 1E-16};
    MPO Sx2 = nmultMPO(Sx, prime(Sx), args_mult);   Sx2.mapPrime(2,1);
    MPO Sx3 = nmultMPO(Sx2, prime(Sx), args_mult);  Sx3.mapPrime(2,1);
    MPO Sx4 = nmultMPO(Sx2, prime(Sx2), args_mult); Sx4.mapPrime(2,1);

    // moments M_k = <(S^x_A)^k>
    complex<double> M1 = innerC(*psi, Sx, *psi);
    complex<double> M2 = innerC(*psi, Sx2, *psi);
    complex<double> M3 = innerC(*psi, Sx3, *psi);
    complex<double> M4 = innerC(*psi, Sx4, *psi);

    return {M1,
            M2 - M1*M1,
            M3 - 3.*M2*M1 + 2.*M1*M1*M1,
            M4 - 4.*M3*M1 - 3.*M2*M2 + 12.*M2*M1*M1 - 6.*M1*M1*M1*M1};
}
