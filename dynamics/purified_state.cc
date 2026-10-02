/**
 * @file purified_state.cc
 * @brief Implementation of purified_state.h (interfaces documented in the header, logic commented here).
 */
#include "purified_state.h"
#include "../dof/spin_half.h"
#include <itensor/all.h>
#include <complex>
#include <functional>
#include <map>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// The 2N sites of the unfolded density matrix are
// | | | | |
// o-o-o-o-o-   (ket, sites N+1..2N)
// |
// o-o-o-o-o-   (bra, sites 1..N, reversed)
// | | | | |
// Tr(rho O) contracts bra site k with ket site 2N+1-k (physical site q = N+1-k), inserting O.

// Contract every pair (k, 2N+1-k), k = 1..N, through pair_operators[k] (indices of the two sites),
// or through a delta (trace) if the pair has no operator.
static Cplx
compute_folded_contraction(MPS* psi, const map<int, ITensor>& pair_operators)
{
    int N = length(*psi)/2;
    ITensor M;
    for(int k = 1 ; k <= N ; k++)
    {
        map<int, ITensor>::const_iterator it = pair_operators.find(k);
        ITensor O = (it != pair_operators.end()) ? it->second : delta(siteIndex(*psi,k), siteIndex(*psi,2*N-k+1));
        if(k == 1) M = (*psi)(k) * O * (*psi)(2*N-k+1);
        else
        {
            M *= (*psi)(k);
            M *= (*psi)(2*N-k+1);
            M *= O;
        }
    }
    return eltC(M);
}


// Operator O (indices s, s' of any site) moved onto the pair of physical site q:
// the input index s goes to the ket site, the output index s' to the bra site.
static ITensor
move_to_pair(const ITensor& O, MPS* psi, const int q)
{
    int N   = length(*psi)/2;
    Index ket = siteIndex(*psi, N+q);
    Index bra = siteIndex(*psi, N-q+1);
    Index col = noPrime(inds(O)[0]);
    Index row = prime(col);

    ITensor Oq = O;
    if(ket != col) Oq *= delta(ket, col);
    if(bra != row) Oq *= delta(bra, row);
    return Oq;
}


// Pauli matrix along direction on the pair of physical site q (zero on non-spin sites)
static ITensor
make_pauli_on_pair(const string& direction, MPS* psi, const int q)
{
    int N   = length(*psi)/2;
    Index ket = siteIndex(*psi, N+q);
    Index bra = siteIndex(*psi, N-q+1);
    if(!hasTags(ket,"Site,S=1/2")) return ITensor(ket, bra);
    return make_pauli_operator(ket, bra, direction);
}


// Tr(rho O_q) (divided by Tr(rho) if compute_normalization) on physical site q, or on all sites if q = -1
static vector<complex<double> >
measure_on_sites(MPS* psi, bool compute_normalization, int q, const function<ITensor(int)>& make_operator)
{
    int N = length(*psi)/2;
    double norm = compute_normalization ? compute_trace_purified(psi) : 1.;

    vector<int> physical_sites;
    if(q == -1) for(int p = 1 ; p <= N ; p++) physical_sites.push_back(p);
    else        physical_sites.push_back(q);

    vector<complex<double> > values;
    for(int p : physical_sites)
        values.push_back(compute_folded_contraction(psi, {{N-p+1, make_operator(p)}}) / norm);
    return values;
}


// folded contraction with a delta on every pair
double
compute_trace_purified(MPS* psi)
{
    return real(compute_folded_contraction(psi, {}));
}


// measure_on_sites with the Pauli matrix of each pair
vector<complex<double> >
measure_magnetization_purified(MPS *psi , string direction, bool compute_normalization, int q)
{
    return measure_on_sites(psi, compute_normalization, q, [&](int p) { return make_pauli_on_pair(direction, psi, p); });
}


// measure_on_sites with O moved onto each pair
vector<complex<double> >
measure_local_operator_purified(MPS *psi , const ITensor O, bool compute_normalization, int q)
{
    return measure_on_sites(psi, compute_normalization, q, [&](int p) { return move_to_pair(O, psi, p); });
}


// folded contraction with O on the pairs of q1 and q2 (and the two one-point functions if connected)
complex<double>
measure_correlation_purified(MPS *psi , const ITensor O, bool compute_normalization, int q1, int q2, bool connected)
{
    int N = length(*psi)/2;
    if(q1 == q2) throw ITError("measure_correlation_purified: q1 == q2 (autocorrelation) is not implemented");
    if(q1 < 1 || q2 < 1 || q1 > N || q2 > N) throw ITError("measure_correlation_purified: sites out of the chain");

    double norm = compute_normalization ? compute_trace_purified(psi) : 1.;
    complex<double> O1O2 = compute_folded_contraction(psi, {{N-q1+1, move_to_pair(O, psi, q1)},
                                                            {N-q2+1, move_to_pair(O, psi, q2)}}) / norm;
    if(connected)
    {
        complex<double> O1 = measure_local_operator_purified(psi, O, compute_normalization, q1)[0];
        complex<double> O2 = measure_local_operator_purified(psi, O, compute_normalization, q2)[0];
        O1O2 -= O1 * O2;
    }
    return O1O2;
}
