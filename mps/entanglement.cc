/**
 * @file entanglement.cc
 * @brief Implementation of entanglement.h (interfaces documented in the header, logic commented here).
 */
#include "entanglement.h"
#include <itensor/all.h>
#include <cmath>
#include <vector>

using namespace std;
using namespace itensor;


// Orthogonality center on site, SVD of the two-site tensor of the bond (site, site+1): the
// eigenvalues of the spectrum are the squared Schmidt values p_k.
double
compute_entanglement_entropy( MPS* psi, int site, bool natural_log )
{
    (*psi).position(site);
    ITensor wf = (*psi)(site) * (*psi)(site+1);
    ITensor U  = (*psi)(site);
    ITensor S, V;
    Spectrum spectrum = svd(wf, U, S, V);

    double entropy = 0.;
    for(double p : spectrum.eigs())
    {
        if(p > 1E-12) entropy += -p * (natural_log ? log(p) : log2(p));
    }
    return entropy;
}


// one SVD per bond, moving the orthogonality center along the chain
vector<double>
compute_entanglement_entropies( MPS* psi, bool natural_log )
{
    vector<double> entropies;
    for(int j = 1 ; j < length(*psi) ; j++) entropies.push_back(compute_entanglement_entropy(psi, j, natural_log));
    return entropies;
}
