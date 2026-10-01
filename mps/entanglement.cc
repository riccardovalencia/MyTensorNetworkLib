/**
 * @file entanglement.cc
 * @brief Implementation of entanglement.h (the functions are documented in the header).
 */
#include "entanglement.h"
#include <itensor/all.h>
#include <cmath>

using namespace std;
using namespace itensor;


double
compute_entanglement_entropy( MPS* psi, int site, bool natural_log )
{
    (*psi).position(site);
    ITensor wf = (*psi)(site) * (*psi)(site+1);
    ITensor U  = (*psi)(site);
    ITensor S, V;
    auto spectrum = svd(wf, U, S, V);

    double entropy = 0.;
    for(auto p : spectrum.eigs())
    {
        if(p > 1E-12) entropy += -p * (natural_log ? log(p) : log2(p));
    }
    return entropy;
}
