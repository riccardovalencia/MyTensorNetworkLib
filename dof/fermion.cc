/**
 * @file fermion.cc
 * @brief Implementation of fermion.h (the functions are documented in the header).
 */
#include "fermion.h"
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


MPS
initial_computational_electron_state(const SiteSet sites, const int Nupfill, const int Ndnfill)
{
    int N = length(sites);

    if(Nupfill + Ndnfill > 2*N || Nupfill > N || Ndnfill > N)
    {
        cerr << "Filling larger than the one it can be hosted.\n";
        exit(-1);
    }

    InitState state = InitState(sites,"0");

    for(int j = 1; j <= Nupfill + Ndnfill; j += 1)
    {
        if(j <= Nupfill && j <= Ndnfill) state.set(j,"UpDn");
        else if(j<=Nupfill) state.set(j,"Up");
        else if(j<=Ndnfill) state.set(j,"Dn");
    }

    return MPS(state); 

}
