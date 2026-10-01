#ifndef MYTN_BOSONS_SEARCH_STATE_ADIABATIC_H
#define MYTN_BOSONS_SEARCH_STATE_ADIABATIC_H

#include <itensor/all.h>
#include <iostream>
#include <string>
#include <sstream>	//for ostringstream
#include <tuple>
#include <iomanip>
#include "observables.h"

using namespace std;
using namespace itensor;


// Load from results_dir the state prepared by adiabatic dressing of a coherent state of amplitude alpha
// (files "sites_size..." and "psi_file_size..." written by the adiabatic drivers). Returns {state, sites}.
tuple<MPS, Boson>
search_state_adiabatic_coherent(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta );

#endif
