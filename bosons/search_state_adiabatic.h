#include <itensor/all.h>
#include <iostream>
#include <string>
#include <sstream>	//for ostringstream
#include <tuple>
#include <iomanip>
#include "observables.h"

using namespace std;
using namespace itensor;


tuple<MPS, Boson>
search_state_adiabatic_coherent(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta );
