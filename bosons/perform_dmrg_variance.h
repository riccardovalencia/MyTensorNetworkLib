#ifndef PERFORM_DMRG_VARIANCE_H
#define PERFORM_DMRG_VARIANCE_H

#include <itensor/all.h>
#include <iostream>
#include <fstream>	//output file
#include <sstream> // for ostringstream
#include <string>
#include <iomanip>

using namespace std;
using namespace itensor;

double
perform_DMRG_variance(double energy_target, MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

void
initialize_excited_state( MPS *ground_state_variance, const SiteSet sites, const int size, const double energy_target, const int cut_off_fock_space);


#endif
