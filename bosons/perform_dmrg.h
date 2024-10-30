#ifndef PERFORM_DMRG_H
#define PERFORM_DMRG_H

#include <itensor/all.h>
#include <iostream>
#include <fstream>	//output file
#include <sstream> // for ostringstream
#include <string>
#include <iomanip>

using namespace std;
using namespace itensor;

int
perform_DMRG(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

int
perform_DMRG_soft(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args);

#endif
