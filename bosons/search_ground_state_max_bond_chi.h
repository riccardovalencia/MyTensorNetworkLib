#ifndef SEARCH_GROUND_STATE_MAX_BOND_CHI
#define SEARCH_GROUND_STATE_MAX_BOND_CHI
//get_data.h

#include <iostream>
#include <string>
#include <tuple>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>

using namespace std;
using namespace itensor;

//----------------------------------------------------------------------
// look for the ground state computed at the end of DMRG calculation

tuple<MPS, Boson, double>
search_ground_state_max_bond_chi(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension);


tuple<MPS, Boson, double>
search_ground_state_max_bond_chi_no_v(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension);

#endif
