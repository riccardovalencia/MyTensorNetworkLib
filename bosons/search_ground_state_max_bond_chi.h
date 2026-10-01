#ifndef MYTN_BOSONS_SEARCH_GROUND_STATE_MAX_BOND_CHI_H
#define MYTN_BOSONS_SEARCH_GROUND_STATE_MAX_BOND_CHI_H
//get_data.h

#include <itensor/all.h>
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

// Look in results_dir (folder structure of perform_DMRG) for the ground state with the largest bond
// dimension available, starting from bond_dimension and scaling it by scaling_bond_dimension; stops at the
// first state whose energy variance is below a fixed tolerance. Returns {state, sites, variance}.
tuple<MPS, Boson, double>
search_ground_state_max_bond_chi(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension);


// As above, for the folder layout without the symmetry-sector subfolder.
tuple<MPS, Boson, double>
search_ground_state_max_bond_chi_no_v(string results_dir, int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension , double scaling_bond_dimension);

#endif
