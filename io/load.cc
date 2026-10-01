/**
 * @file load.cc
 * @brief Implementation of load.h (the functions are documented in the header).
 */
#include "load.h"
#include "../models/bosonic_east_model.h"
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


// "<results_dir><size>_cutoff<lambda>/mmGcbQEM_size..._c<c>/mmGcbQEM_size..._s<s>": folder of the
// ground states of one parameter set (a suffix "_v<version>" is added for the versioned layout)

static string
make_ground_state_dir(const string& results_dir, int size, int lambda, int symmetry_sector, double s, double c)
{
    string prefix = tinyformat::format("mmGcbQEM_size%d_cutoff%d_sector%d", size, lambda, symmetry_sector);
    return tinyformat::format("%s%d_cutoff%d/%s_c%.2f/%s_s%.2f", results_dir, size, lambda, prefix, c, prefix, s);
}


// State with the largest bond dimension bond_dimension * scaling^k found in dir, with its site set and
// energy variance; variance 100 and a placeholder state if dir has no site file.

static tuple<MPS, Boson, double>
load_largest_bond_dimension(const string& dir, int size, int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension, const string& symmetry_sector_dir)
{
    string sites_file = tinyformat::format("%s/sites_file_n0%d", dir, n0);
    if(!fileExists(sites_file))
    {
        Boson sites = Boson(1,{"ConserveQNs",false,"MaxOcc=",1});
        return {randomMPS(sites), sites, 100.};
    }

    Boson sites;
    readFromFile(sites_file, sites);

    auto state_file = [&](int chi) { return tinyformat::format("%s/ground_state_file_n0%d_chi%d", dir, n0, chi); };
    while(fileExists(state_file(bond_dimension))) bond_dimension = int(bond_dimension * scaling_bond_dimension);
    bond_dimension = int(bond_dimension / scaling_bond_dimension);

    MPS psi = randomMPS(sites);
    readFromFile(state_file(bond_dimension), psi);
    psi /= norm(psi);
    cerr << "Opened file : " << state_file(bond_dimension) << endl;

    double variance_H = compute_bosonic_east_model_energy_variance(&psi, sites, size , lambda, n0, symmetry_sector, s, c, symmetry_sector_dir);
    cerr << "variance : " << variance_H << endl;
    return {psi, sites, variance_H};
}


// versions _v1, _v2, ... are tried in order until one has variance below 1E-8

tuple<MPS, Boson, double>
load_ground_state_max_bond_dimension(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension, const string symmetry_sector_dir)
{
    const double tolerance_variance = 1E-8;
    string dir = make_ground_state_dir(results_dir, size, lambda, symmetry_sector, s, c);

    for(int version = 1 ; fileExists(tinyformat::format("%s_v%d/sites_file_n0%d", dir, version, n0)) ; version++)
    {
        auto loaded = load_largest_bond_dimension(tinyformat::format("%s_v%d", dir, version), size, lambda, n0, symmetry_sector, s, c, bond_dimension, scaling_bond_dimension, symmetry_sector_dir);
        if(get<2>(loaded) < tolerance_variance) return loaded;
    }
    Boson sites = Boson(1,{"ConserveQNs",false,"MaxOcc=",1});
    return {randomMPS(sites), sites, 100.};
}


tuple<MPS, Boson, double>
load_ground_state_max_bond_dimension_no_version(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension, const string symmetry_sector_dir)
{
    string dir = make_ground_state_dir(results_dir, size, lambda, symmetry_sector, s, c);
    return load_largest_bond_dimension(dir, size, lambda, n0, symmetry_sector, s, c, bond_dimension, scaling_bond_dimension, symmetry_sector_dir);
}


tuple<MPS, Boson>
load_adiabatic_state(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta )
{   
    string name_state;
    if( state_choice == 0 ) name_state = "super_coherent";
    if( state_choice == 3 ) name_state = "cat_state";

    string file_sites  = tinyformat::format("%s/sites_size%d_cutoff%d_s%.2f_c%.2f_alpha%.2f_adiabatic_state_choice%d_beta%.2f_linear",results_dir,size,cut_off, s,c,alpha.real(),state_choice,beta);
    string file_psi    = tinyformat::format("%s/psi_file_size%d_cutoff%d_s%.2f_c%.2f_alpha%.2f_adiabatic_state_choice%d_beta%.2f_linear",results_dir,size,cut_off, s,c,alpha.real(),state_choice,beta);

    Boson sites;

    if( fileExists( file_sites ) == true && fileExists( file_psi ) == true)
    {
        readFromFile(tinyformat::format("%s",file_sites), sites);        
        MPS psi = randomMPS(sites);
        readFromFile(tinyformat::format("%s",file_psi),psi);
        return  {psi , sites};
    }

    throw ITError("load_adiabatic_state: no file " + file_sites);
}
