/**
 * @file load.cc
 * @brief Implementation of load.h (the functions are documented in the header).
 */
#include "load.h"
#include "../models/bqem.h"
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


tuple<MPS, Boson, double>
search_ground_state_max_bond_chi(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension, const string symmetry_sector_dir)
{
    double tolerance_variance = 1E-8;
    stringstream  name_dir_cutoff;
    
    name_dir_cutoff << results_dir << size << "_cutoff" << lambda ; 


    stringstream name_dir_fixed_c;

    name_dir_fixed_c << name_dir_cutoff.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_c" ;
    name_dir_fixed_c << fixed << setprecision(2) << c ;


    stringstream name_dir_s_prefix ;

    name_dir_s_prefix << name_dir_fixed_c.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_s" ;
    name_dir_s_prefix << fixed << setprecision(2) << s << "_v" ;

    int version = 1 ; 
    stringstream name_dir_s ; 
    name_dir_s << name_dir_s_prefix.str() << version;

    Boson sites;

    while( fileExists( tinyformat::format("%s/sites_file_n0%d",name_dir_s.str(), n0) ) == true ){
        readFromFile(tinyformat::format("%s/sites_file_n0%d",name_dir_s.str(),n0), sites);        
        MPS psi = randomMPS(sites);

        while(fileExists( tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) ) == true)
        {
        bond_dimension = int(bond_dimension * scaling_bond_dimension);
        }
        bond_dimension = int(bond_dimension / scaling_bond_dimension ); 

        readFromFile(tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(), n0 ,bond_dimension),psi);
        psi /= norm(psi);

        cerr << "Opened file : " << tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) << endl;

        double variance_H = compute_variance_H_mmGcbQEM(&psi, sites, size , lambda, n0, symmetry_sector, s, c, symmetry_sector_dir);

        cerr << "variance : " << variance_H << endl;

        if(variance_H < tolerance_variance) return {psi , sites, variance_H};

        version += 1;
        name_dir_s_prefix.str("");
        name_dir_s << name_dir_s_prefix.str() << version;
    }
    double variance_H = 100;
    sites = Boson(1,{"ConserveQNs",false,"MaxOcc=",1});	
    MPS psi = randomMPS(sites);
    return {psi , sites, variance_H};
}


tuple<MPS, Boson, double>
search_ground_state_max_bond_chi_no_v(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension, const string symmetry_sector_dir)
{
    double tolerance_variance = 100000;
    stringstream  name_dir_cutoff;
    
    name_dir_cutoff << results_dir << size << "_cutoff" << lambda ; 


    stringstream name_dir_fixed_c;

    name_dir_fixed_c << name_dir_cutoff.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_c" ;
    name_dir_fixed_c << fixed << setprecision(2) << c ;


    stringstream name_dir_s ;

    name_dir_s << name_dir_fixed_c.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_s" ;
    name_dir_s << fixed << setprecision(2) << s ;


    Boson sites;

    if( fileExists( tinyformat::format("%s/sites_file_n0%d",name_dir_s.str(), n0) ) == true )
    {
        readFromFile(tinyformat::format("%s/sites_file_n0%d",name_dir_s.str(),n0), sites);        
        MPS psi = randomMPS(sites);

        while(fileExists( tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) ) == true)
        {
        bond_dimension = int(bond_dimension * scaling_bond_dimension);
        }
        bond_dimension = int(bond_dimension / scaling_bond_dimension ); 

        readFromFile(tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(), n0 ,bond_dimension),psi);
        psi /= norm(psi);
        cerr << "Opened file : " << tinyformat::format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) << endl;
        double variance_H = compute_variance_H_mmGcbQEM(&psi, sites, size , lambda, n0, symmetry_sector, s, c, symmetry_sector_dir);
        cerr << "variance : " << variance_H << endl;

        if(variance_H < tolerance_variance) return {psi , sites, variance_H};

    }
    double variance_H = 100;
    sites = Boson(1,{"ConserveQNs",false,"MaxOcc=",1});	
    MPS psi = randomMPS(sites);
    return {psi , sites, variance_H};
}


tuple<MPS, Boson>
search_state_adiabatic_coherent(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta )
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

    cerr << "No file found" << endl;
    cerr << file_sites << endl;
    exit(0);
}
