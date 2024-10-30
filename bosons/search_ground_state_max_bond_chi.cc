#include <itensor/all.h>
#include <iostream>
#include <string>
#include <sstream>	//for ostringstream
#include <tuple>
#include <iomanip>
#include "observables.h"

using namespace std;
using namespace itensor;


tuple<MPS, Boson, double>
search_ground_state_max_bond_chi(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension)
{
    double tolerance_variance = 1E-8;
    stringstream  name_dir_cutoff;
    
    // name_dir_cutoff << "/home/ricval/Documenti/Bosonic/GcbQEM_Data/Results_Noise_size" << size << "_cutoff" << lambda ; 
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

    while( fileExists( format("%s/sites_file_n0%d",name_dir_s.str(), n0) ) == true ){
        readFromFile(format("%s/sites_file_n0%d",name_dir_s.str(),n0), sites);        
        MPS psi = randomMPS(sites);

        while(fileExists( format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) ) == true)
        {
        bond_dimension = int(bond_dimension * scaling_bond_dimension);
        }
        bond_dimension = int(bond_dimension / scaling_bond_dimension ); 

        readFromFile(format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(), n0 ,bond_dimension),psi);
        psi /= norm(psi);

        cerr << "Opened file : " << format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) << endl;

        double variance_H = compute_variance_H_mmGcbQEM(&psi, sites, size , lambda, n0, symmetry_sector, s, c);

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
search_ground_state_max_bond_chi_no_v(string results_dir , int size , int lambda, int n0, int symmetry_sector, double s, double c, int bond_dimension, double scaling_bond_dimension)
{
    double tolerance_variance = 100000;
    stringstream  name_dir_cutoff;
    
    // name_dir_cutoff << "/home/ricval/Documenti/Bosonic/GcbQEM_Data/Results_Noise_size" << size << "_cutoff" << lambda ; 
    name_dir_cutoff << results_dir << size << "_cutoff" << lambda ; 


    stringstream name_dir_fixed_c;

    name_dir_fixed_c << name_dir_cutoff.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_c" ;
    name_dir_fixed_c << fixed << setprecision(2) << c ;


    stringstream name_dir_s ;

    name_dir_s << name_dir_fixed_c.str() << "/mmGcbQEM_size" << size << "_cutoff" << lambda << "_sector" << symmetry_sector << "_s" ;
    name_dir_s << fixed << setprecision(2) << s ;

    
    Boson sites;

    if( fileExists( format("%s/sites_file_n0%d",name_dir_s.str(), n0) ) == true )
    {
        readFromFile(format("%s/sites_file_n0%d",name_dir_s.str(),n0), sites);        
        MPS psi = randomMPS(sites);

        while(fileExists( format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) ) == true)
        {
        bond_dimension = int(bond_dimension * scaling_bond_dimension);
        }
        bond_dimension = int(bond_dimension / scaling_bond_dimension ); 

        readFromFile(format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(), n0 ,bond_dimension),psi);
        psi /= norm(psi);
        cerr << "Opened file : " << format("%s/ground_state_file_n0%d_chi%d",name_dir_s.str(),n0,bond_dimension) << endl;
        double variance_H = compute_variance_H_mmGcbQEM(&psi, sites, size , lambda, n0, symmetry_sector, s, c);
        cerr << "variance : " << variance_H << endl;

        if(variance_H < tolerance_variance) return {psi , sites, variance_H};

    }
    double variance_H = 100;
    sites = Boson(1,{"ConserveQNs",false,"MaxOcc=",1});	
    MPS psi = randomMPS(sites);
    return {psi , sites, variance_H};
}

  
