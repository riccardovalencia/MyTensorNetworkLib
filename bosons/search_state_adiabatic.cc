#include <itensor/all.h>
#include <iostream>
#include <string>
#include <sstream>	//for ostringstream
#include <tuple>
#include <iomanip>

using namespace std;
using namespace itensor;


tuple<MPS, Boson>
search_state_adiabatic_coherent(string results_dir , const int size , const int cut_off, const double s, const double c, const complex<double> alpha, const int state_choice, double beta )
{   
    string name_state;
    if( state_choice == 0 ) name_state = "super_coherent";
    if( state_choice == 3 ) name_state = "cat_state";
    // string dir_state = format("%s/mmGcBQEM_adiabatic_%s_size%d_cutoff%d_sector0_c%.2f",results_dir,name_state, size, cut_off,c);
    // dir_state = format("%s/mmGcBQEM_adiabatic_%s_size%d_cutoff%d_sector0_s%.2f", dir_state,name_state, size,cut_off,s);

    string file_sites  = format("%s/sites_size%d_cutoff%d_s%.2f_c%.2f_alpha%.2f_adiabatic_state_choice%d_beta%.2f_linear",results_dir,size,cut_off, s,c,alpha.real(),state_choice,beta);
    string file_psi    = format("%s/psi_file_size%d_cutoff%d_s%.2f_c%.2f_alpha%.2f_adiabatic_state_choice%d_beta%.2f_linear",results_dir,size,cut_off, s,c,alpha.real(),state_choice,beta);

    Boson sites;

    if( fileExists( file_sites ) == true && fileExists( file_psi ) == true)
    {
        readFromFile(format("%s",file_sites), sites);        
        MPS psi = randomMPS(sites);
        readFromFile(format("%s",file_psi),psi);
        return  {psi , sites};
    }

    cerr << "No file found" << endl;
    cerr << file_sites << endl;
    exit(0);
}
