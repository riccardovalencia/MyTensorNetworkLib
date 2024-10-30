#include "perform_dmrg_meanfield.h"
#include <itensor/all.h>

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>

using namespace std;
using namespace itensor;


void
perform_DMRG_meanfield(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{

int size = physical_args.getInt("size");
int cut_off_fock_space = physical_args.getInt("cut_off_fock_space");
int n0 = physical_args.getInt("n0");
double s = physical_args.getReal("s");

stringstream file_name_DMRG , file_name_ener;
file_name_DMRG << "meanfield_delta_energy_size" << size << "_s" << int(s*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
file_name_ener << "meanfield_energy_size" << size << "_s" << int(s*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";



int bond_dimension = 1;
int number_sweep_fixed_bond_dimension = numerical_args.getInt("number_sweep_fixed_bond_dimension");

double precision_dmrg = numerical_args.getReal("precision_dmrg");
double lower_bound_singular_values = numerical_args.getReal("lower_bound_singular_values");
double original_noise = numerical_args.getReal("original_noise");

(*ground_state).position(1);
(*ground_state) /= norm((*ground_state));
double ground_state_energy = inner( (*ground_state) , H , (*ground_state) ); //starting energy

ofstream save_file_DMRG( file_name_DMRG.str() );
ofstream save_file_ener( file_name_ener.str() );
save_file_DMRG << setprecision(set_output_precision);	
save_file_ener << setprecision(set_output_precision);	

int number_sweeps_done;
int number_sweeps_total = 1;
int number_sweeps_in_a_row = 0;

double current_noise = 0.;

do
    {
    double delta_energy;
	double sum_energy;
    number_sweeps_done = 0;
    do
        {
        auto sweeps = Sweeps( number_sweep_fixed_bond_dimension );
        sweeps.maxdim() = bond_dimension;			
        sweeps.cutoff() = lower_bound_singular_values;
		//sweeps.noise() = current_noise;	
        auto [energy,new_ground_state] = dmrg(H, *ground_state, sweeps);	
        *ground_state = new_ground_state;
        delta_energy = energy - ground_state_energy ;
		sum_energy = energy + ground_state_energy;
		cerr << "Delta energy : " << delta_energy << endl;
		cerr << "Sum energy : " << sum_energy << endl;
		cerr << "Fraction : " << delta_energy / sum_energy << endl;
        ground_state_energy = energy;
        save_file_DMRG << number_sweeps_total*number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy <<endl;
        number_sweeps_done += 1;
        number_sweeps_total += 1;
        } while ( abs(delta_energy/sum_energy) > precision_dmrg && number_sweeps_done < 5 );

    if(number_sweeps_done == 1)
		{
		number_sweeps_in_a_row += 1;
		}
    else 
		{
		number_sweeps_in_a_row = 0;
		}
    } while( number_sweeps_in_a_row < 3 );


// MPO H2 = MPO(sites);
// H2 = nmultMPO(prime(H),H,{"MaxDim",2000,"Cutoff",1E-15});

// save_file_ener << H << endl;
// save_file_ener << H2 << endl;

// attempt - I build H2

// MPO H2bis = MPO(sites);


// for( int j=1; j<= size ; j++)
//     {   
//     H2bis.ref(j) = prime(H(j)) * H(j);
//     }

// fai (H - Lambda Id)^2

// double variance3 = inner((*ground_state),H2,(*ground_state));

// double variance2 = inner((*ground_state),H2bis,(*ground_state))  ;

double variance = inner((*ground_state),H,H,(*ground_state))  -   inner((*ground_state),H,(*ground_state)) *  inner((*ground_state),H,(*ground_state));
double variance2 = inner((*ground_state),H,H,(*ground_state))  -  ground_state_energy * ground_state_energy;
 
// cerr << setprecision(set_output_precision);	 
// cerr  << variance << " " << variance2 << " " << variance3 << endl;
//  - ground_state_energy*ground_state_energy;


save_file_ener << "# bond-dimension . energy_inner . variance_inner . energy_DMRG . variance_DMRG" << endl;
save_file_ener << 1 << " " << ground_state_energy << " " << variance << " " <<  inner((*ground_state),H,(*ground_state)) << " "<< variance2  << endl;

save_file_ener.close(); 
save_file_DMRG.close();


}


