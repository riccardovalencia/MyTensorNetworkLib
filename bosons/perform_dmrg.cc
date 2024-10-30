#include "perform_dmrg.h"
#include "observables.h"
#include <itensor/all.h>

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>

using namespace std;
using namespace itensor;


int
perform_DMRG(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{

int size = physical_args.getInt("size");
int cut_off_fock_space = physical_args.getInt("cut_off_fock_space");
int n0 = physical_args.getInt("n0");
double s = physical_args.getReal("s");
double c = physical_args.getReal("c");

stringstream file_name_DMRG , file_name_ener;
file_name_DMRG << "delta_energy_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
file_name_ener << "energy_size" << size << "_s" << int(s*100) <<  "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";


int bond_dimension = numerical_args.getInt("bond_dimension");
int scaling_bond_dimension = numerical_args.getInt("scaling_bond_dimension");
int max_bond_dimension = numerical_args.getInt("max_bond_dimension");
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

double current_noise = original_noise;

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
		sweeps.noise() = current_noise;	
        auto [energy,new_ground_state] = dmrg(H, *ground_state, sweeps);	
        *ground_state = new_ground_state;
        delta_energy = energy - ground_state_energy ;
		sum_energy = energy + ground_state_energy;
        ground_state_energy = energy;

        save_file_DMRG << number_sweeps_total*number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy <<endl;
        
        number_sweeps_done  += 1;
        number_sweeps_total += 1;
        
        } while ( abs(delta_energy/sum_energy) > precision_dmrg && number_sweeps_done < 2 );

    if(number_sweeps_done >=2 ) cerr << "Increased bond dimension because number weeps exceeded 5!!!" << endl;

    if(number_sweeps_done >=2 )
        {
        if( current_noise > current_noise*1E-01) current_noise = 0;
        else
            {
            cerr << "Noise is zero and it doesn't converge in 5 number sweeps! Increasing bond dimension and putting non zero noise" << endl;
            current_noise = original_noise;
            *ground_state /= norm(*ground_state);
            writeToFile(format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);
            bond_dimension = int( bond_dimension * scaling_bond_dimension );
            }
        }

    else{
        *ground_state /= norm(*ground_state);
        writeToFile(format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);

        bond_dimension = int( bond_dimension * scaling_bond_dimension );
        if(number_sweeps_done == 1)
            {
            number_sweeps_in_a_row += 1;
            current_noise *= 1E-02  ;
            }
        else if(number_sweeps_in_a_row >= 1)
            {
            current_noise = 0;
            }
        else 
            {
            number_sweeps_in_a_row = 0;
            current_noise = original_noise;
            }
    }
   
    } while( number_sweeps_in_a_row < 2 && bond_dimension < max_bond_dimension);

if( bond_dimension > max_bond_dimension)
    {
    cerr << "Reached max bond dimension available. Aborted" << endl;
	exit(0);
    }

double variance = inner((*ground_state),H,H,(*ground_state))  -   inner((*ground_state),H,(*ground_state)) *  inner((*ground_state),H,(*ground_state));
double variance2 = inner((*ground_state),H,H,(*ground_state))  -  ground_state_energy * ground_state_energy;
 
save_file_ener << "# bond-dimension . energy_inner . variance_inner . energy_DMRG . variance_DMRG" << endl;
save_file_ener << int(bond_dimension/scaling_bond_dimension) << " " << ground_state_energy << " " << variance << " " <<  inner((*ground_state),H,(*ground_state)) << " "<< variance2  << endl;

save_file_ener.close(); 
save_file_DMRG.close();

bond_dimension = int( bond_dimension / scaling_bond_dimension );

return bond_dimension;

}



int
perform_DMRG_soft(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{

int size = physical_args.getInt("size");
int cut_off_fock_space = physical_args.getInt("cut_off_fock_space");
int n0 = physical_args.getInt("n0");
double s = physical_args.getReal("s");
double c = physical_args.getReal("c");

stringstream file_name_DMRG , file_name_ener;
file_name_DMRG << "delta_energy_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
file_name_ener << "energy_size" << size << "_s" << int(s*100) <<  "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";


int bond_dimension = numerical_args.getInt("bond_dimension");
int scaling_bond_dimension = numerical_args.getInt("scaling_bond_dimension");
int max_bond_dimension = numerical_args.getInt("max_bond_dimension");
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

double current_noise = original_noise;

double delta_energy;
double sum_energy;

do
    {
    delta_energy = 0;
    sum_energy  = 0;
    number_sweeps_done = 0;
    do
        {
        auto sweeps = Sweeps( number_sweep_fixed_bond_dimension );
        sweeps.maxdim() = bond_dimension;			
        sweeps.cutoff() = lower_bound_singular_values;
		sweeps.noise() = current_noise;	
        auto [energy,new_ground_state] = dmrg(H, *ground_state, sweeps);	
        *ground_state = new_ground_state;
        delta_energy = energy - ground_state_energy ;
		sum_energy = energy + ground_state_energy;
        ground_state_energy = energy;

        save_file_DMRG << number_sweeps_total*number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy <<endl;
        
        number_sweeps_done  += 1;
        
        } while ( abs(delta_energy/sum_energy) > precision_dmrg && number_sweeps_done < 3 );

    if(number_sweeps_done >=3 &&  abs(delta_energy/sum_energy) > precision_dmrg)
        {
        if( current_noise > current_noise*1E-01) current_noise = 0;
        else
            {
            cerr << "Noise is zero and it doesn't converge in 5 number sweeps! Increasing bond dimension and putting non zero noise" << endl;
            current_noise = original_noise;
            *ground_state /= norm(*ground_state);
            writeToFile(format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);
            bond_dimension = int( bond_dimension * scaling_bond_dimension );
            }
        }
   
    } while(  abs(delta_energy/sum_energy) > precision_dmrg && bond_dimension < max_bond_dimension);

*ground_state /= norm(*ground_state);
writeToFile(format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);

if( bond_dimension > max_bond_dimension)
    {
    cerr << "Reached max bond dimension available. Aborted" << endl;
	exit(0);
    }

double variance = inner((*ground_state),H,H,(*ground_state))  -   inner((*ground_state),H,(*ground_state)) *  inner((*ground_state),H,(*ground_state));
double variance2 = inner((*ground_state),H,H,(*ground_state))  -  ground_state_energy * ground_state_energy;
 
save_file_ener << "# bond-dimension . energy_inner . variance_inner . energy_DMRG . variance_DMRG" << endl;
save_file_ener << bond_dimension << " " << ground_state_energy << " " << variance << " " <<  inner((*ground_state),H,(*ground_state)) << " "<< variance2  << endl;

save_file_ener.close(); 
save_file_DMRG.close();


return bond_dimension;

}
