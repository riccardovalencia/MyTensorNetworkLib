#include "perform_dmrg_variance.h"
#include "observables.h"
#include <itensor/all.h>

#include <sys/stat.h>
#include <iostream>
#include <string>
#include <fstream>	//output file
#include <sstream>	//for ostringstream
#include <iomanip>
#include <random>

using namespace std;
using namespace itensor;


double
perform_DMRG_variance(double energy_target, MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{

int size = physical_args.getInt("size");
int cut_off_fock_space = physical_args.getInt("cut_off_fock_space");
int n0 = physical_args.getInt("n0");
double s = physical_args.getReal("s");
double c = physical_args.getReal("c");

stringstream file_name_DMRG , file_name_ener;
file_name_DMRG << "excited_states_delta_energy_size" << size << "_s" << int(s*100) << "_c" << int(c*100) << "_cutoff" << cut_off_fock_space << "_n0" << n0 << ".dat";
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

ofstream save_file_DMRG( file_name_DMRG.str() , ios::app);
ofstream save_file_ener( file_name_ener.str() , ios::app);
save_file_DMRG << setprecision(set_output_precision);	
save_file_ener << setprecision(set_output_precision);	

int number_sweeps_done;
int number_sweeps_total = 1;
int number_sweeps_in_a_row = 0;

double current_noise = original_noise;

int array_maxdim [8] = {10,10,20,20,50,50,50,50};
int array_mindim [8] = {10,10,10,20,30,40,50,50};
double array_noise [8] = {current_noise, current_noise*1E-2, current_noise, current_noise*1E-2, current_noise, current_noise*1E-2, 0, 0};

save_file_DMRG << energy_target ;

for(int j=0; j<=7; j++)
	{
	auto sweeps = Sweeps( number_sweep_fixed_bond_dimension );
	sweeps.maxdim() = array_maxdim[j];	
	sweeps.mindim() = array_mindim[j];		
	sweeps.cutoff() = lower_bound_singular_values;
	sweeps.noise() = array_noise[j];	
	auto [energy,new_ground_state] = dmrg(H, *ground_state, sweeps);	
	*ground_state = new_ground_state;

	double delta_energy = energy - ground_state_energy ;
	double sum_energy = energy + ground_state_energy;
    ground_state_energy = energy;

    save_file_DMRG << " " << delta_energy/sum_energy;
        
	}

save_file_DMRG << endl;
// auto sweeps = Sweeps( number_sweep_fixed_bond_dimension );
// sweeps.maxdim() = 10,10,20,20,50,50;	
// sweeps.mindim() = 10,10,10,20,30,50;		
// sweeps.cutoff() = lower_bound_singular_values;
// sweeps.noise() = current_noise, current_noise*1E-2, current_noise, current_noise*1E-2, current_noise, 0;	
// auto [energy,new_ground_state] = dmrg(H, *ground_state, sweeps);	
// *ground_state = new_ground_state;

double variance = inner((*ground_state),H,H,(*ground_state))  -   inner((*ground_state),H,(*ground_state)) *  inner((*ground_state),H,(*ground_state));
double variance2 = inner((*ground_state),H,H,(*ground_state))  -  ground_state_energy * ground_state_energy;
 
save_file_ener << "# bond-dimension . energy_inner . variance_inner . energy_DMRG . variance_DMRG" << endl;
save_file_ener << int(bond_dimension/scaling_bond_dimension) << " " << ground_state_energy << " " << variance << " " <<  inner((*ground_state),H,(*ground_state)) << " "<< variance2  << endl;

save_file_ener.close(); 
save_file_DMRG.close();

bond_dimension = int( bond_dimension / scaling_bond_dimension );

return variance;

}



// initialization state before DMRG calculation over the variances



void
initialize_excited_state( MPS *ground_state_variance, const SiteSet sites, const int size, const double energy_target, const int cut_off_fock_space)
{
    int occupied_sites;
    int number_particles;
    if( energy_target >= 0 ) number_particles = int(energy_target);
    else number_particles = 0;
    // if( energy_target >= 0 ) occupied_sites = int(energy_target);
    // else occupied_sites = 0;

	std::default_random_engine generator;
 	std::uniform_real_distribution<double> distribution(0, size-1);


    // occupied_sites = int(number_particles / cut_off_fock_space);
    // if(occupied_sites>size) occupied_sites = size;

	// for(int i=1; i <= occupied_sites ; i++)

	vector<int> current_occupation;

	for(int i=1; i<=size; i++) current_occupation.push_back(1);

	for(int i=1; i<=size; i++)
	{
		auto si = sites(i);
		auto wf = ITensor(si);
		wf.set(si(1),1);
        for(int j=2; j<= cut_off_fock_space+1; j++) wf.set(si(j),0);
		(*ground_state_variance).set(i,wf);
	}

	for(int i=1; i <= number_particles ; i++)
	{
		
		int target_site;
		do{
			target_site	 = int(distribution(generator));
		}while(current_occupation[target_site]>cut_off_fock_space);


		auto si = sites(target_site + 1);
		auto wf = ITensor(si);

		// cerr << wf << endl;

		// wf.set(si(1),1);
        // wf.set(si(1),0);
        // wf.set(si(2),1);        
		// for(int j=3; j<= cut_off_fock_space+1 ; j++) wf.set(si(j),0);

		// old initialization
        // for(int j=1; j<= cut_off_fock_space; j++) wf.set(si(j),0);
		// wf.set(si(cut_off_fock_space+1),1);

		// new initialization

		// cerr << "current occupation target site : " << current_occupation[target_site] << endl;

		for(int j=1; j<= current_occupation[target_site]; j++) wf.set(si(j),0);
		wf.set(si(current_occupation[target_site]+1),1);
		for(int j=(current_occupation[target_site]+2); j<= cut_off_fock_space+1; j++) wf.set(si(j),0);
        
		cerr << "wf post initialization" << endl;
		cerr << wf << endl;

		current_occupation[target_site] +=1 ;
	

        // cerr << wf << endl;

		(*ground_state_variance).set(target_site + 1,wf);
	}

	// for(int i=occupied_sites+1; i <= size ; i++)
	// {
    //     cerr << "Changing site : " << i << endl;

	// 	auto si = sites(i);
	// 	auto wf = ITensor(si);

    //     wf.set(si(1),1);
    //     for(int j=2; j<= cut_off_fock_space+1; j++) wf.set(si(j),0);

	// 	// wf.set(si(cut_off_fock_space+1),0);

    //     // cerr << wf << endl;
	// 	(*ground_state_variance).set(i,wf);
	// }

    (*ground_state_variance).position(1);

	(*ground_state_variance) /= norm((*ground_state_variance));

	//test occupation number

	cerr << "Measuring occupation number over the initial MPS state guessed. Energy target : " << energy_target << endl;
	vector<double> occupation_number;

	measure_occupation_number( ground_state_variance , sites , size ,  occupation_number );


	for(int j = 1 ; j <= size ; j++)
		{
		cerr << j << " " << occupation_number[j-1] << endl;
		}
	
}