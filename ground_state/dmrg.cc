/**
 * @file dmrg.cc
 * @brief Implementation of dmrg.h (the functions are documented in the header).
 */
#include "dmrg.h"
#include "../dof/boson.h"
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
            writeToFile(tinyformat::format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);
            bond_dimension = int( bond_dimension * scaling_bond_dimension );
            }
        }

    else{
        *ground_state /= norm(*ground_state);
        writeToFile(tinyformat::format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);

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
            writeToFile(tinyformat::format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);
            bond_dimension = int( bond_dimension * scaling_bond_dimension );
            }
        }
   
    } while(  abs(delta_energy/sum_energy) > precision_dmrg && bond_dimension < max_bond_dimension);

*ground_state /= norm(*ground_state);
writeToFile(tinyformat::format("ground_state_file_n0%d_chi%d",n0,bond_dimension),*ground_state);

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


double variance = inner((*ground_state),H,H,(*ground_state))  -   inner((*ground_state),H,(*ground_state)) *  inner((*ground_state),H,(*ground_state));
double variance2 = inner((*ground_state),H,H,(*ground_state))  -  ground_state_energy * ground_state_energy;


save_file_ener << "# bond-dimension . energy_inner . variance_inner . energy_DMRG . variance_DMRG" << endl;
save_file_ener << 1 << " " << ground_state_energy << " " << variance << " " <<  inner((*ground_state),H,(*ground_state)) << " "<< variance2  << endl;

save_file_ener.close(); 
save_file_DMRG.close();


}


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
    int number_particles = (energy_target >= 0) ? int(energy_target) : 0;

	std::default_random_engine generator;
 	std::uniform_real_distribution<double> distribution(0, size-1);

	// start from the vacuum, keeping the link indices of *ground_state_variance
	initial_state_vacuum_state_correct_link( ground_state_variance, sites, size );

	// current_occupation[j] - 1 bosons on site j+1
	vector<int> current_occupation(size, 1);

	// add the bosons one at a time on random sites, without exceeding the cutoff
	for(int i=1; i <= number_particles ; i++)
	{
		int target_site;
		do{
			target_site	 = int(distribution(generator));
		}while(current_occupation[target_site]>cut_off_fock_space);

		put_occupation( ground_state_variance, sites, target_site + 1, current_occupation[target_site] );
		current_occupation[target_site] +=1 ;
	}

    (*ground_state_variance).position(1);
	(*ground_state_variance) /= norm((*ground_state_variance));

	cerr << "Measuring occupation number over the initial MPS state guessed. Energy target : " << energy_target << endl;
	vector<double> occupation_number;
	measure_occupation_number( ground_state_variance , sites , size ,  occupation_number );
	for(int j = 1 ; j <= size ; j++) cerr << j << " " << occupation_number[j-1] << endl;
}


MPS 
fermi_sea_electrons(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, Sweeps sweeps, double min_varH)
{
    int N = length(sites);

    if(Nupfill + Ndnfill > 2*N || Nupfill > N || Ndnfill > N)
    {
        cerr << "Filling larger than the one it can be hosted.\n";
        exit(-1);
    }

    InitState state = InitState(sites,"0");

    for(int j = 1; j <= Nupfill + Ndnfill; j += 1)
    {
        if(j <= Nupfill && j <= Ndnfill) state.set(j,"UpDn");
        else if(j<=Nupfill) state.set(j,"Up");
        else if(j<=Ndnfill) state.set(j,"Dn");
    }

    MPS psi0 = MPS(state); 

    cerr << "Computing ground state...\n";

    auto [energy,psi_gs] = dmrg(H,psi0,sweeps,{"Quiet",true});
    double E2 = inner(H,psi_gs,H,psi_gs) ;
    double E  = inner(psi_gs,H,psi_gs);
    double varE = E2 - E*E;
    
    if(abs(varE) > min_varH)
    {
        cerr << "Variance too large : " << varE << "\n";
        exit(-1);
    }

    // check filling
    AutoMPO ampo = AutoMPO(sites);
    for(int j : range1(N)) ampo += 1, "Nup" , j;
    MPO Nuptot = toMPO(ampo);

    ampo = AutoMPO(sites);
    for(int j : range1(N)) ampo += 1, "Ndn" , j;
    MPO Ndntot = toMPO(ampo);

    double nuptot = real(innerC(psi_gs,Nuptot,psi_gs));
    double ndntot = real(innerC(psi_gs,Ndntot,psi_gs));

     if( abs(Nupfill - nuptot) > 1E-7 ||  abs(Ndnfill - ndntot) > 1E-7 )
    {
        cerr << "Expected : " << Nupfill << " Computed : " << nuptot << "\n";
        cerr << "Expected : " << Ndnfill << " Computed : " << ndntot << "\n";
        exit(-1);
    }


    return psi_gs;

}
