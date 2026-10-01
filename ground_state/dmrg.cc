/**
 * @file dmrg.cc
 * @brief Implementation of dmrg.h (the functions are documented in the header).
 */
#include "dmrg.h"
#include "../dof/boson.h"
#include "../dof/fermion.h"
#include <itensor/all.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace itensor;


// ----------------------------------------------------------
// shared parts of the DMRG drivers

// numerical parameters (see dmrg.h)
struct DmrgParameters
{
    int bond_dimension;
    int scaling_bond_dimension;
    int max_bond_dimension;
    int number_sweep_fixed_bond_dimension;
    double precision_dmrg;
    double lower_bound_singular_values;
    double original_noise;
};

static DmrgParameters
read_dmrg_parameters(Args const& numerical_args)
{
    DmrgParameters p;
    p.bond_dimension                    = numerical_args.getInt("bond_dimension", 1);
    p.scaling_bond_dimension            = numerical_args.getInt("scaling_bond_dimension", 1);
    p.max_bond_dimension                = numerical_args.getInt("max_bond_dimension", 1);
    p.number_sweep_fixed_bond_dimension = numerical_args.getInt("number_sweep_fixed_bond_dimension");
    p.precision_dmrg                    = numerical_args.getReal("precision_dmrg");
    p.lower_bound_singular_values       = numerical_args.getReal("lower_bound_singular_values");
    p.original_noise                    = numerical_args.getReal("original_noise");
    return p;
}


// "<prefix>_size<size>_s<100 s>[_c<100 c>]_cutoff<cutoff>_n0<n0>.dat"
static string
make_output_name(const string& prefix, Args const& physical_args, const bool include_c)
{
    stringstream name;
    name << prefix << "_size" << physical_args.getInt("size") << "_s" << int(physical_args.getReal("s")*100);
    if(include_c) name << "_c" << int(physical_args.getReal("c")*100);
    name << "_cutoff" << physical_args.getInt("cut_off_fock_space") << "_n0" << physical_args.getInt("n0") << ".dat";
    return name.str();
}


// one call of ITensor dmrg with fixed maximal bond dimension and noise; returns the energy
static double
run_dmrg_sweeps(MPS* psi, const MPO& H, const int number_sweeps, const int max_dim, const double cutoff, const double noise, const int min_dim = 1)
{
    auto sweeps = Sweeps(number_sweeps);
    sweeps.maxdim() = max_dim;
    sweeps.mindim() = min_dim;
    sweeps.cutoff() = cutoff;
    sweeps.noise()  = noise;
    auto [energy, new_state] = dmrg(H, *psi, sweeps);
    *psi = new_state;
    return energy;
}


// normalize psi and write it to "ground_state_file_n0<n0>_chi<bond_dimension>"
static void
save_ground_state(MPS* psi, const int n0, const int bond_dimension)
{
    *psi /= norm(*psi);
    writeToFile(tinyformat::format("ground_state_file_n0%d_chi%d",n0,bond_dimension), *psi);
}


// write the final energies and variances; returns the variance <H^2> - <H>^2
static double
write_energy_and_variance(ofstream& file, MPS* psi, const MPO& H, const int bond_dimension, const double energy_dmrg)
{
    double H2       = inner(*psi,H,H,*psi);
    double energy   = inner(*psi,H,*psi);
    double variance = H2 - energy * energy;
    file << "# bond-dimension . energy_DMRG . variance . energy . variance_DMRG" << endl;
    file << bond_dimension << " " << energy_dmrg << " " << variance << " " << energy << " " << H2 - energy_dmrg * energy_dmrg << endl;
    return variance;
}


// normalized psi, orthogonality center on site 1; returns <psi|H|psi>
static double
prepare_initial_state(MPS* psi, const MPO& H)
{
    (*psi).position(1);
    (*psi) /= norm(*psi);
    return inner(*psi, H, *psi);
}


// ----------------------------------------------------------

int
perform_dmrg(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{
    int n0 = physical_args.getInt("n0");
    DmrgParameters p = read_dmrg_parameters(numerical_args);
    int bond_dimension = p.bond_dimension;

    ofstream save_file_DMRG( make_output_name("delta_energy", physical_args, true) );
    ofstream save_file_ener( make_output_name("energy", physical_args, true) );
    save_file_DMRG << setprecision(set_output_precision);
    save_file_ener << setprecision(set_output_precision);

    double ground_state_energy = prepare_initial_state(ground_state, H);

    int number_sweeps_total    = 1;
    int number_sweeps_in_a_row = 0;
    double current_noise       = p.original_noise;

    do
    {
        // at most 2 rounds of sweeps at fixed bond dimension
        double delta_energy, sum_energy;
        int number_sweeps_done = 0;
        do
        {
            double energy = run_dmrg_sweeps(ground_state, H, p.number_sweep_fixed_bond_dimension, bond_dimension, p.lower_bound_singular_values, current_noise);
            delta_energy = energy - ground_state_energy;
            sum_energy   = energy + ground_state_energy;
            ground_state_energy = energy;

            save_file_DMRG << number_sweeps_total*p.number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy << endl;
            number_sweeps_done  += 1;
            number_sweeps_total += 1;
        } while ( abs(delta_energy/sum_energy) > p.precision_dmrg && number_sweeps_done < 2 );

        if(number_sweeps_done >= 2)
        {
            // not converged: first switch the noise off, then increase the bond dimension
            cerr << "Not converged in 2 rounds of sweeps at bond dimension " << bond_dimension << endl;
            if(current_noise > current_noise*1E-01) current_noise = 0;
            else
            {
                current_noise = p.original_noise;
                save_ground_state(ground_state, n0, bond_dimension);
                bond_dimension = int( bond_dimension * p.scaling_bond_dimension );
            }
        }
        else
        {
            save_ground_state(ground_state, n0, bond_dimension);
            bond_dimension = int( bond_dimension * p.scaling_bond_dimension );
            if(number_sweeps_done == 1)
            {
                number_sweeps_in_a_row += 1;
                current_noise *= 1E-02;
            }
            else if(number_sweeps_in_a_row >= 1) current_noise = 0;
            else
            {
                number_sweeps_in_a_row = 0;
                current_noise = p.original_noise;
            }
        }
    } while( number_sweeps_in_a_row < 2 && bond_dimension < p.max_bond_dimension);

    if(bond_dimension > p.max_bond_dimension) throw ITError("perform_dmrg: reached max_bond_dimension without convergence");

    bond_dimension = int( bond_dimension / p.scaling_bond_dimension );
    write_energy_and_variance(save_file_ener, ground_state, H, bond_dimension, ground_state_energy);
    return bond_dimension;
}


int
perform_dmrg_soft(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{
    int n0 = physical_args.getInt("n0");
    DmrgParameters p = read_dmrg_parameters(numerical_args);
    int bond_dimension = p.bond_dimension;

    ofstream save_file_DMRG( make_output_name("delta_energy", physical_args, true) );
    ofstream save_file_ener( make_output_name("energy", physical_args, true) );
    save_file_DMRG << setprecision(set_output_precision);
    save_file_ener << setprecision(set_output_precision);

    double ground_state_energy = prepare_initial_state(ground_state, H);

    int number_sweeps_total = 1;
    double current_noise    = p.original_noise;
    double delta_energy, sum_energy;

    do
    {
        // at most 3 rounds of sweeps at fixed bond dimension
        delta_energy = 0;
        sum_energy   = 0;
        int number_sweeps_done = 0;
        do
        {
            double energy = run_dmrg_sweeps(ground_state, H, p.number_sweep_fixed_bond_dimension, bond_dimension, p.lower_bound_singular_values, current_noise);
            delta_energy = energy - ground_state_energy;
            sum_energy   = energy + ground_state_energy;
            ground_state_energy = energy;

            save_file_DMRG << number_sweeps_total*p.number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy << endl;
            number_sweeps_done += 1;
        } while ( abs(delta_energy/sum_energy) > p.precision_dmrg && number_sweeps_done < 3 );

        if(number_sweeps_done >= 3 && abs(delta_energy/sum_energy) > p.precision_dmrg)
        {
            // not converged: first switch the noise off, then increase the bond dimension
            if(current_noise > current_noise*1E-01) current_noise = 0;
            else
            {
                cerr << "Not converged in 3 rounds of sweeps without noise: increasing the bond dimension" << endl;
                current_noise = p.original_noise;
                save_ground_state(ground_state, n0, bond_dimension);
                bond_dimension = int( bond_dimension * p.scaling_bond_dimension );
            }
        }
    } while( abs(delta_energy/sum_energy) > p.precision_dmrg && bond_dimension < p.max_bond_dimension);

    save_ground_state(ground_state, n0, bond_dimension);
    if(bond_dimension > p.max_bond_dimension) throw ITError("perform_dmrg_soft: reached max_bond_dimension without convergence");

    write_energy_and_variance(save_file_ener, ground_state, H, bond_dimension, ground_state_energy);
    return bond_dimension;
}


void
perform_dmrg_meanfield(MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{
    DmrgParameters p = read_dmrg_parameters(numerical_args);
    const int bond_dimension = 1;

    ofstream save_file_DMRG( make_output_name("meanfield_delta_energy", physical_args, false) );
    ofstream save_file_ener( make_output_name("meanfield_energy", physical_args, false) );
    save_file_DMRG << setprecision(set_output_precision);
    save_file_ener << setprecision(set_output_precision);

    double ground_state_energy = prepare_initial_state(ground_state, H);

    int number_sweeps_total    = 1;
    int number_sweeps_in_a_row = 0;

    // stop after 3 consecutive bond dimensions converged with a single round of sweeps
    do
    {
        double delta_energy, sum_energy;
        int number_sweeps_done = 0;
        do
        {
            double energy = run_dmrg_sweeps(ground_state, H, p.number_sweep_fixed_bond_dimension, bond_dimension, p.lower_bound_singular_values, 0.);
            delta_energy = energy - ground_state_energy;
            sum_energy   = energy + ground_state_energy;
            ground_state_energy = energy;

            save_file_DMRG << number_sweeps_total*p.number_sweep_fixed_bond_dimension << " " << bond_dimension << " " << delta_energy << " " << delta_energy/sum_energy << endl;
            number_sweeps_done  += 1;
            number_sweeps_total += 1;
        } while ( abs(delta_energy/sum_energy) > p.precision_dmrg && number_sweeps_done < 5 );

        number_sweeps_in_a_row = (number_sweeps_done == 1) ? number_sweeps_in_a_row + 1 : 0;
    } while( number_sweeps_in_a_row < 3 );

    write_energy_and_variance(save_file_ener, ground_state, H, bond_dimension, ground_state_energy);
}


double
perform_dmrg_variance(double energy_target, MPS * ground_state , const MPO H, const SiteSet sites, const int set_output_precision, Args const& physical_args, Args const& numerical_args)
{
    DmrgParameters p = read_dmrg_parameters(numerical_args);

    ofstream save_file_DMRG( make_output_name("excited_states_delta_energy", physical_args, true) , ios::app);
    ofstream save_file_ener( make_output_name("energy", physical_args, true) , ios::app);
    save_file_DMRG << setprecision(set_output_precision);
    save_file_ener << setprecision(set_output_precision);

    double ground_state_energy = prepare_initial_state(ground_state, H);

    // fixed schedule: bond dimension 10 -> 50, alternating noise
    double noise = p.original_noise;
    int    max_dims[8] = {10,10,20,20,50,50,50,50};
    int    min_dims[8] = {10,10,10,20,30,40,50,50};
    double noises[8]   = {noise, noise*1E-2, noise, noise*1E-2, noise, noise*1E-2, 0, 0};

    save_file_DMRG << energy_target;
    for(int k = 0 ; k < 8 ; k++)
    {
        double energy = run_dmrg_sweeps(ground_state, H, p.number_sweep_fixed_bond_dimension, max_dims[k], p.lower_bound_singular_values, noises[k], min_dims[k]);
        save_file_DMRG << " " << (energy - ground_state_energy)/(energy + ground_state_energy);
        ground_state_energy = energy;
    }
    save_file_DMRG << endl;

    return write_energy_and_variance(save_file_ener, ground_state, H, int(p.bond_dimension/p.scaling_bond_dimension), ground_state_energy);
}


// initialization state before DMRG calculation over the variances


void
set_excited_state_guess( MPS *ground_state_variance, const SiteSet sites, const int size, const double energy_target, const int cut_off_fock_space)
{
    int number_particles = (energy_target >= 0) ? int(energy_target) : 0;

	std::default_random_engine generator;
 	std::uniform_real_distribution<double> distribution(0, size-1);

	// start from the vacuum, keeping the link indices of *ground_state_variance
	set_vacuum_state( ground_state_variance, sites );

	// current_occupation[j] - 1 bosons on site j+1
	vector<int> current_occupation(size, 1);

	// add the bosons one at a time on random sites, without exceeding the cutoff
	for(int i=1; i <= number_particles ; i++)
	{
		int target_site;
		do{
			target_site	 = int(distribution(generator));
		}while(current_occupation[target_site]>cut_off_fock_space);

		set_site_occupation( ground_state_variance, sites, target_site + 1, current_occupation[target_site] );
		current_occupation[target_site] +=1 ;
	}

    (*ground_state_variance).position(1);
	(*ground_state_variance) /= norm((*ground_state_variance));

}


MPS 
find_fermi_sea(MPO H, const SiteSet sites, const int Nupfill, const int Ndnfill, Sweeps sweeps, double min_varH)
{
    int N = length(sites);
    MPS psi0 = make_electron_product_state(sites, Nupfill, Ndnfill);

    auto [energy,psi_gs] = dmrg(H,psi0,sweeps,{"Quiet",true});
    double E    = inner(psi_gs,H,psi_gs);
    double varE = inner(H,psi_gs,H,psi_gs) - E*E;
    if(abs(varE) > min_varH)
        throw ITError(tinyformat::format("find_fermi_sea: energy variance %.3e larger than %.3e", varE, min_varH));

    // check the filling
    auto total_number = [&](const string& n)
    {
        AutoMPO ampo = AutoMPO(sites);
        for(int j : range1(N)) ampo += 1, n, j;
        return real(innerC(psi_gs, toMPO(ampo), psi_gs));
    };
    if( abs(Nupfill - total_number("Nup")) > 1E-7 || abs(Ndnfill - total_number("Ndn")) > 1E-7 )
        throw ITError("find_fermi_sea: the ground state does not have the requested filling");

    return psi_gs;
}
