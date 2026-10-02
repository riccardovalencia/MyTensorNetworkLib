/**
 * @file thermal.cc
 * @brief Implementation of thermal.h (interfaces documented in the header, logic commented here).
 */
#include "thermal.h"
#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;


// normalized identity on the sites of the AutoMPO (infinite-temperature state)
static MPO
make_infinite_temperature_state(const SiteSet& sites)
{
    AutoMPO ampo_id(sites);
    for(int j = 1 ; j <= length(sites) ; j++) ampo_id += "Id", j;
    MPO rho = toMPO(ampo_id, {"Exact=", true});
    rho /= trace(rho);
    return rho;
}


// rho -> e^{-dbeta H} rho e^{-dbeta H} with two MPO products (nmultMPO(A, prime(B)) = B A), normalized
static void
apply_imaginary_time_step(MPO& rho, const MPO& expH, const Args& args)
{
    rho = nmultMPO(rho, prime(expH), args);
    rho.mapPrime(2, 1);
    rho = nmultMPO(expH, prime(rho), args);
    rho.mapPrime(2, 1);
    rho /= trace(rho);
}


// Cool down from beta = 0 in steps of 2 dbeta, recording beta and the energy. The step that crosses the
// target energy is redone with the fraction of dbeta given by linear interpolation of the energy, so
// that the final energy misses the target only by the curvature of E(beta) within one step.
ThermalState
find_thermal_state(const AutoMPO& hamiltonian, const double target_energy, const double dbeta, const Args& args)
{
    MPO H    = toMPO(hamiltonian);
    MPO expH = toExpH(hamiltonian, dbeta);
    double max_beta = args.getReal("MaxBeta", 100.);
    Args args_mult  = {"Cutoff", args.getReal("Cutoff", 1E-14), "MaxDim", args.getInt("MaxDim", 1000)};

    ThermalState thermal = {make_infinite_temperature_state(hamiltonian.sites()), 0., 0., {}, {}};
    thermal.energy = real(traceC(thermal.rho, H));
    if(target_energy > thermal.energy)
        throw ITError(tinyformat::format("find_thermal_state: target energy %g above the infinite-temperature energy %g", target_energy, thermal.energy));

    thermal.betas.push_back(0.);
    thermal.energies.push_back(thermal.energy);
    while(true)
    {
        if(thermal.beta > max_beta)
            throw ITError(tinyformat::format("find_thermal_state: energy %g not reached at beta = %g", target_energy, max_beta));
        MPO rho = thermal.rho;
        apply_imaginary_time_step(rho, expH, args_mult);
        double energy = real(traceC(rho, H));

        if(energy <= target_energy)
        {
            // partial last step
            double fraction = (thermal.energy - target_energy) / (thermal.energy - energy);
            rho = thermal.rho;
            apply_imaginary_time_step(rho, toExpH(hamiltonian, fraction * dbeta), args_mult);
            thermal.rho    = rho;
            thermal.beta  += 2 * fraction * dbeta;
            thermal.energy = real(traceC(rho, H));
            thermal.betas.push_back(thermal.beta);
            thermal.energies.push_back(thermal.energy);
            return thermal;
        }
        thermal.rho    = rho;
        thermal.beta  += 2 * dbeta;
        thermal.energy = energy;
        thermal.betas.push_back(thermal.beta);
        thermal.energies.push_back(thermal.energy);
    }
}
