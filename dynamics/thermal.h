/**
 * @file thermal.h
 * @brief Thermal (Gibbs) states as MPOs, prepared by imaginary-time evolution of the
 *        infinite-temperature state.
 *
 * Typical use: the state a closed system thermalizes to after a quench from |psi_0> is the Gibbs
 * state with the same energy, find_thermal_state(H, <psi_0|H|psi_0>, dbeta).
 */
#ifndef MYTN_DYNAMICS_THERMAL_H
#define MYTN_DYNAMICS_THERMAL_H

#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;

/** @brief Result of find_thermal_state. */
struct ThermalState
{
    MPO    rho;                 ///< rho(beta) = e^{-beta H} / Tr e^{-beta H} (indices s, s')
    double beta;                ///< inverse temperature reached
    double energy;              ///< Tr(rho H)
    vector<double> betas;       ///< inverse temperature after every step (starting from 0)
    vector<double> energies;    ///< Tr(rho H) after every step
};

/**
 * @brief Thermal state with energy target_energy: starting from rho(0) ~ Id, rho -> e^{-dbeta H} rho e^{-dbeta H}
 *        (first-order MPO of e^{-dbeta H}, toExpH) until Tr(rho H) reaches target_energy.
 *
 * The step that crosses the target is redone with the fraction of dbeta obtained by linear interpolation
 * of the energy, so the energy matches the target up to O(dbeta^2); the remaining error is the
 * first-order error of the imaginary-time steps.
 * @param hamiltonian   Terms of H.
 * @param target_energy Energy (not density) of the thermal state, e.g. <psi_0|H|psi_0> after a quench.
 * @param dbeta         Imaginary-time step (beta grows by 2 dbeta per step).
 * @param args          Truncation of the MPO products, "Cutoff" and "MaxDim"; "MaxBeta" (default 100) stops
 *                      the search for unreachable (too low) energies.
 * @throws ITError if target_energy is above the infinite-temperature energy Tr(H)/Tr(Id) (negative
 *         temperature) or not reached before MaxBeta.
 */
ThermalState
find_thermal_state(const AutoMPO& hamiltonian, const double target_energy, const double dbeta, const Args& args);

#endif
