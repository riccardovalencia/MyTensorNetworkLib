/**
 * @file dmrg.h
 * @brief Ground states with DMRG: bond-dimension ramp, noise, convergence check and random restarts.
 *
 * find_ground_state runs two-site DMRG (ITensor dmrg) with a fixed schedule:
 * - sweep k has maximal bond dimension min(max_dim, 10 * 2^(k-1));
 * - noise * 10^-(k-1) is added in the first 4 sweeps (sweeps 1-4), none afterwards;
 * - the local eigensolver (Davidson) does up to 10 iterations per bond (ITensor's default of 2 leaves
 *   DMRG stuck in mixtures of quasi-degenerate states, e.g. near a symmetry-breaking transition);
 * - DMRG stops at the first sweep with the full max_dim, no noise and
 *   |E_k - E_(k-1)| < tolerance * max(1, |E_k|), or after max_sweeps sweeps.
 * The search is repeated number_restarts times from different initial states and the state with
 * the lowest energy is returned: restarts get DMRG out of local minima (metastable states,
 * symmetry-broken states, small bond dimensions), which a single run cannot detect.
 */
#ifndef MYTN_GROUND_STATE_DMRG_H
#define MYTN_GROUND_STATE_DMRG_H

#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;

/** @brief Parameters of find_ground_state (the defaults suit short-range chains of up to ~100 sites). */
struct DmrgParameters
{
    int    max_dim         = 200;     ///< maximal bond dimension
    double cutoff          = 1E-12;   ///< truncation: discarded weight of each SVD
    double noise           = 1E-6;    ///< noise of the first sweep (0: no noise), divided by 10 at each of the next 3 sweeps
    double tolerance       = 1E-10;   ///< convergence: relative energy change between two sweeps
    int    max_sweeps      = 30;      ///< maximal number of sweeps of each run
    int    number_restarts = 1;       ///< independent runs from random states (without conserved quantities)
    int    seed            = 0;       ///< seed of the random initial states (0: different at every call)
};

/** @brief Result of find_ground_state. */
struct GroundState
{
    MPS    psi;                       ///< lowest-energy state found (normalized)
    double energy;                    ///< <psi|H|psi>
    double variance;                  ///< <psi|H^2|psi> - <psi|H|psi>^2 (zero for an exact eigenstate)
    bool   converged;                 ///< the returned run met the tolerance within max_sweeps
    vector<double> restart_energies;  ///< final energy of every run, in order
};

/**
 * @brief DMRG parameters from an input file: keys dmrg_max_dim, dmrg_cutoff, dmrg_noise,
 *        dmrg_tolerance, dmrg_max_sweeps, dmrg_restarts, dmrg_seed (missing keys keep the defaults).
 */
DmrgParameters
read_dmrg_parameters(InputGroup& input);

/**
 * @brief Ground state of H without conserved quantities, starting each run from a random MPS of
 *        bond dimension min(max_dim, 10).
 * @throws ITError if the site indices of H carry quantum numbers (use the InitState overload).
 */
GroundState
find_ground_state(const MPO& H, const DmrgParameters& parameters = DmrgParameters());

/**
 * @brief Ground state of H in the quantum-number sector of the product state initial_state.
 *
 * DMRG preserves the quantum numbers, so every run starts from initial_state and number_restarts
 * is ignored; the noise is what lets the bond dimension grow out of the product state.
 */
GroundState
find_ground_state(const MPO& H, const InitState& initial_state, const DmrgParameters& parameters = DmrgParameters());

/**
 * @brief Ground state of an electron Hamiltonian H at fixed filling (sites with conserved
 *        quantum numbers), from the product state of make_electron_product_state.
 * @param H        Hamiltonian (e.g. make_tight_binding_mpo).
 * @param sites    Electron site set.
 * @param Nupfill  Number of up electrons.
 * @param Ndnfill  Number of down electrons.
 * @throws ITError if the ground state does not have the requested filling.
 */
GroundState
find_fermi_sea(const MPO& H, const SiteSet& sites, const int Nupfill, const int Ndnfill, const DmrgParameters& parameters = DmrgParameters());

#endif
