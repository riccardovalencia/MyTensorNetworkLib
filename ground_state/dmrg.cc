/**
 * @file dmrg.cc
 * @brief Implementation of dmrg.h (interfaces documented in the header, logic commented here).
 */
#include "dmrg.h"
#include "../dof/fermion.h"
#include <itensor/all.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

using namespace std;
using namespace itensor;

static const int initial_dim         = 10;   // bond dimension of the first sweep and of the random states
static const int number_noisy_sweeps = 4;
static const int davidson_iterations = 10;   // ITensor's default (2) leaves DMRG stuck in mixtures of quasi-degenerate states


// ----------------------------------------------------------
// one DMRG run

// sweep k: max_dim min(max_dim, initial_dim 2^(k-1)), noise * 10^-(k-1) for k <= number_noisy_sweeps,
// at most davidson_iterations of the local eigensolver
static Sweeps
make_sweep_schedule(const DmrgParameters& parameters)
{
    Sweeps sweeps(parameters.max_sweeps);
    for(int k = 1 ; k <= parameters.max_sweeps ; k++)
    {
        int dim = initial_dim;
        for(int doubling = 1 ; doubling < k && dim < parameters.max_dim ; doubling++) dim *= 2;
        sweeps.setmaxdim(k, min(dim, parameters.max_dim));
        sweeps.setcutoff(k, parameters.cutoff);
        sweeps.setniter(k, davidson_iterations);
        sweeps.setnoise(k, k <= number_noisy_sweeps ? parameters.noise * pow(10., 1-k) : 0.);
    }
    return sweeps;
}


// stops DMRG at the first sweep with the final bond dimension, no noise and a converged energy
class ConvergenceObserver : public DMRGObserver
{
    public:

    ConvergenceObserver(const MPS& psi, const DmrgParameters& parameters)
        : DMRGObserver(psi), max_dim_(parameters.max_dim), tolerance_(parameters.tolerance) { }

    void measure(Args const&) override { }

    bool checkDone(Args const& args) override
    {
        double energy  = args.getReal("Energy");
        bool final_schedule = args.getInt("MaxDim") == max_dim_ && args.getReal("Noise") == 0.;
        converged = final_schedule && abs(energy - last_energy_) < tolerance_ * max(1., abs(energy));
        last_energy_ = energy;
        return converged;
    }

    bool converged = false;

    private:

    int    max_dim_;
    double tolerance_;
    double last_energy_ = numeric_limits<double>::quiet_NaN();
};


struct DmrgRun
{
    MPS    psi;
    double energy;
    bool   converged;
};


// <psi|H|psi> (real part; psi normalized)
static double
compute_energy(const MPO& H, const MPS& psi)
{
    return real(innerC(psi, H, psi));
}


// one ITensor dmrg call from psi with the sweep schedule, stopped early by the observer
static DmrgRun
run_dmrg(const MPO& H, MPS psi, const DmrgParameters& parameters)
{
    ConvergenceObserver observer(psi, parameters);
    dmrg(psi, H, make_sweep_schedule(parameters), observer, {"Silent", true});
    return {psi, compute_energy(H, psi), observer.converged};
}


// reject parameters that would give no sweep or no run
static void
check_parameters(const DmrgParameters& parameters)
{
    if(parameters.max_dim < 1 || parameters.max_sweeps < 1 || parameters.number_restarts < 1)
        throw ITError("find_ground_state: max_dim, max_sweeps and number_restarts must be positive");
}


// ----------------------------------------------------------
// initial states and result

// unprimed site indices of H
static IndexSet
find_site_indices(const MPO& H)
{
    vector<Index> sites;
    for(int j : range1(length(H)))
        for(const Index& s : siteInds(H, j))
            if(primeLevel(s) == 0) sites.push_back(s);
    return IndexSet(sites);
}


// normalized MPS with gaussian random tensors of bond dimension bond_dim
static MPS
make_random_state(const IndexSet& sites, const int bond_dim, mt19937& generator)
{
    normal_distribution<double> normal(0., 1.);
    MPS psi(sites, bond_dim);
    for(int j : range1(length(psi))) psi.ref(j).generate([&]() { return normal(generator); });
    psi.position(1);
    psi.normalize();
    return psi;
}


// result of the chosen run, with the variance <H^2> - <H>^2 computed once at the end
static GroundState
make_ground_state(const MPO& H, const DmrgRun& run, const vector<double>& restart_energies)
{
    double variance = real(innerC(H, run.psi, H, run.psi)) - run.energy * run.energy;
    return {run.psi, run.energy, variance, run.converged, restart_energies};
}


// ----------------------------------------------------------

// every key falls back to the default of DmrgParameters
DmrgParameters
read_dmrg_parameters(InputGroup& input)
{
    DmrgParameters p;
    p.max_dim         = input.getInt ("dmrg_max_dim",    p.max_dim);
    p.cutoff          = input.getReal("dmrg_cutoff",     p.cutoff);
    p.noise           = input.getReal("dmrg_noise",      p.noise);
    p.tolerance       = input.getReal("dmrg_tolerance",  p.tolerance);
    p.max_sweeps      = input.getInt ("dmrg_max_sweeps", p.max_sweeps);
    p.number_restarts = input.getInt ("dmrg_restarts",   p.number_restarts);
    p.seed            = input.getInt ("dmrg_seed",       p.seed);
    return p;
}


// number_restarts runs, each from a new random state drawn from the same seeded generator;
// keep the run with the lowest energy
GroundState
find_ground_state(const MPO& H, const DmrgParameters& parameters)
{
    check_parameters(parameters);
    if(hasQNs(H)) throw ITError("find_ground_state: H conserves quantum numbers, pass an InitState to fix the sector");

    IndexSet sites = find_site_indices(H);
    mt19937 generator(parameters.seed != 0 ? parameters.seed : random_device{}());

    DmrgRun best = {MPS(), numeric_limits<double>::infinity(), false};
    vector<double> restart_energies;
    for(int restart = 0 ; restart < parameters.number_restarts ; restart++)
    {
        DmrgRun run = run_dmrg(H, make_random_state(sites, min(initial_dim, parameters.max_dim), generator), parameters);
        restart_energies.push_back(run.energy);
        if(run.energy < best.energy) best = run;
    }
    return make_ground_state(H, best, restart_energies);
}


// a single run from the product state (restarts would start from the same state)
GroundState
find_ground_state(const MPO& H, const InitState& initial_state, const DmrgParameters& parameters)
{
    check_parameters(parameters);
    DmrgRun run = run_dmrg(H, MPS(initial_state), parameters);
    return make_ground_state(H, run, {run.energy});
}


// DMRG in the sector of the filled product state, then check <N_up> and <N_dn>
GroundState
find_fermi_sea(const MPO& H, const SiteSet& sites, const int Nupfill, const int Ndnfill, const DmrgParameters& parameters)
{
    GroundState ground_state = find_ground_state(H, make_electron_init_state(sites, Nupfill, Ndnfill), parameters);

    auto total_number = [&](const string& n)
    {
        AutoMPO ampo = AutoMPO(sites);
        for(int j : range1(length(sites))) ampo += 1, n, j;
        return real(innerC(ground_state.psi, toMPO(ampo), ground_state.psi));
    };
    if( abs(Nupfill - total_number("Nup")) > 1E-7 || abs(Ndnfill - total_number("Ndn")) > 1E-7 )
        throw ITError("find_fermi_sea: the ground state does not have the requested filling");

    return ground_state;
}
