#include <itensor/all.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <random>
#include "mytn.h"

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// Example: closed dynamics of a 1D Rydberg chain via TEBD.
//
// H = \sum_j Delta_j n_j + \sum_j Omega_j X_j + \sum_j V_j n_j n_{j+1} + (next-nearest-neighbour terms)
//
// Interactions decay as 1/r^6 and are kept up to next-nearest neighbours.
// Atoms sit at alternating distances d1, d2, so that V1 = 1/d1^6 and V2 = 1/d2^6.
// Detunings are set to the anti-blockade condition Delta = -V on even/odd sites.
// The initial state is a "kink": the first M atoms in the Rydberg state, the rest in the ground state.
// Model in Eq. 1 of https://arxiv.org/abs/2309.12392 (sigmax > 0 reproduces Fig. S1 there).
//
// Optional spatial disorder (e.g. finite temperature in the traps): each atom is displaced from its ideal
// position by gaussian noise of width sigmax along the chain, sigmay = sigmax and sigmaz = 5*sigmax
// (transverse trap 5 times weaker). The actual interactions V_j are computed from the displaced positions.
// sigmax = 0 gives the clean chain. seed selects the disorder realization.
//
// Usage: ./rydberg_chain_tebd input.txt
//   input parameters (with defaults in the code): N, M, V2, Omega, T, dt, maxDim, sigmax, seed
// Output: data/<file_root>.txt     -> t, fidelity with initial state, half-chain entropy, max bond dimension
//         data/<file_root>_nj.txt  -> t, Rydberg density n_j on each site
//         data/<file_root>_Vj.txt  -> couplings V_j between sites j and j+1

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    // Hamiltonian and simulation parameters (missing entries take these defaults)
    int    N      = input.getInt("N", 12);
    int    M      = input.getInt("M", 2);      // number of initial consecutive excitations
    double V2     = input.getReal("V2", 2.);
    double Omega  = input.getReal("Omega", 0.1);
    double T      = input.getReal("T", 10.);
    double dt     = input.getReal("dt", 0.05);
    int    maxDim = input.getInt("maxDim", 64);
    double sigmax = input.getReal("sigmax", 0.);     // disorder on atomic positions (units of d1)
    int    seed   = input.getInt("seed", 1);      // disorder realization

    double V1        = 1.;
    double cut_off   = 1E-12;
    int steps_measure = 10;
    int total_steps  = int(T / dt);

    Args TEBD_args = {"Cutoff=", cut_off, "MaxDim=", maxDim};

    // ---------------------------------
    // Sites and initial state (1: Rydberg - 0: ground)

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", false});

    string initial_state = string(M, '1') + string(N - M, '0');   // kink: "11000..."

    MPS psi    = make_product_state(sites, initial_state);
    MPS psi_t0 = psi;

    // ---------------------------------
    // Atomic positions -> interactions V_j between site j and j+1

    double d1 = pow(1/V1, 1./6);
    double d2 = pow(1/V2, 1./6);

    vector<vector<double> > rj;
    double x = 0.;
    for(int j = 0 ; j < N ; j++)
    {
        rj.push_back({x, 0., 0.});
        x += (j % 2 == 0) ? d1 : d2;
    }

    if(sigmax > 0)
    {
        default_random_engine generator;
        generator.seed(seed);
        normal_distribution<double> noise_x(0, sigmax);
        normal_distribution<double> noise_y(0, sigmax);
        normal_distribution<double> noise_z(0, 5*sigmax);

        for(vector<double>& r : rj)
        {
            r[0] += noise_x(generator);
            r[1] += noise_y(generator);
            r[2] += noise_z(generator);
        }
    }

    vector<double> Vj = compute_power_law_couplings(rj, 6.);

    vector<double> Deltaj, Omegaj;
    for(int j : range1(N))
    {
        Omegaj.push_back(Omega);
        Deltaj.push_back(j % 2 == 0 ? -V1 : -V2);  // anti-blockade
    }

    vector<TebdGate> gates = make_rydberg_gates_nnn(sites, Deltaj, Omegaj, Vj, dt);

    // ---------------------------------
    // Output files

    fs::create_directories("data");
    string file_root = tinyformat::format("data/rydberg_N%d_M%d_V2_%.2f_Om_%.3f_D%d", N, M, V2, Omega, maxDim);
    if(sigmax > 0) file_root += tinyformat::format("_sigmax%.5f_seed%d", sigmax, seed);

    // couplings V_j (also used by rydberg_chain_tebd_exact_diagonalization.py)
    ofstream save_file_V(file_root + "_Vj.txt");
    save_file_V << setprecision(16);
    for(double v : Vj) save_file_V << v << "\n";
    save_file_V.close();

    ofstream save_file(file_root + ".txt");
    save_file << "# t . fidelity . entropy . MaxD\n";
    save_file << setprecision(12);

    ofstream save_file_nj(file_root + "_nj.txt");
    save_file_nj << "# t";
    for(int j : range1(N)) save_file_nj << " . " << j;
    save_file_nj << "\n";
    save_file_nj << setprecision(12);

    // ---------------------------------
    // Time evolution

    for(int k = 0 ; k <= total_steps ; k++)
    {
        double t = k*dt;

        if(k % steps_measure == 0)
        {
            double fidelity = pow(abs(innerC(psi_t0, psi)), 2);
            double EE       = compute_entanglement_entropy(&psi, N/2);

            save_file << t << " " << fidelity << " " << EE << " " << maxLinkDim(psi) << "\n";

            vector<double> mz = measure_magnetization(&psi, sites, "z");
            save_file_nj << t;
            for(double m : mz) save_file_nj << " " << (1-m)/2.;
            save_file_nj << "\n";

            cerr << "t = " << t << "  maxD = " << maxLinkDim(psi) << "  S(N/2) = " << EE << "\n";
        }

        if(k == total_steps) break;

        for(TebdGate g : gates) psi = apply_gate(psi, g.gate(), g.sites(), TEBD_args);

        psi.position(1);
        psi.normalize();
    }

    save_file.close();
    save_file_nj.close();
    return 0;
}
