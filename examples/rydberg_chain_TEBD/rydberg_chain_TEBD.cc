#include <itensor/all.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include "../../spin_boson.h"

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
// (Clean version, without disorder, of the model in Eq. 1 of https://arxiv.org/abs/2309.12392)
//
// Usage: ./rydberg_chain_TEBD [N] [M] [V2] [Omega] [T] [dt] [maxDim]
// Output: data/<file_root>.txt     -> t, fidelity with initial state, half-chain entropy, max bond dimension
//         data/<file_root>_nj.txt  -> t, Rydberg density n_j on each site

int main(int argc, char* argv[])
{
    // Hamiltonian and simulation parameters (defaults allow running with no arguments)
    int    N      = argc > 1 ? atoi(argv[1]) : 12;
    int    M      = argc > 2 ? atoi(argv[2]) : 2;      // number of initial consecutive excitations
    double V2     = argc > 3 ? atof(argv[3]) : 2.;
    double Omega  = argc > 4 ? atof(argv[4]) : 0.1;
    double T      = argc > 5 ? atof(argv[5]) : 10.;
    double dt     = argc > 6 ? atof(argv[6]) : 0.05;
    int    maxDim = argc > 7 ? atoi(argv[7]) : 64;

    double V1        = 1.;
    double cut_off   = 1E-12;
    int steps_measure = 10;
    int total_steps  = int(T / dt);

    Args TEBD_args = {"Cutoff=", cut_off, "MaxDim=", maxDim};

    // ---------------------------------
    // Sites and initial state (1: Rydberg - 0: ground)

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", false});

    vector<int> initial_state;
    for(int j : range1(N)) initial_state.push_back(j <= M ? 1 : 0);

    MPS psi    = initial_computational_state(sites, initial_state);
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

    vector<double> Vj = compute_potential(rj, 6.);

    vector<double> Deltaj, Omegaj;
    for(int j : range1(N))
    {
        Omegaj.push_back(Omega);
        Deltaj.push_back(j % 2 == 0 ? -V1 : -V2);  // anti-blockade
    }

    vector<MyBondGate> gates = gates_rydberg_up_to_VNNN(sites, Deltaj, Omegaj, Vj, dt);

    // ---------------------------------
    // Output files

    fs::create_directories("data");
    string file_root = tinyformat::format("data/rydberg_N%d_M%d_V2_%.2f_Om_%.3f_D%d", N, M, V2, Omega, maxDim);

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
            double EE       = entanglement_entropy(&psi, N/2);

            save_file << t << " " << fidelity << " " << EE << " " << maxLinkDim(psi) << "\n";

            vector<double> mz = measure_magnetization(&psi, sites, "z");
            save_file_nj << t;
            for(double m : mz) save_file_nj << " " << (1-m)/2.;
            save_file_nj << "\n";

            cerr << "t = " << t << "  maxD = " << maxLinkDim(psi) << "  S(N/2) = " << EE << "\n";
        }

        if(k == total_steps) break;

        for(MyBondGate g : gates) psi = apply_gate(psi, g.gate(), g.jn(), TEBD_args);

        psi.position(1);
        psi.normalize();
    }

    save_file.close();
    save_file_nj.close();
    return 0;
}
