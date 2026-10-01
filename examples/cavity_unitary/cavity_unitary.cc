#include <itensor/all.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "mytn.h"

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// Closed dynamics of N spin-1/2 coupled to a single cavity mode a (pure state, no losses):
//   H = omega0 a^dag a + h S^z + (g/sqrt(N)) * coupling,
//   coupling = (a + a^dag) S^x  ("dicke")  or  S^+ a + S^- a^dag  ("tavis").
// The boson is site 1 of the MPS and the spins sites 2..N+1; the boson interacts with every spin
// by travelling through the chain with swap gates. The initial state is the cavity vacuum times
// all spins in the coherent state (theta, phi = 0). Benchmarked with exact diagonalization.
//
// Usage: ./cavity_unitary input.txt
// Output (data/): <root>.txt with t, fidelity |<psi(0)|psi(t)>|^2, <S^x>/N, <S^z>/N (Pauli units),
//                 <a^dag a>/N, maxD.

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N        = input.getInt("N", 3);
    int    max_occ  = input.getInt("max_occ", 6);
    double omega0   = input.getReal("omega0", 1.);
    double h        = input.getReal("h", 0.5);
    double g        = input.getReal("g", 1.5);
    double theta    = input.getReal("theta", 0.5) * M_PI;   // in units of pi
    double T        = input.getReal("T", 50.);
    double dt       = input.getReal("dt", 0.005);
    double cut_off  = input.getReal("cut_off", 1E-8);
    int    maxDim   = input.getInt("maxDim", 50);
    string coupling = input.getString("coupling", "dicke");

    int steps_measure = 10;
    int total_steps   = int(T / dt);

    SiteSet sites = make_spin_boson_sites(N+1, max_occ);
    MPS psi    = make_spin_boson_state(sites, 0, theta, 0.);
    MPS psi_t0 = psi;

    // local gates, then the photon-matter gates: the boson travels through the chain with swaps
    vector<TebdGate> gates = make_light_matter_gates(sites, omega0, h, g/sqrt(N), dt, "short-range", coupling);
    for(TebdGate gate : make_light_matter_gates(sites, omega0, h, g/sqrt(N), dt, "long-range", coupling)) gates.push_back(gate);
    Args args = {"Cutoff=", cut_off, "MaxDim=", maxDim};

    fs::create_directories("data");
    string root = tinyformat::format("data/cavity_unitary_%s_N%d_maxocc%d_omega%.2f_h%.2f_g%.2f", coupling, N, max_occ, omega0, h, g);
    ofstream out(root + ".txt");
    out << "# t . fidelity . <S^x>/N . <S^z>/N . <a^dag a>/N . maxD\n" << setprecision(8);

    for(int k = 0 ; k < total_steps ; k++)
    {
        double t = (k+1)*dt;

        psi = tebd_step(psi, gates, args);
        psi.position(1);
        psi.normalize();

        if((k+1) % steps_measure != 0) continue;

        double fidelity = pow(abs(innerC(psi_t0, psi)), 2);
        vector<double> mx = measure_magnetization(&psi, sites, "x");   // mx[0] = <a^dag a>
        vector<double> mz = measure_magnetization(&psi, sites, "z");
        double Sx = 0., Sz = 0.;
        for(int j = 1 ; j <= N ; j++) { Sx += mx[j]; Sz += mz[j]; }

        out << t << " " << fidelity << " " << Sx/N << " " << Sz/N << " " << mx[0]/N << " " << maxLinkDim(psi) << endl;
        cerr << t << " " << maxLinkDim(psi) << endl;
    }
    return 0;
}
