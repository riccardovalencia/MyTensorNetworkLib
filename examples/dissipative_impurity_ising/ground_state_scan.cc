#include <itensor/all.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include "mytn.h"

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// Ground states of the Ising chain of impurity_dynamics (without the impurity),
//   H = -sum_j Z_j Z_{j+1} + Jxx sum_j X_j X_{j+1} + Jzzz sum_j Z_j Z_{j+2} + hx sum_j X_j,
// for hx = hx_min, hx_min + dhx, ..., hx_max, with DMRG.
//
// Output (data/): <root>.txt with hx, energy, energy variance <H^2> - <H>^2 and the magnetization
// <X>, <Y>, <Z> of the central site.
//
// Usage: ./ground_state_scan input.txt
//   input parameters (with defaults in the code): N, Jxx, Jzzz, hx_min, hx_max, dhx

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N      = input.getInt("N", 20);
    double Jxx    = input.getReal("Jxx", 0.2);
    double Jzzz   = input.getReal("Jzzz", 0.);
    double hx_min = input.getReal("hx_min", 0.5);
    double hx_max = input.getReal("hx_max", 1.5);
    double dhx    = input.getReal("dhx", 0.05);

    double Jzz = -1.;
    int j_meas = N/2;

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", false});

    auto sweeps = Sweeps(20);
    sweeps.maxdim() = 10,10,10,20,20,40,40,100,200,200;
    sweeps.cutoff() = 1E-14;
    sweeps.noise()  = 0;

    fs::create_directories("data");
    string root = tinyformat::format("data/ground_state_scan_N%d_Jxx%.3f_Jzzz%.3f", N, Jxx, Jzzz);
    ofstream out(root + ".txt");
    out << setprecision(8) << "# hx . E . <H^2>-<H>^2 . <X_N/2> . <Y_N/2> . <Z_N/2>\n";

    for(double hx = hx_min ; hx <= hx_max + 1E-9 ; hx += dhx)
    {
        auto ampo = AutoMPO(sites);
        for(int j = 1 ; j <= N ; j++)  ampo += 2 * hx, "Sx", j;
        for(int j = 1 ; j < N ; j++)
        {
            ampo += 4 * Jzz, "Sz", j, "Sz", j+1;
            ampo += 4 * Jxx, "Sx", j, "Sx", j+1;
        }
        for(int j = 1 ; j < N-1 ; j++) ampo += 4 * Jzzz, "Sz", j, "Sz", j+2;
        MPO H = toMPO(ampo);

        auto [energy, psi] = dmrg(H, randomMPS(sites), sweeps, {"Quiet", true});
        double variance = inner(psi, H, H, psi) - energy*energy;

        out << hx << " " << energy << " " << variance;
        for(string d : {"x", "y", "z"}) out << " " << measure_magnetization(&psi, sites, d)[j_meas-1];
        out << endl;
        cerr << "hx = " << hx << "  E = " << energy << "\n";
    }
    return 0;
}
