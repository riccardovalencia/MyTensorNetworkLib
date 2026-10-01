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
// for hx = hx_min, hx_min + dhx, ..., hx_max, with DMRG (find_ground_state). With hz = 0 the ground
// state is twice degenerate in the ordered phase (small hx) and DMRG returns one of the two states:
// fix dmrg_seed to make the output reproducible.
//
// Output (data/): <root>.txt with hx, energy, energy variance <H^2> - <H>^2 and the magnetization
// <X>, <Y>, <Z> of the central site.
//
// Usage: ./ground_state_scan input.txt
//   input parameters (with defaults in the code): N, Jxx, Jzzz, hx_min, hx_max, dhx, hz,
//   and the DMRG parameters dmrg_* (read_dmrg_parameters, ground_state/dmrg.h)

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
    double hz     = input.getReal("hz", 0.);             // symmetry-breaking longitudinal field
    DmrgParameters dmrg_parameters = read_dmrg_parameters(input);

    double Jzz = -1.;
    int j_meas = N/2;

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", false});

    fs::create_directories("data");
    string root = tinyformat::format("data/ground_state_scan_N%d_Jxx%.3f_Jzzz%.3f", N, Jxx, Jzzz);
    ofstream out(root + ".txt");
    out << setprecision(8) << "# hx . E . <H^2>-<H>^2 . <X_N/2> . <Y_N/2> . <Z_N/2>\n";

    for(double hx = hx_min ; hx <= hx_max + 1E-9 ; hx += dhx)
    {
        MPO H = make_spin_chain_mpo(sites, {Jxx, 0., Jzz}, {0., 0., Jzzz}, {hx, 0., hz});
        GroundState ground_state = find_ground_state(H, dmrg_parameters);
        MPS psi = ground_state.psi;

        out << hx << " " << ground_state.energy << " " << ground_state.variance;
        for(string d : {"x", "y", "z"}) out << " " << measure_magnetization(&psi, sites, d)[j_meas-1];
        out << endl;
        cerr << "hx = " << hx << "  E = " << ground_state.energy << "\n";
    }
    return 0;
}
