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

// Quench in the Ising chain with longitudinal (hx) and transverse (hz) fields,
//   H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ],
// starting from a product state polarized along x, and full counting statistics of the
// magnetization S^x_A of a block A of l = 1..N/2 sites centered in the chain (arXiv:2005.01679).
//
// Every t_measure the program writes to data/:
//   <root>_entropy.txt    t, S_1, ..., S_{N-1}  (entanglement entropy across each bond, natural log)
//   <root>_gf_t<t>.txt    theta, Re G_1, Im G_1, ..., Re G_{N/2}, Im G_{N/2}
// with the generating function G_l(theta) = <exp(i theta S^x_A)> for a block of l sites.
//
// Usage: ./ising_quench input.txt
//   input parameters (with defaults in the code): N, J, hx, hz, T, dt, maxDim, state
//        state = up (all |+x>, default), down (all |-x>) or wall (domain wall)

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N      = input.getInt("N", 16);
    double J      = input.getReal("J", 1.);
    double hx     = input.getReal("hx", 0.1);
    double hz     = input.getReal("hz", 1.);
    double T      = input.getReal("T", 5.);
    double dt     = input.getReal("dt", 0.01);
    int    maxDim = input.getInt("maxDim", 128);
    string state  = input.getString("state", "up");

    double t_measure    = 0.5;    // time between measurements
    int    numberPoints = 100;    // values of theta in [-pi, pi)
    int    maxLength    = N/2;    // largest block

    Args args = {"Cutoff=", 1E-16, "MaxDim=", maxDim, "Verbose=", false};

    // ---------------------------------
    // Initial product state along x

    SpinHalf sites = SpinHalf(N, {"ConserveQNs=", false});

    string config;
    if(state == "up")        config = string(N, '0');
    else if(state == "down") config = string(N, '1');
    else if(state == "wall") config = string(N/2, '0') + string(N - N/2, '1');
    else { cerr << "Unknown state " << state << " (use up, down or wall)\n"; return 1; }

    MPS psi = make_product_state(sites, config, "x");

    // ---------------------------------
    // Second-order Trotter gates: forward sweep with dt/2, then the reversed sweep

    vector<TebdGate> gates = make_ising_gates(sites, N, J, hx, hz, dt);

    // ---------------------------------
    // Output

    fs::create_directories("data");
    string root = tinyformat::format("data/ising_quench_N%d_J%.2f_hx%.2f_hz%.2f_D%d_%s", N, J, hx, hz, maxDim, state);

    ofstream out_entropy(root + "_entropy.txt");
    out_entropy << setprecision(10) << "# t . S_1 . ... . S_{N-1}\n";

    vector<double> theta = make_theta_grid(numberPoints);

    // ---------------------------------
    // Time evolution

    int n_measure         = int(T / t_measure + 1E-9);
    int steps_per_measure = int(t_measure / dt + 1E-9);
    for(int n = 0 ; n <= n_measure ; n++)
    {
        double t = n * t_measure;
        for(int k = 0 ; n > 0 && k < steps_per_measure ; k++)
        {
            psi = tebd_step(psi, gates, args);
            psi.position(1);
            psi.normalize();
        }

        out_entropy << t;
        for(int b = 1 ; b < N ; b++) out_entropy << " " << compute_entanglement_entropy(&psi, b, true);
        out_entropy << endl;

        vector<vector<complex<double> > > G;

        for(int l = 1 ; l <= maxLength ; l++) G.push_back(compute_generating_function(&psi, sites, l, theta));
        write_generating_function(tinyformat::format("%s_gf_t%.2f.txt", root, t), theta, G);

        cerr << "t = " << t << "  maxD = " << maxLinkDim(psi) << "  S(N/2) = " << compute_entanglement_entropy(&psi, N/2, true) << "\n";
    }

    return 0;
}
