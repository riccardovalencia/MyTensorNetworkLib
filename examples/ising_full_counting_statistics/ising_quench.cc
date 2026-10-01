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
// Usage: ./ising_quench [N] [J] [hx] [hz] [T] [dt] [maxDim] [state]
//        state = up (all |+x>, default), down (all |-x>) or wall (domain wall)

// Write the generating function of blocks of 1..G_re.size() sites, one line per theta.
void
write_generating_function(const string file, const vector<double>& theta,
                          const vector<vector<double> >& G_re, const vector<vector<double> >& G_im)
{
    ofstream out(file);
    out << setprecision(10) << "# theta . Re G_l . Im G_l  (l = 1, 2, ...)\n";
    for(size_t k = 0 ; k < theta.size() ; k++)
    {
        out << theta[k];
        for(size_t l = 0 ; l < G_re.size() ; l++) out << " " << G_re[l][k] << " " << G_im[l][k];
        out << "\n";
    }
}

int main(int argc, char* argv[])
{
    int    N      = argc > 1 ? atoi(argv[1]) : 16;
    double J      = argc > 2 ? atof(argv[2]) : 1.;
    double hx     = argc > 3 ? atof(argv[3]) : 0.1;
    double hz     = argc > 4 ? atof(argv[4]) : 1.;
    double T      = argc > 5 ? atof(argv[5]) : 5.;
    double dt     = argc > 6 ? atof(argv[6]) : 0.01;
    int    maxDim = argc > 7 ? atoi(argv[7]) : 128;
    string state  = argc > 8 ? argv[8] : "up";

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

    MPS psi = initial_computational_state(sites, config, "x");

    // ---------------------------------
    // Second-order Trotter gates: forward sweep with dt/2, then the reversed sweep

    vector<BondGate> gates;
    for(int b = 1 ; b <= N-1 ; b++)
    {
        ITensor hterm;
        build_single_step(&hterm, sites, N, J, hx, hz, b);
        gates.push_back(BondGate(sites, b, b+1, BondGate::tReal, dt/2., hterm));
    }
    for(int b = N-1 ; b >= 1 ; b--) { BondGate g = gates[b-1]; gates.push_back(g); }

    // ---------------------------------
    // Output

    fs::create_directories("data");
    string root = tinyformat::format("data/ising_quench_N%d_J%.2f_hx%.2f_hz%.2f_D%d_%s", N, J, hx, hz, maxDim, state);

    ofstream out_entropy(root + "_entropy.txt");
    out_entropy << setprecision(10) << "# t . S_1 . ... . S_{N-1}\n";

    vector<double> theta = {-M_PI};
    for(int k = 0 ; k < numberPoints-1 ; k++) theta.push_back(theta.back() + theta_step(k, numberPoints));

    // ---------------------------------
    // Time evolution

    int n_measure = int(T / t_measure + 1E-9);
    for(int n = 0 ; n <= n_measure ; n++)
    {
        double t = n * t_measure;
        if(n > 0) gateTEvol(gates, t_measure, dt, psi, args);

        out_entropy << t;
        for(int b = 1 ; b < N ; b++) out_entropy << " " << entanglement_entropy(&psi, N, b);
        out_entropy << endl;

        vector<vector<double> > G_re(maxLength), G_im(maxLength);
        for(int l = 1 ; l <= maxLength ; l++)
            generating_function_sim_size(G_re[l-1], G_im[l-1], l-1, N, numberPoints, &psi, sites);
        write_generating_function(tinyformat::format("%s_gf_t%.2f.txt", root, t), theta, G_re, G_im);

        cerr << "t = " << t << "  maxD = " << maxLinkDim(psi) << "  S(N/2) = " << entanglement_entropy(&psi, N, N/2) << "\n";
    }

    return 0;
}
