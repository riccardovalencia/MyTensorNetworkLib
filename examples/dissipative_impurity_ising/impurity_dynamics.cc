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

// Ising chain with a dephasing impurity on the first site (arXiv:2404.04255):
//   H = -sum_j Z_j Z_{j+1} + Jxx sum_j X_j X_{j+1} + Jzzz sum_j Z_j Z_{j+2} + hx sum_j X_j,
//   L = sqrt(gamma) Z_1.
// The density matrix is purified on 2N sites: the bra on sites 1..N (mirrored, evolving with -H)
// and the ket on N+1..2N (evolving with +H), so that the impurity sits on the central bond:
//
//   o-o-o-o-o-   (ket)
//   |
//   o-o-o-o-o-   (bra)
//
// 1. Ground state of H with DMRG, purified into |rho> = |psi><psi|.
// 2. Lindblad evolution up to Tness (relaxation towards the stationary state); every t_measure
//    the program writes the bond dimension, Tr(rho) and the profile <X_j> (normalized by Tr(rho)).
// 3. Autocorrelation <Z_1(t) Z_1(0)> in the stationary state up to time T (quantum regression
//    theorem: Z_1 is applied to the ket and the state is evolved further).
//
//
// Output (data/): <root>.txt (t, maxD, Tr rho), <root>_xj.txt (t, <X_1>, ..., <X_N>),
//                 <root>_Tness<Tness>_z1z1.txt (t, Re, Im, |.| of <Z_1(t) Z_1(0)>)
//
// Usage: ./impurity_dynamics input.txt
//   input parameters (with defaults in the code): N, hx, Jxx, Jzzz, gamma, Tness, T, dt, max_dim

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N      = input.getInt("N", 10);
    double hx     = input.getReal("hx", 0.5);
    double Jxx    = input.getReal("Jxx", 0.2);
    double Jzzz   = input.getReal("Jzzz", 0.);
    double gamma  = input.getReal("gamma", 0.5);
    double Tness  = input.getReal("Tness", 5.);
    double T      = input.getReal("T", 5.);
    double dt     = input.getReal("dt", 0.05);
    int    max_dim = input.getInt("max_dim", 128);
    double cut_off     = input.getReal("cut_off", 1E-14);     // SVD truncation
    double t_measure   = input.getReal("t_measure", 0.2);     // profile measurements during the relaxation
    double t_corr      = input.getReal("t_corr", 0.05);       // autocorrelation measurements
    int    dmrg_sweeps = input.getInt("dmrg_sweeps", 20);     // ground-state search
    double hz          = input.getReal("hz", 0.);             // symmetry-breaking field, ground-state search only

    double Jzz       = -1.;
    bool   dissipative = abs(gamma) > 1E-10;

    Args args = {"Cutoff=", cut_off, "MaxDim=", max_dim, "Verbose=", false, "Normalize=", false};

    // ---------------------------------
    // Ground state of H

    SiteSet sites_phys = SpinHalf(N, {"ConserveQNs=", false});

    auto ampo = AutoMPO(sites_phys);
    for(int j = 1 ; j < N ; j++)
    {
        ampo += 4 * Jzz, "Sz", j, "Sz", j+1;
        ampo += 4 * Jxx, "Sx", j, "Sx", j+1;
    }
    for(int j = 1 ; j < N-1 ; j++) ampo += 4 * Jzzz, "Sz", j, "Sz", j+2;
    for(int j = 1 ; j <= N ; j++)  ampo += 2 * hx, "Sx", j;
    for(int j = 1 ; j <= N ; j++)  ampo += 2 * hz, "Sz", j;   // only for the ground-state search
    MPO H = toMPO(ampo);

    auto sweeps = Sweeps(dmrg_sweeps);
    sweeps.maxdim() = 10,10,10,20,20,40,40,100,200,200;
    sweeps.cutoff() = 1E-14;
    sweeps.noise()  = 0;
    auto [energy, psi] = dmrg(H, randomMPS(sites_phys), sweeps, {"Quiet", true});
    cerr << "Ground-state energy: " << energy << "\n";

    // ---------------------------------
    // Purified state on 2N sites: bra (mirrored, conjugated) on 1..N, ket on N+1..2N

    SiteSet sites = SpinHalf(2*N, {"ConserveQNs=", false});
    MPS rho = randomMPS(sites);
    insert_state(&rho, psi, 1,   true,  true);
    insert_state(&rho, psi, N+1, false, false);

    // ---------------------------------
    // Gates: dephasing Z on the impurity (first site of bra and ket), coherent part with dt/2
    // on both sides of the dissipative step

    vector<ITensor> Lj = {2. * op(sites, "Sz", N), 2. * op(sites, "Sz", N+1)};
    vector<double> J_NN  = {Jxx, 0., Jzz};
    vector<double> J_NNN = {0., 0., Jzzz};
    vector<double> h     = {hx, 0., 0.};
    double dt_coherent = dissipative ? dt/2. : dt;

    // two-site gates for nearest-neighbour couplings only, three-site gates otherwise
    vector<TebdGate> gates = (Jzzz == 0.) ? make_spin_impurity_gates(sites, J_NN, h, Lj, gamma, dt_coherent)
                                          : make_spin_impurity_nnn_gates(sites, J_NN, J_NNN, h, dt_coherent);

    vector<DissipativeGate> gates_D;
    if(dissipative) gates_D = make_impurity_dissipative_gates(sites, Lj, gamma, dt);

    // ---------------------------------
    // Output

    fs::create_directories("data");
    string root = tinyformat::format("data/impurity_N%d_Jxx%.3f_Jzzz%.3f_hx%.3f_gamma%.3f_dt%.4f_D%d", N, Jxx, Jzzz, hx, gamma, dt, max_dim);

    ofstream out(root + ".txt");
    out << setprecision(14) << "# t . maxD . Tr(rho)\n";
    ofstream out_xj(root + "_xj.txt");
    out_xj << setprecision(14) << "# t . <X_1> . ... . <X_N>\n";

    // ---------------------------------
    // 2. Relaxation up to Tness

    int steps_measure = max(1, int(t_measure / dt + 1E-9));
    int total_steps   = int(Tness / dt + 1E-9);
    for(int k = 1 ; k <= total_steps ; k++)
    {
        rho = tebd_step(rho, gates, gates_D, args);
        if(k % steps_measure != 0) continue;

        double t = k*dt;
        out << t << " " << maxLinkDim(rho) << " " << compute_trace_purified(&rho) << endl;
        out_xj << t;
        for(complex<double> x : measure_magnetization_purified(&rho, "x", true)) out_xj << " " << x.real();
        out_xj << endl;
        cerr << "t = " << t << "  maxD = " << maxLinkDim(rho) << "\n";
    }

    // ---------------------------------
    // 3. Autocorrelation <Z_1(t) Z_1(0)>: apply Z on the first ket site and keep evolving

    rho /= compute_trace_purified(&rho);
    ITensor A = rho(N+1) * 2 * op(sites, "Sz", N+1);
    A.mapPrime(1, 0);
    rho.set(N+1, A);

    ofstream out_corr(tinyformat::format("%s_Tness%.1f_z1z1.txt", root, Tness));
    out_corr << setprecision(14) << "# t . Re . Im . abs of <Z_1(t) Z_1(0)>\n";

    int steps_corr = max(1, int(t_corr / dt + 1E-9));
    total_steps    = int(T / dt + 1E-9);
    for(int k = 1 ; k <= total_steps ; k++)
    {
        rho = tebd_step(rho, gates, gates_D, args);
        if(k % steps_corr != 0) continue;

        complex<double> c = measure_magnetization_purified(&rho, "z", false, 1)[0];
        out_corr << k*dt << " " << c.real() << " " << c.imag() << " " << abs(c) << endl;
    }

    return 0;
}
