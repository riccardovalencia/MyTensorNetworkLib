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
// With Jzzz = 0 the coherent part uses two-site gates, otherwise three-site gates.
//
// Output (data/): <root>.txt (t, maxD, Tr rho), <root>_xj.txt (t, <X_1>, ..., <X_N>),
//                 <root>_Tness<Tness>_z1z1.txt (t, Re, Im, |.| of <Z_1(t) Z_1(0)>)
//
// Usage: ./impurity_dynamics [N] [hx] [Jxx] [Jzzz] [gamma] [Tness] [T] [dt] [maxDim]

int main(int argc, char* argv[])
{
    int    N      = argc > 1 ? atoi(argv[1]) : 10;
    double hx     = argc > 2 ? atof(argv[2]) : 0.5;
    double Jxx    = argc > 3 ? atof(argv[3]) : 0.2;
    double Jzzz   = argc > 4 ? atof(argv[4]) : 0.;
    double gamma  = argc > 5 ? atof(argv[5]) : 0.5;
    double Tness  = argc > 6 ? atof(argv[6]) : 5.;
    double T      = argc > 7 ? atof(argv[7]) : 5.;
    double dt     = argc > 8 ? atof(argv[8]) : 0.05;
    int    maxDim = argc > 9 ? atoi(argv[9]) : 128;

    double Jzz       = -1.;
    double cut_off   = 1E-14;
    double t_measure = 0.2;     // profile measurements during the relaxation
    double t_corr    = 0.05;    // autocorrelation measurements
    bool   dissipative = abs(gamma) > 1E-10;

    Args args = {"Cutoff=", cut_off, "MaxDim=", maxDim, "Verbose=", false, "Normalize=", false};

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
    MPO H = toMPO(ampo);

    auto sweeps = Sweeps(20);
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

    vector<BondGate>   gates_2sites;
    vector<MyBondGate> gates_3sites;
    if(Jzzz == 0.) gates_2sites = gates_coherent_part_spin_dissipative_impurity_model(sites, J_NN, h, Lj, gamma, dt_coherent);
    else           gates_3sites = gates_coherent_part_spin_dissipative_nnn_interactions_impurity_model(sites, J_NN, J_NNN, h, dt_coherent);

    vector<MyBondGateDiss> gates_D;
    if(dissipative) gates_D = gates_dissipative_impurity(sites, Lj, gamma, dt);

    auto coherent_step = [&](MPS& state)
    {
        if(Jzzz == 0.) gateTEvol(gates_2sites, dt, dt, state, args);
        else for(MyBondGate g : gates_3sites) state = apply_gate(state, g.gate(), g.jn(), args);
    };
    auto time_step = [&](MPS& state)
    {
        coherent_step(state);
        if(dissipative)
        {
            for(MyBondGateDiss g : gates_D) state = apply_dissipative_gate(state, g, args);
            coherent_step(state);
        }
    };

    // ---------------------------------
    // Output

    fs::create_directories("data");
    string root = tinyformat::format("data/impurity_N%d_Jxx%.3f_Jzzz%.3f_hx%.3f_gamma%.3f_dt%.4f_D%d", N, Jxx, Jzzz, hx, gamma, dt, maxDim);

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
        time_step(rho);
        if(k % steps_measure != 0) continue;

        double t = k*dt;
        out << t << " " << maxLinkDim(rho) << " " << compute_norm_purified_impurity(&rho) << endl;
        out_xj << t;
        for(complex<double> x : measure_magnetization_impurity_first_site(&rho, "x", true)) out_xj << " " << x.real();
        out_xj << endl;
        cerr << "t = " << t << "  maxD = " << maxLinkDim(rho) << "\n";
    }

    // ---------------------------------
    // 3. Autocorrelation <Z_1(t) Z_1(0)>: apply Z on the first ket site and keep evolving

    rho /= compute_norm_purified_impurity(&rho);
    ITensor A = rho(N+1) * 2 * op(sites, "Sz", N+1);
    A.mapPrime(1, 0);
    rho.set(N+1, A);

    ofstream out_corr(tinyformat::format("%s_Tness%.1f_z1z1.txt", root, Tness));
    out_corr << setprecision(14) << "# t . Re . Im . abs of <Z_1(t) Z_1(0)>\n";

    int steps_corr = max(1, int(t_corr / dt + 1E-9));
    total_steps    = int(T / dt + 1E-9);
    for(int k = 1 ; k <= total_steps ; k++)
    {
        time_step(rho);
        if(k % steps_corr != 0) continue;

        complex<double> c = measure_magnetization_impurity_first_site(&rho, "z", false, 1)[0];
        out_corr << k*dt << " " << c.real() << " " << c.imag() << " " << abs(c) << endl;
    }

    return 0;
}
