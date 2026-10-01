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

// N spin-1/2 (e.g. Rydberg atoms) coupled to a single lossy cavity mode a:
//   H = omega0 a^dag a + h S^z + (g/sqrt(N)) * coupling + V sum_j n_j n_{j+1},   n = (1 - Z)/2,
//   coupling = (a + a^dag) S^x  ("dicke")  or  S^+ a + S^- a^dag  ("tavis"),
//   d rho/dt = -i[H, rho] + kappa (a rho a^dag - 1/2 {a^dag a, rho}).
// g is given in units of the mean-field critical coupling g_c = sqrt((|h| - V)(omega0^2 + kappa^2/4)/(2 omega0)).
//
// The density matrix is purified on 2(N+1) sites, bra (mirrored) and ket around the two bosons:
//   s_N ... s_1 b | b s_1 ... s_N.
// Each time step applies the local gates and the photon-matter gates (the boson travels through
// the spins with swap gates), the photon losses on the central bond, and again the coherent gates.
// The initial state is the vacuum times all spins in the coherent state theta = 0.9 pi, phi = 0.
//
// Usage: ./leaky_cavity input.txt   (see input.txt; missing entries take the defaults below)
// Output (data/): <root>_obs.txt (t, Tr rho, <X_1>, <Z_1>, <a^dag a>, maxD),
//                 <root>_xj.txt, <root>_zj.txt (t, <X_j> / <Z_j> for j = 1..N); values divided by Tr rho.

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");
    int    N        = input.getInt("N", 3);
    int    max_occ  = input.getInt("max_occ", 2);
    double h        = input.getReal("h", 1.);
    double g_ratio  = input.getReal("g", 1.6);
    double V        = input.getReal("V", 0.);
    double kappa    = input.getReal("kappa", 1.);
    double T        = input.getReal("T", 15.);
    double dt       = input.getReal("dt", 0.01);
    double cut_off  = input.getReal("cut_off", 1E-14);
    int    maxDim   = input.getInt("maxDim", 1024);
    string coupling = input.getString("coupling", "dicke");

    double omega0 = 1.;   // energy unit
    double gc = sqrt(0.5 * (abs(h) - V) * (omega0*omega0 + kappa*kappa/4.) / omega0);
    double g  = g_ratio * gc;
    bool dissipative = kappa > 1E-10;
    Args args = {"Cutoff=", cut_off, "MaxDim=", maxDim};

    int steps_measure = (dt < 0.01) ? int(0.01/dt) : 1;
    int total_steps   = int(T / dt);

    // ---------------------------------
    // Initial state |0> (x) |theta,phi>^N, purified on the doubled chain

    SiteSet sites_single = custom_spin_boson(N+1, max_occ);
    MPS psi = initialize_spin_boson_state(sites_single, 0, 0.9 * M_PI, 0.);

    SiteSet sites = custom_spin_boson_doubling(N+1, max_occ);
    MPS rho = randomMPS(sites);
    insert_state(&rho, psi, 1,   true,  true);    // bra on sites 1..N+1 (mirrored, conjugated)
    insert_state(&rho, psi, N+2, false, false);   // ket on sites N+2..2N+2
    rho /= compute_norm_purified_impurity(&rho);

    // ---------------------------------
    // Operators: photon number (boson), Pauli X and Z (spin), given on the indices of a physical site

    Index b = sites(N+2), s = sites(1);
    ITensor Nb = ITensor(b, prime(b));
    for(int d = 1 ; d <= dim(b) ; d++) Nb.set(b(d), prime(b)(d), d-1);
    ITensor X = ITensor(s, prime(s));
    X.set(s(1), prime(s)(2), 1.);
    X.set(s(2), prime(s)(1), 1.);
    ITensor Z = ITensor(s, prime(s));
    Z.set(s(1), prime(s)(1),  1.);
    Z.set(s(2), prime(s)(2), -1.);

    // photon losses: a on the boson of the bra and of the ket (central bond)
    vector<ITensor> Lj;
    for(int j : {N+1, N+2})
    {
        Index sj = sites(j);
        ITensor A = ITensor(sj, prime(sj));
        for(int d = 1 ; d < dim(sj) ; d++) A.set(sj(d+1), prime(sj)(d), sqrt(d));
        Lj.push_back(A);
    }

    // ---------------------------------
    // Gates: built on the physical chain and copied on ket and bra

    double dt_coherent = dissipative ? dt/2. : dt;
    vector<BondGate> local_single = gates_photon_matter(sites_single, omega0, h, g/sqrt(N), dt_coherent, "short-range", coupling, V);
    vector<BondGate> pm_single    = gates_photon_matter(sites_single, omega0, h, g/sqrt(N), dt_coherent, "long-range",  coupling);
    vector<MyBondGate> gates_local = doubling_space_gates(local_single, sites_single, sites);
    vector<MyBondGate> gates_pm    = doubling_space_gates(pm_single,    sites_single, sites);

    vector<MyBondGateDiss> gates_D;
    if(dissipative) gates_D = gates_dissipative_impurity(sites, Lj, kappa, dt);

    auto coherent_step = [&]()
    {
        for(MyBondGate gate : gates_local) rho = apply_local_gate_purified(rho, gate, args);
        for(MyBondGate gate : gates_pm)    rho = apply_photon_matter_gate_purified(rho, gate, args);
    };

    // ---------------------------------
    // Output

    fs::create_directories("data");
    string root = tinyformat::format("data/leaky_cavity_%s_N%d_maxocc%d_h%.2f_gratio%.2f_V%.2f_kappa%.2f_D%d",
                                     coupling, N, max_occ, h, g_ratio, V, kappa, maxDim);
    ofstream out(root + "_obs.txt"), out_x(root + "_xj.txt"), out_z(root + "_zj.txt");
    out   << setprecision(8) << "# t . Tr(rho) . <X_1> . <Z_1> . <a^dag a> . maxD\n";
    out_x << setprecision(8) << "# t . <X_1> . ... . <X_N>\n";
    out_z << setprecision(8) << "# t . <Z_1> . ... . <Z_N>\n";

    // physical site q: 1 = boson, j+1 = spin j
    auto expectation = [&](const ITensor& O, int q) { return real(measure_local_obs_impurity_first_site(&rho, O, false, q)[0]); };

    // ---------------------------------
    // Time evolution

    for(int k = 0 ; k < total_steps ; k++)
    {
        double t = (k+1)*dt;

        coherent_step();
        if(dissipative)
        {
            for(MyBondGateDiss gate : gates_D) rho = apply_dissipative_gate(rho, gate, args);
            coherent_step();
        }

        if(k % steps_measure != 0) continue;

        double norm = compute_norm_purified_impurity(&rho);
        out << t << " " << norm << " " << expectation(X, 2)/norm << " " << expectation(Z, 2)/norm
            << " " << expectation(Nb, 1)/norm << " " << maxLinkDim(rho) << endl;
        out_x << t;
        out_z << t;
        for(int j = 1 ; j <= N ; j++)
        {
            out_x << " " << expectation(X, j+1)/norm;
            out_z << " " << expectation(Z, j+1)/norm;
        }
        out_x << endl;
        out_z << endl;
        cerr << t << " " << norm << "\n";
    }
    return 0;
}
