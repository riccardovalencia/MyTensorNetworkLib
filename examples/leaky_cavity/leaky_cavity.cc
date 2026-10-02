#include <itensor/all.h>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "mytn.h"

using namespace std;
using namespace itensor;

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
// The initial state is the vacuum times all spins in the coherent state (theta, phi = 0).
//
// Usage: ./leaky_cavity input.txt   (see input.txt; missing entries take the defaults below)
// Output (data/<run>/): observables.txt (t, Tr rho, <X_1>, <Z_1>, <a^dag a>, maxD),
//                      xj.txt, zj.txt (t, <X_j> / <Z_j> for j = 1..N); values divided by Tr rho.

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
    int    max_dim   = input.getInt("max_dim", 1024);
    string coupling = input.getString("coupling", "dicke");
    double theta    = input.getReal("theta", 0.9) * M_PI;   // initial spin state, units of pi
    double t_measure = input.getReal("t_measure", 0.01);    // time between measurements

    double omega0 = 1.;   // energy unit
    double gc = sqrt(0.5 * (abs(h) - V) * (omega0*omega0 + kappa*kappa/4.) / omega0);
    double g  = g_ratio * gc;
    bool dissipative = kappa > 1E-10;
    Args args = {"Cutoff=", cut_off, "MaxDim=", max_dim};

    int steps_measure = compute_steps_per_measure(t_measure, dt);
    int total_steps   = int(T / dt);

    // ---------------------------------
    // Initial state |0> (x) |theta,phi>^N, purified on the doubled chain

    SiteSet sites_single = make_spin_boson_sites(N+1, max_occ);
    MPS psi = make_spin_boson_state(sites_single, 0, theta, 0.);

    SiteSet sites = make_purified_spin_boson_sites(N+1, max_occ);
    MPS rho = randomMPS(sites);
    insert_state(&rho, psi, 1,   true,  true);    // bra on sites 1..N+1 (mirrored, conjugated)
    insert_state(&rho, psi, N+2, false, false);   // ket on sites N+2..2N+2
    rho /= compute_trace_purified(&rho);

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
    // local gates, then the photon-matter gates (the boson travels through the spins with swaps)
    vector<TebdGate> gates_single = make_light_matter_gates(sites_single, omega0, h, g/sqrt(N), dt_coherent, "short-range", coupling, V);
    for(TebdGate gate : make_light_matter_gates(sites_single, omega0, h, g/sqrt(N), dt_coherent, "long-range", coupling)) gates_single.push_back(gate);
    vector<TebdGate> gates = make_purified_gates(gates_single, sites_single, sites);

    vector<DissipativeGate> gates_D;
    if(dissipative) gates_D = make_impurity_dissipative_gates(sites, Lj, kappa, dt);

    // ---------------------------------
    // Output

    string dir = make_run_directory("data", tinyformat::format("leaky_cavity_%s_N%d_maxocc%d_h%.2f_gratio%.2f_V%.2f_kappa%.2f_theta%g_T%g_dt%g_D%d",
                                                               coupling, N, max_occ, h, g_ratio, V, kappa, theta / M_PI, T, dt, max_dim), argv[1]);
    ofstream out(dir + "observables.txt"), out_x(dir + "xj.txt"), out_z(dir + "zj.txt");
    out   << setprecision(8) << "# t . Tr(rho) . <X_1> . <Z_1> . <a^dag a> . maxD\n";
    out_x << setprecision(8) << "# t . <X_1> . ... . <X_N>\n";
    out_z << setprecision(8) << "# t . <Z_1> . ... . <Z_N>\n";

    // physical site q: 1 = boson, j+1 = spin j
    function<double(const ITensor&, int)> expectation = [&](const ITensor& O, int q) { return real(measure_local_operator_purified(&rho, O, false, q)[0]); };

    // ---------------------------------
    // Time evolution

    for(int k = 0 ; k < total_steps ; k++)
    {
        double t = (k+1)*dt;

        rho = tebd_step(rho, gates, gates_D, args);

        if((k+1) % steps_measure != 0) continue;

        double norm = compute_trace_purified(&rho);
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
