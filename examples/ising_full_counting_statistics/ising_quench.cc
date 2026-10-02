#include <itensor/all.h>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "mytn.h"

using namespace std;
using namespace itensor;

// Quench in the Ising chain with longitudinal (hx) and transverse (hz) fields,
//   H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ],
// starting from a product state polarized along x, and full counting statistics of the
// magnetization S^x_A of a block A of l = 1..max_block_size sites centered in the chain (arXiv:2005.01679).
//
// Does the state thermalize? The generating function G_l(theta) = <exp(i theta S^x_A)> along the TEBD
// evolution is compared with the one of the thermal state rho ~ exp(-beta H) with the energy of the
// initial state, <psi(0)|H|psi(0)> (find_thermal_state, imaginary-time evolution of the identity).
// For hx != 0 the chain is not integrable and G_l(theta, t) is expected to approach the thermal one
// (up to finite-size fluctuations); for hx = 0 it relaxes to a generalized Gibbs ensemble instead.
//
// Output (data/<run>/):
//   thermal.txt              beta and energy density of the thermal state, energy density of psi(0)
//   thermal_gf.txt           theta, Re G_1, Im G_1, ... of the thermal state
//   thermal_xj.txt           j, <X_j> in the thermal state
//   and every t_measure:
//   xj.txt                   t, <X_1>, ..., <X_N>  (magnetization along x of each spin)
//   entropy.txt              t, S_1, ..., S_{N-1}  (entanglement entropy across each bond, natural log)
//   gf_t<t>.txt              theta, Re G_1, Im G_1, ..., Re G_L, Im G_L  (L = max_block_size)
//   distance_to_thermal.txt  t, D_1, ..., D_L with D_l(t) = max_theta |G_l(theta, t) - G_l^thermal(theta)|
//
// Usage: ./ising_quench input.txt
//   input parameters (with defaults in the code): N, J, hx, hz, state, T, dt, t_measure, max_dim,
//   cut_off, number_points, max_block_size, dbeta and thermal_max_dim (thermal state)
//   state = up (all |+x>, default), down (all |-x>) or wall (domain wall)


int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N       = input.getInt("N", 16);
    double J       = input.getReal("J", 1.);
    double hx      = input.getReal("hx", 0.1);
    double hz      = input.getReal("hz", 1.);
    string state   = input.getString("state", "up");
    double T       = input.getReal("T", 5.);
    double dt      = input.getReal("dt", 0.01);
    double t_measure       = input.getReal("t_measure", 0.5);       // time between measurements
    int    max_dim         = input.getInt("max_dim", 128);          // TEBD
    double cut_off         = input.getReal("cut_off", 1E-16);       // TEBD
    int    number_points   = input.getInt("number_points", 100);    // values of theta in [-pi, pi)
    int    max_block_size  = input.getInt("max_block_size", N/2);   // largest block
    double dbeta           = input.getReal("dbeta", 0.001);         // thermal state: imaginary-time step
    int    thermal_max_dim = input.getInt("thermal_max_dim", 1000); // thermal state: MPO products

    SpinHalf sites = SpinHalf(N, {"ConserveQNs=", false});
    MPS psi = make_product_state(sites, make_standard_config(N, state), "x");
    vector<double> theta = make_theta_grid(number_points);

    string dir = make_run_directory("data", tinyformat::format("ising_quench_N%d_J%.2f_hx%.2f_hz%.2f_%s_T%g_dt%g_D%d_dbeta%g",
                                                               N, J, hx, hz, state, T, dt, max_dim, dbeta), argv[1]);

    // ---------------------------------
    // Thermal state with the energy of the initial state

    AutoMPO H_terms = make_spin_chain_terms(sites, {-J, 0., 0.}, {0., 0., 0.}, {-J*hx, 0., -J*hz});
    double  energy  = real(innerC(psi, toMPO(H_terms), psi));
    ThermalState thermal = find_thermal_state(H_terms, energy, dbeta, {"Cutoff", 1E-14, "MaxDim", thermal_max_dim});
    vector<vector<complex<double> > > G_thermal = compute_block_generating_functions(&thermal.rho, sites, max_block_size, theta);

    write_generating_function(dir + "thermal_gf.txt", theta, G_thermal);
    write_site_values(dir + "thermal_xj.txt", measure_magnetization(thermal.rho, sites, "x"), 13);
    ofstream out_thermal(dir + "thermal.txt");
    out_thermal << setprecision(13) << "# beta . energy density (thermal state) . energy density of psi(0)\n"
                << thermal.beta << " " << thermal.energy / N << " " << energy / N << endl;
    cerr << "thermal state: beta = " << thermal.beta << "  energy density = " << thermal.energy / N
         << "  (psi(0): " << energy / N << ")\n";

    // ---------------------------------
    // TEBD evolution (second-order Trotter steps), compared with the thermal state

    Args args = {"Cutoff=", cut_off, "MaxDim=", max_dim};
    vector<TebdGate> gates = make_ising_gates(sites, J, hx, hz, dt);

    ofstream out_entropy(dir + "entropy.txt");
    out_entropy << setprecision(10) << "# t . S_1 . ... . S_{N-1}\n";
    ofstream out_xj(dir + "xj.txt");
    out_xj << setprecision(10) << "# t . <X_1> . ... . <X_N>\n";
    ofstream out_distance(dir + "distance_to_thermal.txt");
    out_distance << setprecision(10) << "# t";
    for(int l = 1 ; l <= max_block_size ; l++) out_distance << " . D_" << l;
    out_distance << "\n";

    int n_measure         = int(T / t_measure + 1E-9);
    int steps_per_measure = compute_steps_per_measure(t_measure, dt);
    for(int n = 0 ; n <= n_measure ; n++)
    {
        double t = n * t_measure;
        for(int k = 0 ; n > 0 && k < steps_per_measure ; k++)
        {
            psi = tebd_step(psi, gates, args);
            psi.position(1);
            psi.normalize();
        }

        write_row(out_xj, t, measure_magnetization(&psi, sites, "x"));
        write_row(out_entropy, t, compute_entanglement_entropies(&psi, true));

        vector<vector<complex<double> > > G = compute_block_generating_functions(&psi, sites, max_block_size, theta);
        write_generating_function(tinyformat::format("%sgf_t%.2f.txt", dir, t), theta, G);

        vector<double> distances = compute_generating_function_distances(G, G_thermal);
        write_row(out_distance, t, distances);

        cerr << "t = " << t << "  maxD = " << maxLinkDim(psi) << "  S(N/2) = " << compute_entanglement_entropy(&psi, N/2, true)
             << "  D_L = " << distances.back() << "\n";
    }
    return 0;
}
