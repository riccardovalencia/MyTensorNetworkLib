#include <itensor/all.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "mytn.h"

using namespace std;
using namespace itensor;

// Thermal state of the Ising chain H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ] with the same
// energy as the product state |+x...+x>, i.e. the state the quench of ising_quench thermalizes to,
// and the generating function of the block magnetization in that state (arXiv:2005.01679).
//
// rho(beta) = exp(-beta H/2) rho(0) exp(-beta H/2) / Tr(...) is obtained from the infinite-temperature
// state rho(0) ~ Id by imaginary-time steps (first order in dbeta), until the energy density reaches
// that of |+x...+x>, -J ((N-1)/N + hx).
//
// Output (data/<run>/): energy.txt (beta, energy density) and gf.txt (same format as the
// generating functions of ising_quench).
//
// Usage: ./ising_thermal input.txt
//   input parameters (with defaults in the code): N, J, hx, hz, dbeta

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N     = input.getInt("N", 16);
    double J     = input.getReal("J", 1.);
    double hx    = input.getReal("hx", 0.1);
    double hz    = input.getReal("hz", 1.);
    double dbeta = input.getReal("dbeta", 0.001);
    int    number_points  = input.getInt("number_points", 100);    // values of theta in [-pi, pi)
    int    max_block_size = input.getInt("max_block_size", N/2);   // largest block
    int    max_dim        = input.getInt("max_dim", 1000);         // MPO products
    double cut_off        = input.getReal("cut_off", 1E-14);

    Args args_mult = {"MaxDim", max_dim, "Cutoff", cut_off};

    double energy_target = -J * ((N - 1.) / N + hx);   // energy density of |+x...+x>

    // ---------------------------------
    // Hamiltonian and exp(-dbeta H)

    SpinHalf sites = SpinHalf(N, {"ConserveQNs=", false});

    AutoMPO ampo(sites);
    for(int j = 1 ; j < N ; j++) ampo += -4. * J, "Sx", j, "Sx", j+1;
    for(int j = 1 ; j <= N ; j++)
    {
        ampo += -2. * J * hx, "Sx", j;
        ampo += -2. * J * hz, "Sz", j;
    }
    MPO H    = toMPO(ampo);
    MPO expH = toExpH(ampo, dbeta);

    // ---------------------------------
    // Infinite-temperature state, normalized

    AutoMPO ampo_id(sites);
    for(int j = 1 ; j <= N ; j++) ampo_id += "Id", j;
    MPO rho = toMPO(ampo_id, {"Exact=", true});
    rho /= trace(rho);

    // ---------------------------------
    // Imaginary-time evolution: rho -> expH rho expH, until the target energy is reached

    string dir = make_run_directory("data", tinyformat::format("ising_thermal_N%d_J%.2f_hx%.2f_hz%.2f_dbeta%g_D%d", N, J, hx, hz, dbeta, max_dim), argv[1]);

    ofstream out_energy(dir + "energy.txt");
    out_energy << setprecision(13) << "# beta . energy density\n";

    double beta   = 0.;
    double energy = trace(rho, H) / N;
    out_energy << beta << " " << energy << endl;

    while(energy > energy_target)
    {
        nmultMPO(rho, prime(expH), rho, args_mult);
        rho.mapPrime(2, 1);
        nmultMPO(expH, prime(rho), rho, args_mult);
        rho.mapPrime(2, 1);
        rho /= trace(rho);

        beta  += 2 * dbeta;
        energy = trace(rho, H) / N;
        out_energy << beta << " " << energy << endl;
    }
    cerr << "beta = " << beta << "  energy density = " << energy << "  (target " << energy_target << ")\n";

    // ---------------------------------
    // Generating function of the block magnetization

    vector<double> theta = make_theta_grid(number_points);

    vector<vector<complex<double> > > G;

    for(int l = 1 ; l <= max_block_size ; l++) G.push_back(compute_generating_function(&rho, sites, l, theta));
    write_generating_function(dir + "gf.txt", theta, G);

    return 0;
}
