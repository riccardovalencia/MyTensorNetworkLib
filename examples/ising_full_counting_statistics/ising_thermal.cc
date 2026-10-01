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

// Thermal state of the Ising chain H = -J sum_j [ X_j X_{j+1} + hx X_j + hz Z_j ] with the same
// energy as the product state |+x...+x>, i.e. the state the quench of ising_quench thermalizes to,
// and the generating function of the block magnetization in that state (arXiv:2005.01679).
//
// rho(beta) = exp(-beta H/2) rho(0) exp(-beta H/2) / Tr(...) is obtained from the infinite-temperature
// state rho(0) ~ Id by imaginary-time steps (first order in dbeta), until the energy density reaches
// that of |+x...+x>, -J ((N-1)/N + hx).
//
// Output (data/): <root>_energy.txt (beta, energy density) and <root>_gf.txt (same format as the
// generating functions of ising_quench).
//
// Usage: ./ising_thermal input.txt
//   input parameters (with defaults in the code): N, J, hx, hz, dbeta

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
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N     = input.getInt("N", 16);
    double J     = input.getReal("J", 1.);
    double hx    = input.getReal("hx", 0.1);
    double hz    = input.getReal("hz", 1.);
    double dbeta = input.getReal("dbeta", 0.001);

    int numberPoints = 100;
    int maxLength    = N/2;
    Args args_mult   = {"MaxDim", 1000, "Cutoff", 1E-14};

    double energy_target = -J * ((N - 1.) / N + hx);   // energy density of |+x...+x>

    // ---------------------------------
    // Hamiltonian and exp(-dbeta H)

    SpinHalf sites = SpinHalf(N, {"ConserveQNs=", false});

    auto ampo = AutoMPO(sites);
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

    auto ampo_id = AutoMPO(sites);
    for(int j = 1 ; j <= N ; j++) ampo_id += "Id", j;
    MPO rho = toMPO(ampo_id, {"Exact=", true});
    rho /= trace(rho);

    // ---------------------------------
    // Imaginary-time evolution: rho -> expH rho expH, until the target energy is reached

    fs::create_directories("data");
    string root = tinyformat::format("data/ising_thermal_N%d_J%.2f_hx%.2f_hz%.2f", N, J, hx, hz);

    ofstream out_energy(root + "_energy.txt");
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

    vector<double> theta = {-M_PI};
    for(int k = 0 ; k < numberPoints-1 ; k++) theta.push_back(theta.back() + theta_step(k, numberPoints));

    vector<vector<double> > G_re(maxLength), G_im(maxLength);
    for(int l = 1 ; l <= maxLength ; l++)
        generating_function_sim_size(G_re[l-1], G_im[l-1], l-1, N, numberPoints, &rho, sites);
    write_generating_function(root + "_gf.txt", theta, G_re, G_im);

    return 0;
}
