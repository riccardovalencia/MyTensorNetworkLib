#include <itensor/all.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "mytn.h"

using namespace std;
using namespace itensor;

// Ground state of a spin-1/2 chain with nearest- and next-nearest-neighbour couplings,
//   H = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1})
//     + sum_j (J2x X_j X_{j+2} + J2y Y_j Y_{j+2} + J2z Z_j Z_{j+2}) + sum_j (hx X_j + hz Z_j),
// with find_ground_state (DMRG with random restarts). With conserve_sz = 1 the search is restricted
// to the sector S^z = 0 of the Neel state (needs hx = 0, Jx = Jy, J2x = J2y and N even).
//
// Output (data/<run>/, run name with suffix _sz for conserve_sz = 1):
//                 ground_state.txt, ground_state_restarts.txt  energy, variance, convergence flag;
//                                    energy of every DMRG run (write_ground_state)
//                 profile.txt        j, <Z_j>, <Z_j Z_{j+1}>, entanglement entropy (log2) of the bond (j, j+1)
//
// Usage: ./spin_chain_ground_state input.txt
//   input parameters (with defaults in the code): N, Jx, Jy, Jz, J2x, J2y, J2z, hx, hz, conserve_sz,
//   and the DMRG parameters dmrg_max_dim, dmrg_cutoff, dmrg_noise, dmrg_tolerance, dmrg_max_sweeps,
//   dmrg_restarts, dmrg_seed (ground_state/dmrg.h)

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N   = input.getInt("N", 16);
    double Jx  = input.getReal("Jx", 1.);
    double Jy  = input.getReal("Jy", 1.);
    double Jz  = input.getReal("Jz", 1.);
    double J2x = input.getReal("J2x", 0.);
    double J2y = input.getReal("J2y", 0.);
    double J2z = input.getReal("J2z", 0.);
    double hx  = input.getReal("hx", 0.);
    double hz  = input.getReal("hz", 0.);
    bool   conserve_sz = input.getInt("conserve_sz", 0) == 1;
    DmrgParameters dmrg_parameters = read_dmrg_parameters(input);

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", conserve_sz});
    MPO H = make_spin_chain_mpo(sites, {Jx, Jy, Jz}, {J2x, J2y, J2z}, {hx, 0., hz});

    GroundState ground_state;
    if(conserve_sz)
    {
        InitState neel(sites);
        for(int j = 1 ; j <= N ; j++) neel.set(j, j % 2 == 1 ? "Up" : "Dn");
        ground_state = find_ground_state(H, neel, dmrg_parameters);
    }
    else ground_state = find_ground_state(H, dmrg_parameters);

    MPS psi = ground_state.psi;
    cerr << "E = " << ground_state.energy << "  variance = " << ground_state.variance
         << "  maxD = " << maxLinkDim(psi) << (ground_state.converged ? "" : "  (not converged)") << "\n";

    // ---------------------------------
    // Output

    string run = tinyformat::format("spin_chain_N%d_J%.3f_%.3f_%.3f_J2%.3f_%.3f_%.3f_hx%.3f_hz%.3f", N, Jx, Jy, Jz, J2x, J2y, J2z, hx, hz);
    string dir = make_run_directory("data", conserve_sz ? run + "_sz" : run, argv[1]);

    write_ground_state(dir, ground_state);

    ofstream out_profile(dir + "profile.txt");
    out_profile << setprecision(14) << "# j . <Z_j> . <Z_j Z_{j+1}> . S(j,j+1)\n";
    vector<double> z = measure_magnetization(&psi, sites, "z");
    for(int j = 1 ; j <= N ; j++)
    {
        double zz = 0., entropy = 0.;
        if(j < N)
        {
            ITensor Z_j = make_magnetization_operator(sites(j), "z"), Z_k = make_magnetization_operator(sites(j+1), "z");
            zz      = real(measure_two_point_function(&psi, sites, Z_j, Z_k, j, j+1));
            entropy = compute_entanglement_entropy(&psi, j);
        }
        out_profile << j << " " << z[j-1] << " " << zz << " " << entropy << endl;
    }
    return 0;
}
