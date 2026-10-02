#include <itensor/all.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include "mytn.h"

using namespace std;
using namespace itensor;
namespace fs = std::filesystem;

// <psi_0| Z_j(t) |psi_0> of a spin-1/2 chain,
//   H = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1}) + sum_j (hx X_j + hy Y_j + hz Z_j),
// from a product state, computed in the Schroedinger picture (evolve the state, measure Z_j) and in
// the Heisenberg picture (evolve the MPO of Z_j, take its expectation value on psi_0), with the
// propagator chosen by `propagator`:
//   - "gates": second-order TEBD gates (make_spin_chain_gates), tebd_step / heisenberg_step;
//   - "mpo":   the first-order MPO W = toExpH(H, i dt), applyMPO / heisenberg_step.
// With the same propagator the two pictures agree up to truncation. A field hy != 0 makes H complex,
// so that U^dag Z U and U Z U^dag (evolution backwards in time) give different results.
//
// Output (data/): <root>.txt with t, <Z_j> in the Schroedinger and in the Heisenberg picture, largest
// bond dimension of the state and of the operator.
//
// Usage: ./heisenberg_picture input.txt
//   input parameters (with defaults in the code): N, Jx, Jy, Jz, hx, hy, hz, config, basis, site, T,
//   dt, t_measure, max_dim, cut_off, propagator

int main(int argc, char* argv[])
{
    if(argc != 2) { cerr << "Usage: " << argv[0] << " input.txt\n"; return 1; }

    InputGroup input = InputGroup(argv[1], "input");

    int    N         = input.getInt("N", 10);
    double Jx        = input.getReal("Jx", 1.);
    double Jy        = input.getReal("Jy", 0.);
    double Jz        = input.getReal("Jz", 0.5);
    double hx        = input.getReal("hx", 0.4);
    double hy        = input.getReal("hy", 0.3);
    double hz        = input.getReal("hz", 0.6);
    string config    = input.getString("config", "0110100101");
    string basis     = input.getString("basis", "x");
    int    site      = input.getInt("site", 5);
    double T         = input.getReal("T", 2.);
    double dt        = input.getReal("dt", 0.02);
    double t_measure = input.getReal("t_measure", 0.1);
    int    max_dim   = input.getInt("max_dim", 256);
    double cut_off   = input.getReal("cut_off", 1E-14);
    string propagator = input.getString("propagator", "gates");

    Args args = {"Cutoff=", cut_off, "MaxDim=", max_dim};

    SiteSet sites = SpinHalf(N, {"ConserveQNs=", false});
    MPS psi0 = make_product_state(sites, config, basis);

    AutoMPO ampo_z(sites);
    ampo_z += 2., "Sz", site;
    MPO Z = toMPO(ampo_z);

    if(propagator != "gates" && propagator != "mpo") { cerr << "propagator must be gates or mpo\n"; return 1; }
    vector<double> J = {Jx, Jy, Jz}, h = {hx, hy, hz};
    vector<TebdGate> gates = make_spin_chain_gates(sites, J, h, dt);
    MPO W = toExpH(make_spin_chain_terms(sites, J, {0., 0., 0.}, h), Cplx_i * dt);

    // Schroedinger (psi) and Heisenberg (O) pictures
    MPS psi = psi0;
    MPO O   = Z;

    fs::create_directories("data");
    string root = tinyformat::format("data/heisenberg_%s_N%d_site%d_J%.2f_%.2f_%.2f_h%.2f_%.2f_%.2f_dt%.4f", propagator, N, site, Jx, Jy, Jz, hx, hy, hz, dt);
    ofstream out(root + ".txt");
    out << setprecision(14) << "# t . <Z_j> (Schroedinger) . <Z_j> (Heisenberg) . maxD state . maxD operator\n";

    int steps_measure = compute_steps_per_measure(t_measure, dt);
    int total_steps   = int(T / dt + 1E-9);
    for(int k = 0 ; k <= total_steps ; k++)
    {
        if(k > 0 && propagator == "gates")
        {
            psi = tebd_step(psi, gates, args);
            O   = heisenberg_step(O, gates, args);
        }
        else if(k > 0)
        {
            psi = applyMPO(W, psi, args);
            psi.noPrime("Site");   // applyMPO returns primed site indices
            O   = heisenberg_step(O, W, args);
        }
        if(k % steps_measure != 0) continue;

        out << k*dt << " " << real(innerC(psi, Z, psi)) << " " << real(innerC(psi0, O, psi0))
            << " " << maxLinkDim(psi) << " " << maxLinkDim(O) << endl;
        cerr << "t = " << k*dt << "  maxD state = " << maxLinkDim(psi) << "  maxD operator = " << maxLinkDim(O) << "\n";
    }
    return 0;
}
