/**
 * @file time_evolution.cc
 * @brief Implementation of time_evolution.h (the functions are documented in the header).
 */
#include "time_evolution.h"
#include "lindblad.h"
#include <itensor/all.h>
#include <vector>

using namespace std;
using namespace itensor;


MPS
tebd_step(MPS psi, const vector<TebdGate>& gates, const vector<DissipativeGate>& dissipative_gates, const Args args)
{
    psi = apply_gates(psi, gates, args);
    if(dissipative_gates.empty()) return psi;

    for(DissipativeGate gate : dissipative_gates) psi = apply_dissipative_gate(psi, gate, args);
    return apply_gates(psi, gates, args);
}


MPS
tebd_step(MPS psi, const vector<TebdGate>& gates, const Args args)
{
    return apply_gates(psi, gates, args);
}
