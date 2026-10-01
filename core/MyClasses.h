#ifndef MYTN_CORE_MYCLASSES_H
#define MYTN_CORE_MYCLASSES_H

// Gate containers used by the TEBD routines of the library.
// Sites are 1-indexed, as in ITensor.

#include <itensor/all.h>

using namespace std;
using namespace itensor;

// Unitary gate exp(-i dt h) acting on the (consecutive) sites jn.
// The gate is built at construction: h must carry the unprimed/primed site indices of jn.
class MyBondGate
{
private:
    ITensor gate_;
    vector<int> jn_;
    SiteSet sites_;
public:
    // sites : site set of the MPS the gate acts on
    // j     : sites the gate acts on (e.g. {j,j+1} or {j,j+1,j+2})
    // dt    : time step
    // h     : local Hamiltonian term
	MyBondGate(const SiteSet sites, vector<int> j, const double dt, const ITensor h);
    ITensor gate();        // the gate exp(-i dt h)
    vector<int> jn();      // sites the gate acts on
    void modify_gate(ITensor new_gate);   // overwrite the stored gate (e.g. after re-indexing)
};

// Dissipative gate for the vectorized (bra-ket doubled) Lindblad evolution.
// Stores the first-order term dt*h, acting on sites jket of the ket and jbra of the bra.
class MyBondGateDiss
{
private:
    ITensor gate_;
    vector<int> jnket_;
    vector<int> jnbra_;
    SiteSet sites_;
public:
	MyBondGateDiss(const SiteSet sites, vector<int> jket, vector<int> jbra, const double dt, const ITensor h);
    ITensor gate();          // dt*h (linear approximation of the exponential)
    vector<int> jnket();     // sites acted on in the ket
    vector<int> jnbra();     // sites acted on in the bra
};

// Pair of local operators (Ti on site i, Tj on site j) with a rate gamma.
// Used to build two-site (non-local) jump operators, see gates_nearest_neighbour_local_lindbland.
class MyTrainITensor
{
private:
    ITensor Ti_;
    ITensor Tj_;
    int i_;
    int j_;
    double gamma_;

public:
    // T = {Ti, Tj}, j = {i, j}
	MyTrainITensor(vector<ITensor> T, vector<int> j, double gamma);
    ITensor Ti();
    ITensor Tj();
    int i();
    int j();
    double gamma();
};

#endif
