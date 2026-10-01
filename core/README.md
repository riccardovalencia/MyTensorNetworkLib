# core

Model-independent building blocks, used by all the other modules. Include with `#include "core.h"`.

| Header | Content |
|---|---|
| `MyClasses.h` | `MyBondGate` (unitary gate `exp(-i dt h)` on 1-3 consecutive sites), `MyBondGateDiss` (first-order dissipative gate on the purified state), `MyTrainITensor` (pair of local operators for two-site jump operators). |
| `apply_gate.h` | `apply_gate`: apply a 1-, 2- or 3-site gate to an MPS with SVD truncation. |
| `entanglement.h` | `entanglement_entropy(psi, site)`: von Neumann entropy (base 2) across bond (site, site+1); `entanglement_entropy(psi, N, site)`: same in natural log. |
| `initial_state.h` | `initial_computational_state(sites, "0011", basis)`: spin-1/2 product states from a string of 0/1 in the z (default), x or y basis; `insert_state` / `insert_QN_state` to embed an MPS into another (used to build purified states). |
| `state_manipulation.h` | `swap_gate` (move a site along the MPS), `from_MPS_to_MPDO` (density matrix as MPO), `exctract_reduced_density_matrix`. |

A typical TEBD loop:

```cpp
vector<MyBondGate> gates = gates_rydberg_up_to_VNNN(sites, Deltaj, Omegaj, Vj, dt);  // spin_boson module
for(int step = 0; step < nsteps; step++)
{
    for(MyBondGate g : gates) psi = apply_gate(psi, g.gate(), g.jn(), {"Cutoff=", 1E-12, "MaxDim=", 64});
    psi.position(1);
    psi.normalize();
}
```
