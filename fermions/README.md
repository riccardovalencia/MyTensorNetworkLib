# fermions

Spinful fermions (ITensor `Electron` sites). Include with `#include "fermions.h"`.

| Header | Content |
|---|---|
| `H_fermionic.h` | `H_tight_binding_electrons`: MPO of the tight-binding chain with on-site fields. |
| `H_single_particle.h` | Single-particle matrices of quadratic models (homogeneous chain, impurity model). |
| `initial_state_fermions.h` | Computational-basis electron states and Fermi sea via DMRG. |

Free-fermion TEBD gates (`gates_free_spinful_fermions`) and the Kondo-like impurity gates are in [spin_boson](../spin_boson).
