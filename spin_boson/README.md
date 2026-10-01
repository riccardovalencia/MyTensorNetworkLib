# spin_boson

Spin-1/2 chains, Rydberg arrays and spin-boson (cavity QED) systems, closed and open. Include with `#include "spin_boson.h"`.

| Header | Content |
|---|---|
| `TEBD.h` | Second-order Trotter gates for: generic nearest-neighbour spin chains (`gates_spin_model`), Rydberg chains up to nearest (`gates_rydberg_up_to_VNN`) and next-nearest neighbour (`gates_rydberg_up_to_VNNN`) interactions, PXP, local fields, free spinful fermions; local Lindblad dissipators on the purified state; `TEBD_lindbland_time_evolve`. |
| `TEBD_long_range.h` | Collective light-matter models (one boson coupled to all spins): `gates_photon_matter` (Dicke / Tavis-Cummings, optional spin-spin interactions), `gates_tavis_cummings`, Lindblad evolution with long-range gates. |
| `TEBD_edge_dissipation.h` | Impurity problems with dissipation on the first physical site, on the unfolded density matrix; `doubling_space_gates` maps physical gates to the bra-ket chain. |
| `custom_siteset.h` | `custom_spin_boson` (boson + N-1 spins) and its doubled bra-ket version. |
| `initial_state.h` | Fock state x spin-coherent states. |
| `observables.h` | Magnetizations / occupations, kinks, density-density correlations; traces and local observables on purified states. |
| `parameters.h` | `compute_potential`: couplings `1/r^alpha` from atomic positions. |
| `MyMPO.h` | MPO of the PXP model. |

Examples: [rydberg_chain_TEBD](../examples/rydberg_chain_TEBD), [rydberg_in_leaky_cavity](../examples/rydberg_in_leaky_cavity), [collective_light_matter_unitary](../examples/collective_light_matter_unitary), [collective_light_matter_dissipative](../examples/collective_light_matter_dissipative).
