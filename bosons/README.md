# bosons

Bosonic quantum east model (bQEM) and bosonic chains. Include with `#include "bosons.h"`.

| Header | Content |
|---|---|
| `build_hamiltonian.h` | MPOs of the bQEM and variants (on-site interactions, hopping, fixed or free `n0`). |
| `TEBD.h` | TEBD gates for closed and open (Lindblad) bQEM dynamics, adiabatic protocols. |
| `perform_dmrg*.h` | DMRG drivers with increasing bond dimension (ground state, mean field, excited states). |
| `initial_state.h` | Fock, coherent, squeezed and cat states; super-bosonic states built from bQEM ground states. |
| `observables.h` | Occupations, Fock-space projectors, covariance matrices, squeezing, variance of H. |
| `scalar_product_*.h` | Overlaps between states with different Fock cutoffs or symmetry sectors. |
| `search_*.h` | Load previously computed states from disk. |
| `external_file.h`, `get_data.h` | Input/output and command-line parsing for the drivers. |

Notes:
- `entanglement_entropy(psi, site)` (in [core](../core)) is in base 2. The former bosons-only copy used the natural logarithm: use `entanglement_entropy(psi, N, site)` for that.
- `initial_state_vacuum_state` and `initialize_excited_state` set the site tensors without link indices (ITensor v2 style); prefer `initial_state_vacuum_state_correct_link`.
