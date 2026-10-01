# spins

Ising chain in longitudinal (`hx`) and transverse (`hz`) magnetic fields, and full counting statistics (FCS) of the subsystem magnetization. Include with `#include "transverse_field_ising_chain.h"`.

| Header | Content |
|---|---|
| `TEBD.h` | `build_single_step`: bond term of `H = -J sum_j [X_j X_{j+1} + hx X_j + hz Z_j]`. |
| `full_counting_statistics.h` | Generating function `G(theta) = <exp(i theta S^x_A)>` on a block A centered in the chain, for pure states (MPS) and mixed states (MPO); first four cumulants (`measuring_moments`). |
| `observables.h` | `measure_mx_mz`: print `<X_j>`, `<Z_j>`. |
| `external_file.h`, `get_data.h` | File names, input/output and command-line parsing for the drivers. |

Initial product states are built with `initial_computational_state` from [core](../core), e.g.
`initial_computational_state(sites, "00001111", "x")` for a domain wall along x.
