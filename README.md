# MyTensorNetworkLib

MyTensorNetworkLib is a C++ library for tensor network methods built on top of the C++ ITensor (v3) library. It provides customized methods for one-dimensional tensor networks (Matrix Product States and Matrix Product Operators): preparation of specific initial states, equilibrium properties (DMRG), and dynamics of closed and open systems via the Time Evolving Block Decimation (TEBD) algorithm. It handles spin, bosonic and fermionic degrees of freedom.

## Prerequisites

- C++ ITensor v3 library, see the [installation guide](https://itensor.org/docs.cgi?vers=cppv3&page=install).
- A C++17-compliant compiler.
- Optional: [Doxygen](https://www.doxygen.nl) to generate the HTML documentation.

## Structure

The library is organized by role. Each folder contains pairs `name.h` / `name.cc`; `mytn.h` includes everything.

| Folder | Header | Content |
|---|---|---|
| [mps](mps/) | `gates.h` | Gate containers (`MyBondGate`, `MyBondGateDiss`, `MyTrainITensor`) and `apply_gate`. |
| | `entanglement.h` | Von Neumann entanglement entropy. |
| | `mps_tools.h` | Embedding an MPS into another (`insert_state`), swap gates, density matrices, two-point functions. |
| [dof](dof/) | `spin_half.h` | Spin-1/2 product states (`initial_computational_state`) and local observables. |
| | `boson.h` | Bosonic Fock, coherent, squeezed and cat states; occupations, Fock-space projectors, squeezing. |
| | `fermion.h` | Product states of spinful fermions. |
| | `spin_boson.h` | Boson + spins site sets (cavity QED) and their purified (bra-ket) version; initial states. |
| [models](models/) | `spin_chain.h` | Nearest-neighbour spin chains, Ising chain in longitudinal and transverse fields. |
| | `rydberg.h` | Rydberg arrays (up to next-nearest-neighbour interactions), PXP model. |
| | `light_matter.h` | Dicke and Tavis-Cummings models (one boson coupled to all spins). |
| | `impurity.h` | Impurity models on the purified density matrix (Kondo-like, dissipative spin impurity). |
| | `tight_binding.h` | Free fermions: MPOs, single-particle matrices, TEBD gates. |
| | `bqem.h`, `bqem_states.h` | Bosonic quantum east model: Hamiltonians, TEBD gates, super-bosonic states, dressed operators. |
| [dynamics](dynamics/) | `lindblad.h` | Dissipative gates for Lindblad dynamics on the vectorized density matrix. |
| | `time_evolution.h` | Time-evolution drivers for purified states (coherent + dissipative gates). |
| | `purified_state.h` | Trace and observables of a density matrix stored as a purified MPS. |
| | `adiabatic.h` | Adiabatic ramps (bosonic quantum east model). |
| [ground_state](ground_state/) | `dmrg.h` | DMRG drivers with increasing bond dimension, excited states, Fermi sea. |
| [analysis](analysis/) | `full_counting_statistics.h` | Generating function and cumulants of the subsystem magnetization. |
| [io](io/) | `output.h` | File names and output files. |
| | `load.h` | Loading stored states. |
| | `input.h` | Command-line parsing for the drivers. |
| [examples](examples/) | | Complete simulations using the library. |
| [legacy](legacy/) | | Old code kept for reference, not compiled. |

## Documentation

Every public function is documented in its header with Doxygen comments (`@brief`, `@param`, `@return`, `@warning`). To generate the HTML documentation in `docs/html`:

```bash
doxygen Doxyfile
```

## Build

The library is compiled into a static library with the compiler and flags of your ITensor installation:

```bash
git clone https://github.com/riccardovalencia/MyTensorNetworkLib.git
cd MyTensorNetworkLib
make LIBRARY_DIR=/path/to/itensor        # -> lib/libmytn.a
make debug LIBRARY_DIR=/path/to/itensor  # -> lib/libmytn-g.a (debug symbols)
```

`LIBRARY_DIR` is the ITensor folder containing `options.mk`; its default is set at the top of the [Makefile](Makefile).

## Usage

Include `mytn.h` (or single headers) with `-I/path/to/MyTensorNetworkLib` and link `lib/libmytn.a` before the ITensor libraries:

```cpp
#include "mytn.h"

int main()
{
    int N = 10;
    auto sites = SpinHalf(N, {"ConserveQNs=", false});
    MPS psi = initial_computational_state(sites, "1100000000");     // |down down up up ...>_z

    vector<double> Delta(N, -1.), Omega(N, 0.1), V(N-1, 1.);
    vector<MyBondGate> gates = gates_rydberg_up_to_VNN(sites, Delta, Omega, V, 0.05);
    for(int step = 0; step < 100; step++)
    {
        for(MyBondGate g : gates) psi = apply_gate(psi, g.gate(), g.jn(), {"Cutoff=", 1E-12, "MaxDim=", 64});
        psi.position(1);
        psi.normalize();
    }
    vector<double> mz = measure_magnetization(&psi, sites, "z");
}
```

The [examples](examples/) contain Makefiles that compile a driver and link it against the library (rebuilding the library when needed).

## Conventions

- Sites are 1-indexed, as in ITensor.
- Spin-1/2: `|0> = |up_z>`, `|1> = |down_z>`. The excitation (Rydberg) projector is `n = (1 - Z)/2 = |down_z><down_z|`.
- X, Y, Z denote Pauli matrices; S^a = sigma^a/2 (ITensor `"Sx"`, `"Sy"`, `"Sz"`).
- Gate lists returned by the TEBD builders implement one second-order Trotter step (forward sweep with dt/2, then the reversed sweep).
- Functions taking `MPS*` modify the state in place (at least its orthogonality center).
- Open systems are simulated on the vectorized (purified) density matrix: the bra occupies sites 1..N in reversed order and the ket sites N+1..2N, so that the dissipative site sits on the central bond.
- Headers contain `using namespace std; using namespace itensor;` for compatibility with existing drivers.

## Migration from the previous layout

Drivers written for the previous version need `#include "mytn.h"` instead of `spin_boson.h`, `bosons.h`, `fermions.h`, `transverse_field_ising_chain.h` or `core.h`, and the following renames:

| Old name | New name |
|---|---|
| `gates_local_lindbland` | `gates_local_lindblad` |
| `gates_nearest_neighbour_local_lindbland` | `gates_nearest_neighbour_local_lindblad` |
| `gates_local_nsites_lindbland` | `gates_local_nsites_lindblad` |
| `TEBD_lindbland_time_evolve` | `TEBD_lindblad_time_evolve` |
| `MPO_lindbland_time_evolve` | `MPO_lindblad_time_evolve` |
| `TEBD_long_range_int_lindbland_time_evolve` | `TEBD_long_range_int_lindblad_time_evolve` |
| `compute_norm_purifed_impurity(_QN)` | `compute_norm_purified_impurity(_QN)` |
| `exctract_reduced_density_matrix` | `extract_reduced_density_matrix` |
| `generaring_function_sim_size` | `generating_function_sim_size` |
| `weigth_coherent_state`, `weigth_squeezed_state` | `weight_coherent_state`, `weight_squeezed_state` |
| `initial_computational_state(sites, vector<int>)` | `initial_computational_state(sites, "0101...")` |
| `initial_state_all_UP/DOWN(sites, &psi, N)` | `psi = initial_computational_state(sites, "00..."/"11...", "x")` |
| `initial_state_DOMAIN_WALL(sites, &psi, N)` | `psi = initial_computational_state(sites, "0..01..1", "x")` |
| `initial_state_vacuum_state(&psi, sites, size)` (no link indices) | `initial_state_vacuum_state_correct_link(&psi, sites, size)` |
| `compute_variance_H_mmGcbQEM(..., s, c)` | `compute_variance_H_mmGcbQEM(..., s, c, symmetry_sector_dir)` |
| `scalar_product_different_n0/_cutoff(..., s, c)` | `scalar_product_different_n0/_cutoff(..., s, c, symmetry_sector_dir)` |
| `search_ground_state_max_bond_chi(_no_v)(..., scaling_bond_dimension)` | `search_ground_state_max_bond_chi(_no_v)(..., scaling_bond_dimension, symmetry_sector_dir)` |

Behaviour fixes with respect to the previous version: `measure_magnetization` and the purified-state measurements
return `<sigma^y>` with the correct sign; `initialize_spin_boson_state` uses the azimuthal angle `phi`;
`H_tight_binding_electrons` adds the on-site fields on every site; `get_data_TEBD` reads `total_time` from `argv[8]`.

## License

[MIT](https://choosealicense.com/licenses/mit/)
