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
| [mps](mps/) | `gates.h` | Gate containers (`TebdGate`, `DissipativeGate`, `OperatorPair`) and `apply_gate`. |
| | `entanglement.h` | Von Neumann entanglement entropy. |
| | `mps_tools.h` | Embedding an MPS into another (`insert_state`), swap gates, density matrices, two-point functions. |
| [dof](dof/) | `spin_half.h` | Spin-1/2 product states (`make_product_state`) and local observables. |
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
    MPS psi = make_product_state(sites, "1100000000");     // |down down up up ...>_z

    vector<double> Delta(N, -1.), Omega(N, 0.1), V(N-1, 1.);
    vector<TebdGate> gates = make_rydberg_gates_nn(sites, Delta, Omega, V, 0.05);
    for(int step = 0; step < 100; step++)
    {
        for(TebdGate g : gates) psi = apply_gate(psi, g.gate(), g.sites(), {"Cutoff=", 1E-12, "MaxDim=", 64});
        psi.position(1);
        psi.normalize();
    }
    vector<double> mz = measure_magnetization(&psi, sites, "z");
}
```

The [examples](examples/) are built all at once by `examples/Makefile` (`cd examples && make`), which links them against the library (rebuilding it when needed); `make run` runs each of them on its sample input.

## Conventions

- Names: files and functions in snake_case (acronyms in lowercase, e.g. `perform_dmrg`), classes in PascalCase.
- Sites are 1-indexed, as in ITensor.
- Spin-1/2: `|0> = |up_z>`, `|1> = |down_z>`. The excitation (Rydberg) projector is `n = (1 - Z)/2 = |down_z><down_z|`.
- X, Y, Z denote Pauli matrices; S^a = sigma^a/2 (ITensor `"Sx"`, `"Sy"`, `"Sz"`).
- Gate lists returned by the TEBD builders implement one second-order Trotter step (forward sweep with dt/2, then the reversed sweep).
- Functions taking `MPS*` modify the state in place (at least its orthogonality center).
- Open systems are simulated on the vectorized (purified) density matrix: the bra occupies sites 1..N in reversed order and the ket sites N+1..2N, so that the dissipative site sits on the central bond.
- Headers contain `using namespace std; using namespace itensor;` for compatibility with existing drivers.

## Migration from the previous layout

Drivers written for previous versions need `#include "mytn.h"` instead of `spin_boson.h`, `bosons.h`, `fermions.h`, `transverse_field_ising_chain.h` or `core.h`, and the following renames:

| Old name | New name |
|---|---|
| `bQEM_dressed_operator` | `make_dressed_operator` |
| `build_file_TEBD` | `build_file_tebd` |
| `build_TEBD_dt_step_H` | `build_tebd_dt_step` |
| `build_TEBD_dt_step_H_open` | `build_tebd_dt_step_open` |
| `build_totalSx` | `make_block_sx_mpo` |
| `compute_norm_purifed_impurity` | `compute_trace_purified` |
| `compute_norm_purifed_impurity_QN` | `compute_trace_purified_qn` |
| `compute_variance_H_mmGcbQEM` | `compute_bosonic_east_model_energy_variance` |
| `exp_H_mmGcbQEM_n0_notfixed` | `make_bosonic_east_model_evolution_mpo` |
| `exctract_reduced_density_matrix` | `compute_reduced_density_matrix` |
| `from_MPS_to_MPDO` | `make_density_matrix_mpo` |
| `from_MPS_to_MPDO_v2` | `make_density_matrix_mpo_fused` |
| `gates_coherent_part_spin_dissipative_NNN_interactions_impurity_model` | `make_spin_impurity_nnn_gates` |
| `gates_local_lindbland` | `make_local_dissipative_gates` |
| `gates_local_nsites_lindbland` | `make_multisite_dissipative_gates` |
| `gates_nearest_neighbour_local_lindbland` | `make_two_site_dissipative_gates` |
| `gates_rydberg_up_to_VNN` | `make_rydberg_gates_nn` |
| `gates_rydberg_up_to_VNNN` | `make_rydberg_gates_nnn` |
| `gates_rydberg_up_to_VNNN_deprecated` | `make_rydberg_gates_nnn_deprecated` |
| `generaring_function_sim_size` | `generating_function_sim_size` |
| `get_data_DMRG` | `get_data_dmrg` |
| `get_data_TEBD` | `get_data_tebd` |
| `H_mmGcbQEM` | `make_bosonic_east_model_mpo` |
| `H_mmGcbQEM_minus` | `make_bosonic_east_model_mpo_minus` |
| `H_mmGcbQEM_n0_not_fixed` | `make_bosonic_east_model_mpo_n0_not_fixed` |
| `H_mmGcbQEM_n0_untouched` | `make_bosonic_east_model_mpo_n0_untouched` |
| `H_mmGcbQEM_onsite` | `make_bosonic_east_model_mpo_onsite` |
| `H_mmGcbQEM_onsite_hopping` | `make_bosonic_east_model_mpo_onsite_hopping` |
| `H_mmGcbQEM_onsite_nonext` | `make_bosonic_east_model_mpo_onsite_nonext` |
| `H_mmGcbQEM_with_drift` | `make_bosonic_east_model_mpo_with_drift` |
| `H_number_conserving_fermions` | `make_single_particle_hamiltonian` |
| `H_number_conserving_fermions_impurity` | `make_single_particle_hamiltonian_impurity` |
| `H_tight_binding_electrons` | `make_tight_binding_mpo` |
| `insert_QN_state` | `insert_qn_state` |
| `MPO_lindbland_time_evolve` | `mpo_lindblad_time_evolve` |
| `perform_DMRG` | `perform_dmrg` |
| `perform_DMRG_meanfield` | `perform_dmrg_meanfield` |
| `perform_DMRG_soft` | `perform_dmrg_soft` |
| `perform_DMRG_variance` | `perform_dmrg_variance` |
| `print_input_DMRG` | `write_dmrg_input` |
| `print_input_DMRG_hopping` | `write_dmrg_input_hopping` |
| `TEBD_lindbland_time_evolve` | `tebd_lindblad_time_evolve` |
| `TEBD_long_range_int_lindbland_time_evolve` | `tebd_long_range_int_lindblad_time_evolve` |
| `weigth_coherent_state` | `coherent_state_amplitude` |
| `weigth_squeezed_state` | `squeezed_state_amplitude` |
| `make_product_state(sites, vector<int>)` | `make_product_state(sites, "0101...")` |
| `initial_state_all_UP/DOWN(sites, &psi, N)` | `psi = make_product_state(sites, "00..."/"11...", "x")` |
| `initial_state_DOMAIN_WALL(sites, &psi, N)` | `psi = make_product_state(sites, "0..01..1", "x")` |
| `initial_state_vacuum_state(&psi, sites, size)` (no link indices) | `set_vacuum_state(&psi, sites, size)` |
| `compute_variance_H_mmGcbQEM(..., s, c)` | `compute_bosonic_east_model_energy_variance(..., s, c, symmetry_sector_dir)` |
| `compute_overlap_different_n0/_cutoff(..., s, c)` | `compute_overlap_different_n0/_cutoff(..., s, c, symmetry_sector_dir)` |
| `load_ground_state_max_bond_dimension(_no_v)(..., scaling_bond_dimension)` | `load_ground_state_max_bond_dimension(_no_v)(..., scaling_bond_dimension, symmetry_sector_dir)` |

Behaviour fixes with respect to the previous version: `measure_magnetization` and the purified-state measurements
return `<sigma^y>` with the correct sign; `make_spin_boson_state` uses the azimuthal angle `phi`;
`make_tight_binding_mpo` adds the on-site fields on every site; `get_data_tebd` reads `total_time` from `argv[8]`;
`make_purified_gates` maps the gates to the correct sites of the bra-ket chain.

## License

[MIT](https://choosealicense.com/licenses/mit/)
