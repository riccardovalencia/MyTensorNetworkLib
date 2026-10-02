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
| [mps](mps/) | `gates.h` | Gate containers (`TebdGate`, `DissipativeGate`, `OperatorPair`), `apply_gate(s)`, symmetric Trotter sweeps (`make_symmetric_sweep`), weights of terms shared by overlapping gates (`count_gates_containing`). |
| | `entanglement.h` | Von Neumann entanglement entropy. |
| | `mps_tools.h` | Embedding an MPS into another (`insert_state`), swaps, product-state tensors (`set_site_tensor`), local and two-point expectation values, density matrices. |
| [dof](dof/) | `spin_half.h` | Spin-1/2 product states (`make_product_state`), Pauli operators, magnetization and correlations. |
| | `boson.h` | Bosonic Fock, coherent, squeezed and cat states; occupations, Fock-space probabilities, squeezing. |
| | `fermion.h` | Product states of spinful fermions (`InitState` and MPS). |
| | `spin_boson.h` | Boson + spins site sets (cavity QED) and their purified (bra-ket) version; initial states. |
| [models](models/) | `spin_chain.h` | Short-range spin chains: MPO (or AutoMPO terms) with nearest- and next-nearest-neighbour couplings (`make_spin_chain_mpo`, `make_spin_chain_terms`), bond terms and TEBD gates, Ising chain in longitudinal and transverse fields. |
| | `rydberg.h` | Rydberg arrays (up to next-nearest-neighbour interactions), PXP model. |
| | `light_matter.h` | Dicke and Tavis-Cummings models (one boson coupled to all spins). |
| | `impurity.h` | Impurity models on the purified density matrix (Kondo-like, dissipative spin impurity). |
| | `tight_binding.h` | Free fermions: MPOs, single-particle matrices, bond terms and TEBD gates. |
| | `bosonic_east_model.h`, `bosonic_east_model_states.h` | Bosonic quantum east model: Hamiltonians, TEBD gates, super-bosonic states, dressed operators. |
| [dynamics](dynamics/) | `time_evolution.h` | One TEBD step (`tebd_step`) for pure states and purified density matrices; measurement intervals. |
| | `heisenberg.h` | Operators in the Heisenberg picture: one step O -> U^dag O U with TEBD gates or a propagator MPO (`heisenberg_step`), adjoint MPOs. |
| | `lindblad.h` | Dissipative gates for Lindblad dynamics on the vectorized density matrix; gates on the bra-ket chain. |
| | `purified_state.h` | Trace and observables of a density matrix stored as a purified MPS. |
| | `adiabatic.h` | Adiabatic ramps (bosonic quantum east model). |
| [ground_state](ground_state/) | `dmrg.h` | Ground states with DMRG (`find_ground_state`): bond-dimension ramp, noise, convergence check, random restarts; Fermi sea. |
| [analysis](analysis/) | `full_counting_statistics.h` | Generating function and cumulants of the subsystem magnetization. |
| [io](io/) | `output.h` | Run folders with a copy of the input (`make_run_directory`), output files: site profiles, tables, generating functions, ground-state energies and convergence (`write_ground_state`). |
| [examples](examples/) | | Complete simulations, each with an exact-diagonalization check (see [examples/README.md](examples/README.md)). |
| [tests](tests/) | | Integration tests: tensor networks against exact diagonalization (see [tests/README.md](tests/README.md)). |
| [legacy](legacy/) | | Old code kept for reference, not compiled. |

## Documentation

Every public function is documented in its header with Doxygen comments (`@brief`, `@param`, `@return`, `@throws`, `@warning`). To generate the HTML documentation in `docs/html`:

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
    vector<TebdGate> gates = make_rydberg_gates_nn(sites, Delta, Omega, V, 0.05);   // one step dt = 0.05
    for(int step = 0; step < 100; step++)
    {
        psi = tebd_step(psi, gates, {"Cutoff=", 1E-12, "MaxDim=", 64});
        psi.position(1);
        psi.normalize();
    }
    vector<double> mz = measure_magnetization(&psi, sites, "z");
}
```

Ground states come from `find_ground_state`, which needs only the MPO; its parameters have defaults
(see [ground_state/dmrg.h](ground_state/dmrg.h)):

```cpp
auto sites = SpinHalf(N, {"ConserveQNs=", false});
MPO H = make_spin_chain_mpo(sites, {1., 1., 1.}, {0.3, 0.3, 0.3}, {0., 0., 0.1});   // J1-J2 chain in a field hz

DmrgParameters parameters;          // max_dim 200, cutoff 1E-12, noise 1E-6, tolerance 1E-10, 30 sweeps
parameters.number_restarts = 4;     // 4 runs from random states, the lowest energy is kept
GroundState ground_state = find_ground_state(H, parameters);
// ground_state.psi, .energy, .variance, .converged, .restart_energies; write_ground_state(dir, ground_state) saves them
```

With conserved quantum numbers, pass the product state that fixes the sector instead:
`find_ground_state(H, InitState(...), parameters)`. In the programs, `read_dmrg_parameters(input)` reads
the parameters from the input file (keys `dmrg_max_dim`, `dmrg_cutoff`, `dmrg_noise`, `dmrg_tolerance`,
`dmrg_max_sweeps`, `dmrg_restarts`, `dmrg_seed`).

The [examples](examples/) are built all at once by `examples/Makefile` (`cd examples && make`), which links them against the library (rebuilding it when needed); `make run` runs each of them on its sample input and `make compare` compares them with exact diagonalization.

## Tests

[tests/](tests/) contains integration tests: each example is run on a small input and compared with exact diagonalization, the time-step convergence order of the TEBD schemes is checked, the DMRG ground-state search is tested where DMRG can get stuck (frustration, conserved quantities, quasi-degenerate ground states, metastable states, truncation), and operators evolved in the Heisenberg picture are compared with states evolved in the Schroedinger picture.

```bash
cd tests && make test LIBRARY_DIR=/path/to/itensor     # ~70 s; needs numpy, scipy and quimb
```

## Conventions

- Names: files, functions and variables in snake_case, classes in PascalCase; physics symbols keep their usual
  spelling (`N`, `J`, `hx`, `Omega`). Functions start with a verb: `make_*` returns a new object, `set_*` overwrites
  a state, `apply_*`/`evolve_*` act on a state, `measure_*` returns expectation values, `compute_*` other derived
  quantities, `write_*`/`load_*` handle files. Pure mathematical helpers are nouns (`factorial`, `coherent_state_amplitude`).
- Types are written explicitly (no `auto`): lambdas are stored in `std::function`, the results of `svd` are
  unpacked with `std::tie`.
- Errors (invalid arguments, missing files) throw `ITError`; the library never calls `exit`.
- Sites are 1-indexed, as in ITensor.
- Spin-1/2: `|0> = |up_z>`, `|1> = |down_z>`. The excitation (Rydberg) projector is `n = (1 - Z)/2 = |down_z><down_z|`.
- X, Y, Z denote Pauli matrices; S^a = sigma^a/2 (ITensor `"Sx"`, `"Sy"`, `"Sz"`).
- Gate lists returned by the TEBD builders implement one second-order Trotter step: a sweep with dt/2 followed by
  the reversed sweep (`make_symmetric_sweep`). On-site terms shared by several gates are divided among them
  (`count_gates_containing`).
- Functions taking `MPS*` modify the state in place (at least its orthogonality center).
- Open systems are simulated on the vectorized (purified) density matrix: the bra occupies sites 1..N in reversed order and the ket sites N+1..2N, so that the dissipative site sits on the central bond.
- Headers contain `using namespace std; using namespace itensor;` for compatibility with existing drivers.

## Migration from the previous versions

Drivers written for previous versions need `#include "mytn.h"` instead of `spin_boson.h`, `bosons.h`, `fermions.h`,
`transverse_field_ising_chain.h` or `core.h`, and the following renames. Functions that used to fill a vector passed
by reference (`vector<BondGate>& gates`, `ITensor* hterm`) now return it.

The bosonic east model functions take the couplings that enter the Hamiltonian, J = e^{-s} and U = 1 - 2c,
instead of s and c: pass `exp(-s), 1 - 2*c` where you passed `s, c` (and `1 - 2*c` where you passed `c` next to J).
This concerns the MPOs (`make_bosonic_east_model_mpo*`, except `_onsite` and `_onsite_nonext`), the bond terms and
gates, `make_bosonic_east_model_evolution_mpo`, `evolve_adiabatic_*_ramp`, `make_super_bosonic_*_state` and
`make_dressed_operator`. The functions that read or name data files (`compute_bosonic_east_model_energy_variance`,
`compute_overlap_*`) still take s and c.

| Old name | New name |
|---|---|
| `adiabatic_transformation_linear_protocol` / `_tanh_protocol` | `evolve_adiabatic_linear_ramp` / `evolve_adiabatic_tanh_ramp` |
| `bQEM_dressed_operator` | `make_dressed_operator` |
| `build_single_step` (bosonic east model) | `make_bosonic_east_model_bond_hamiltonian` |
| `build_single_step` (Ising chain) | `make_ising_gates` (whole Trotter step) |
| `build_single_step_n0` | `make_bosonic_east_model_bond_hamiltonian_n0` |
| `build_single_step_jumps` | `make_bosonic_east_model_bond_hamiltonian_dephasing` |
| `build_TEBD_dt_step_H` | `make_bosonic_east_model_gates` |
| `build_totalSx` | `make_block_sx_mpo` |
| `coherent_state_site_j` / `coherent_state_all_sites` | `set_coherent_state_on_site` / `set_coherent_state_all_sites` |
| `compute_norm_purifed_impurity`, `compute_norm_purifed_impurity_QN` | `compute_trace_purified` |
| `compute_potential` | `compute_power_law_couplings` |
| `compute_two_point` | `measure_two_point_function` |
| `compute_variance_H_mmGcbQEM(..., s, c)` | `compute_bosonic_east_model_energy_variance(..., s, c, symmetry_sector_dir)` |
| `custom_spin_boson` / `custom_spin_boson_doubling` | `make_spin_boson_sites` / `make_purified_spin_boson_sites` |
| `doubling_space_gates` | `make_purified_gates` |
| `entanglement_entropy` | `compute_entanglement_entropy` (log2; natural log with `natural_log = true`) |
| `exctract_reduced_density_matrix` | `compute_reduced_density_matrix` |
| `exp_H_mmGcbQEM_n0_notfixed` | `make_bosonic_east_model_evolution_mpo` (without `n0`) |
| `expectation_value_sigma_x(_square, _n)`, `expectation_value_n_sigma_x` | `measure_sigma_x(_squared, _n)`, `measure_n_sigma_x` |
| `fermi_sea_electrons(H, sites, Nup, Ndn, sweeps, min_varH)` | `find_fermi_sea(H, sites, Nup, Ndn, parameters)` (returns a `GroundState`: check its `variance`) |
| `from_MPS_to_MPDO`, `from_MPS_to_MPDO_v2` | `make_density_matrix_mpo` |
| `gates_coherent_part_spin_dissipative_impurity_model` | `make_spin_impurity_gates` (without `Lj`, `gamma`) |
| `gates_coherent_part_spin_dissipative_NNN_interactions_impurity_model` | `make_spin_impurity_nnn_gates` |
| `gates_coherent_unfolded_kondo_impurity_model` | `make_kondo_impurity_gates` |
| `gates_dissipative_impurity` / `_high_pade` | `make_impurity_dissipative_gates` / `make_impurity_dissipative_gates_pade` |
| `gates_free_spinful_fermions` | `make_free_fermion_gates` |
| `gates_local_lindbland` | `make_local_dissipative_gates` |
| `gates_local_nsites_lindbland` | `make_multisite_dissipative_gates` |
| `gates_nearest_neighbour_local_lindbland` | `make_two_site_dissipative_gates` |
| `gates_photon_matter` | `make_light_matter_gates` |
| `gates_pxp` / `mpo_pxp` | `make_pxp_gates` / `make_pxp_mpo` |
| `gates_rydberg_up_to_VNN` | `make_rydberg_gates_nn` |
| `gates_rydberg_up_to_VNNN`, `gates_rydberg_up_to_VNNN_deprecated` | `make_rydberg_gates_nnn` |
| `gates_spin_eff_model_BondGate` | `make_spin_chain_effective_gates` |
| `gates_spin_local_field` | `make_local_field_gates` |
| `gates_spin_model`, `gates_spin_model_BondGate` | `make_spin_chain_gates` |
| `gates_tavis_cummings` | `make_light_matter_gates(..., "tavis")` |
| `generaring_function_sim_size`, `measure_generating_function` | `compute_generating_function` |
| `H_mmGcbQEM(_minus, _n0_not_fixed, _n0_untouched, _onsite, _onsite_hopping, _onsite_nonext, _with_drift)` | `make_bosonic_east_model_mpo(...)` with the same suffix |
| `H_number_conserving_fermions(_impurity)` | `make_single_particle_hamiltonian(_impurity)` |
| `H_tight_binding_electrons` | `make_tight_binding_mpo` |
| `initial_computational_electron_state` | `make_electron_product_state` |
| `initial_computational_state`, `initial_state_all_UP/DOWN`, `initial_state_DOMAIN_WALL` | `make_product_state(sites, "0101...", basis)` |
| `initial_state_vacuum_state(_correct_link)` / `initial_state_all_one_state(_correct_link)` | `set_vacuum_state` / `set_unit_filling_state` |
| `initial_state_n0_excitation(_pinned)` | `set_fock_excitation` |
| `initial_state_cat_state_site_j` / `squeezed_state_site_j` / `kink_state` / `put_occupation` | `set_cat_state_on_site` / `set_squeezed_state_on_site` / `set_kink_state` / `set_site_occupation` |
| `initialize_spin_boson_state` | `make_spin_boson_state` |
| `insert_QN_state` | `insert_state` |
| `max_projector_at_cutoff` | `compute_max_cutoff_probability` |
| `measure_correlation_impurity_first_site` / `measure_local_obs_impurity_first_site` / `measure_magnetization_impurity_first_site` | `measure_correlation_purified` / `measure_local_operator_purified` / `measure_magnetization_purified` |
| `measure_correlations` / `measure_kink` / `measure_mx_mz` | `measure_density_correlations` / `measure_kink_number` / `print_magnetization` |
| `measure_covariance_matrix_number_operator` | `measure_number_covariance` |
| `measure_delta_x` / `measure_delta_p` / `measure_imbalance` | `measure_variance_x` / `measure_variance_p` / `compute_imbalance` |
| `measure_projector_all_sites(..., occupation_number)` | `measure_fock_probabilities(...)` |
| `measure_square_occupation_number` | `measure_occupation_number_squared` |
| `measuring_moments` | `compute_cumulants` |
| `MPO_lindbland_time_evolve`, `TEBD_lindbland_time_evolve`, `TEBD_long_range_int_lindbland_time_evolve` | `tebd_step` (one step; the loop and the output are in the program) |
| `MyBondGate` / `MyBondGateDiss` / `MyTrainITensor` | `TebdGate` / `DissipativeGate` / `OperatorPair` |
| `perform_DMRG`, `perform_DMRG_soft` | `find_ground_state(H, parameters)` (no output files: write `energy` and `variance` in the program) |
| `perform_DMRG_meanfield` | `find_ground_state(H, parameters)` with `max_dim = 1` (and several restarts) |
| `perform_DMRG_variance` | `find_ground_state` of the MPO (H - E)^2 (`nmultMPO`) |
| `printing_generating_function` | `write_generating_function` |
| `scalar_product_different_n0` / `_cutoff` | `compute_overlap_different_n0` / `compute_overlap_different_cutoffs` (with `symmetry_sector_dir`) |
| `super_bosonic_state` / `_coherent_state` / `_squeezed_state` | `make_super_bosonic_state` / `_coherent_state` / `_squeezed_state` |
| `swap_gate` | `swap_sites` |
| `theta_step` | `make_theta_grid` |
| `weigth_coherent_state` / `weigth_squeezed_state` | `coherent_state_amplitude` / `squeezed_state_amplitude` |

Removed without replacement: the drivers' file and argument helpers (`build_file_*`, `get_data*`, `print_info`,
`print_input`, `print_matrix`, `print_occupation_number*`, `print_projector_fockspace*`; the examples read InputGroup
files and write with `write_site_values` / `write_site_table`), `build_single_step_jumps_v2` and
`build_TEBD_dt_step_H_open` (deprecated since 2022; use `make_bosonic_east_model_gates(..., "open", gamma)`),
`gates_coherent_unfolded_kondo_impurity_model_energy_basis` (it needs gates on non-consecutive sites), and
`initialize_excited_state` (initial guess of `perform_DMRG_variance`; `find_ground_state` starts from random states).
The DMRG input/output tied to the bosonic east model runs (`print_input_DMRG(_hopping)`, `search_ground_state_max_bond_chi(_no_v)`,
`search_state_adiabatic_coherent`, which read and wrote one fixed folder layout) are replaced by the generic
`write_ground_state` (energy, variance, convergence and the energy of every restart of `find_ground_state`); states
can be stored and read back with ITensor's `writeToFile` / `readFromFile`.

Behaviour fixes with respect to the previous versions: `measure_magnetization` and the purified-state measurements
return `<sigma^y>` with the correct sign, and `measure_local_operator_purified` no longer transposes the operator
when measuring all sites; `make_spin_boson_state` uses the azimuthal angle `phi`; `make_tight_binding_mpo` adds the
on-site fields on every site; `make_purified_gates` maps the gates to the correct sites of the bra-ket chain; the
two-site dissipators (`make_two_site_dissipative_gates`, `make_multisite_dissipative_gates`) use the full identity on
two sites; the three-site gates (`make_rydberg_gates_nnn`, `make_spin_impurity_nnn_gates`) split the on-site terms
correctly for N <= 5; `compute_overlap_different_n0/_cutoffs` return |overlap|^2; `gates_tavis_cummings` exchanged a and a^dag;
`make_pauli_operator`, `make_magnetization_operator` and `make_identity_operator` also work on sites with conserved
quantum numbers (the input index is `dag`-ed); `compute_bosonic_east_model_energy_variance` also works on complex
states (e.g. after a time evolution).

## License

[MIT](https://choosealicense.com/licenses/mit/)
