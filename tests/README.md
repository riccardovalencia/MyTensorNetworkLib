# Tests

Integration tests of MyTensorNetworkLib: the tensor-network (TN) programs of [../examples](../examples/)
are run on small inputs and compared with exact diagonalization (ED).

## Run

```bash
cd tests
make test LIBRARY_DIR=/path/to/itensor     # builds the library and the examples, then runs all the tests
make test TEST=test_integration.TimeStepConvergence             # one class
make test TEST=test_integration.TensorNetworkVsExactDiagonalization.test_leaky_cavity   # one test
```

Requirements: the ITensor v3 installation used for the library, Python 3.11+ with numpy, scipy and quimb
(the same as the ED scripts of the examples). The whole suite takes about 45 seconds.
Set `MYTN_KEEP_TEST_DATA=1` to keep the temporary folders with the TN and ED outputs of each test.

## How a test works

1. The input file `inputs/<case>.txt` (same format as the sample inputs of the examples, smaller
   parameters) is copied into a temporary folder, with some parameters possibly overridden (e.g. `dt`).
2. The TN program `examples/<example>/<program>` runs there and writes its results to `data/`.
3. The ED script `examples/<example>/<program>_exact_diagonalization.py` is imported and its
   `main(input_file)` solves the same problem exactly, matches the rows of the TN files (same time
   or theta) and returns the largest |TN - ED| of every observable.
4. The test checks that the largest difference is below a tolerance.

The TN results are not exact: they carry the discretization error of the time step (and, for the
thermal state, of the imaginary-time step). The tolerances are about three times the errors observed
when the tests were written, so that a regression is caught while round-off and DMRG noise are not.

## Tests

`TensorNetworkVsExactDiagonalization`: one test per model and regime.

| Test | Input | What it exercises | max \|TN - ED\| | Tolerance |
|---|---|---|---|---|
| `test_rydberg_chain` | `rydberg_chain.txt` | three-site gates in three layers (`make_rydberg_gates_nnn`) | 1.2e-5 | 5e-5 |
| `test_ising_quench` | `ising_quench.txt` | spin-chain gates (`make_ising_gates`), entanglement entropy, full counting statistics | 1.8e-5 | 5e-5 |
| `test_cavity_dicke` | `cavity_dicke.txt` | boson-spin swap gates, Dicke coupling (`make_light_matter_gates`) | 2.0e-3 | 5e-3 |
| `test_cavity_tavis_cummings` | `cavity_tavis.txt` | same, Tavis-Cummings coupling | 3.2e-5 | 1e-4 |
| `test_ising_thermal` | `ising_thermal.txt` | imaginary-time evolution of a purified MPO, generating functions of a mixed state | 1.4e-4 | 5e-4 |
| `test_leaky_cavity` | `leaky_cavity.txt` | Lindblad dynamics: `make_purified_gates`, local dissipators, purified-state observables | 4.5e-4 | 1.5e-3 |
| `test_impurity_integrable` | `impurity_integrable.txt` | impurity gates and dissipator on the central bond, regression theorem; exact free-fermion solution | 1.5e-3 | 5e-3 |
| `test_impurity_xx_coupling` | `impurity_xx.txt` | same with an XX coupling (dense Lindblad solution) | 2.9e-3 | 1e-2 |
| `test_impurity_next_nearest_neighbour` | `impurity_nnn.txt` | three-site impurity gates (`make_spin_impurity_nnn_gates`) on N = 5 | 3.1e-3 | 1e-2 |

`TimeStepConvergence`: the TN error is a discretization error with the expected order. Each test
runs the same input with time steps dt and dt/2 and checks the ratio of the errors (at least 80% of
2^order).

| Test | Input, dt | Expected order | Ratio |
|---|---|---|---|
| `test_second_order_trotter` | `ising_quench.txt`, 0.02 | 2 (symmetric sweep, `make_symmetric_sweep`) | 4.00 |
| `test_first_order_dissipative_step` | `impurity_nnn.txt`, 0.02 | 1 (first-order dissipative gates) | 2.01 |
| `test_first_order_local_long_range_splitting` | `cavity_dicke.txt`, 0.005 | 1 (local and photon-matter sweeps in sequence) | 2.01 |

`GroundStateSearch`: the DMRG driver `find_ground_state` ([../ground_state/dmrg.h](../ground_state/dmrg.h)),
through the example `spin_chain_ground_state`, on short-range spin-1/2 chains of 12-18 sites in the
regimes where DMRG can get stuck. ED is a Lanczos ground state (in the sector S^z = 0 when S^z is
conserved). Compared: the energy and the profile `<Z_j>`, `<Z_j Z_{j+1}>`, entanglement entropy S(j, j+1).

| Test | Input | What it exercises | Energy error | Profile error | Tolerances (E, profile) |
|---|---|---|---|---|---|
| `test_frustrated_j1_j2_chain` | `ground_state_j1j2.txt` | J1-J2 Heisenberg chain (J2 = 0.3) in a field, no conserved quantities, random initial states, 2 restarts | 8.7e-11 | 1.2e-9 | 1e-8, 1e-6 |
| `test_conserved_magnetization` | `ground_state_xxz_sz.txt` | XXZ chain with next-nearest-neighbour couplings, conserved S^z: start from the Neel state (`InitState` overload) | 9.4e-11 | 1.8e-9 | 1e-8, 1e-6 |
| `test_quasi_degenerate_ground_state` | `ground_state_quasi_degenerate.txt` | Ising chain near its transition, gap 8.5e-4: with ITensor's default 2 Davidson iterations per bond DMRG stops in a mixture of the two lowest states (energy error 4.5e-4) | 1.8e-10 | 6e-4 (`<Z_j>`) | 1e-8, 2e-3 |
| `test_restarts_escape_metastable_state` | `ground_state_metastable.txt` | ferromagnetic Ising chain magnetized against a weak field: ~1/3 of the single runs end in this metastable state (variance ~1e-11, so only the energy reveals it); the test finds such a seed and checks that 8 restarts from it return the ground state | 7.2e-12 | 6.0e-11 | 1e-8, 1e-6 |
| `test_bond_dimension_convergence` | `ground_state_critical_ising.txt` | critical transverse-field Ising chain with `max_dim` = 4, 8, 16: the energy error must drop by more than 100 at each doubling | 1.5e-4, 1.9e-7, 4.5e-11 | | ratios > 100, last < 1e-8 |

The DMRG tolerances are not set by the observed errors (~1e-10, round-off and cutoff) but by the
errors of the failures they are meant to catch, which are orders of magnitude larger. Replacing the
10 Davidson iterations with ITensor's default 2, or returning the first restart instead of the
lowest one, makes `test_quasi_degenerate_ground_state` and `test_restarts_escape_metastable_state` fail.

Notes on the inputs:
- The impurity tests with `Jxx` or `Jzzz` use `hx = 1.5`: the ground state is then gapped, so DMRG
  finds a unique initial state (for `hx = 0.5` the two lowest states are almost degenerate and the TN
  and ED initial states may differ).
- `t_measure` and `t_corr` must be multiples of `dt` also when `dt` is overridden.

## Adding a test

1. Add an input file to `inputs/` (small system: the ED cost grows exponentially).
2. If the example has no ED script yet, write `<program>_exact_diagonalization.py` next to the program,
   with a `main(input_file)` returning the dictionary of `save_and_compare` results (see
   [../examples/exact_diagonalization_tools.py](../examples/exact_diagonalization_tools.py)).
3. Add a method to `TensorNetworkVsExactDiagonalization` (or `GroundStateSearch`) in [test_integration.py](test_integration.py):

   ```python
   def test_my_model(self):
       """What the test exercises."""
       comparisons = run_and_compare('my_example', 'my_program', 'my_model.txt')
       self.assert_matches_exact_diagonalization(comparisons, tolerance)
   ```

   and choose the tolerance from the error observed with `make test` (about three times larger).
