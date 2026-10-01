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
(the same as the ED scripts of the examples). The whole suite takes about 15 seconds.
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
3. Add a method to `TensorNetworkVsExactDiagonalization` in [test_integration.py](test_integration.py):

   ```python
   def test_my_model(self):
       """What the test exercises."""
       comparisons = run_and_compare('my_example', 'my_program', 'my_model.txt')
       self.assert_matches_exact_diagonalization(comparisons, tolerance)
   ```

   and choose the tolerance from the error observed with `make test` (about three times larger).
