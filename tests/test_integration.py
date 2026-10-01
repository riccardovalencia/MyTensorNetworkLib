"""Integration tests: tensor-network (TN) examples against exact diagonalization (ED).

Each test runs a TN program of examples/ on a small input of tests/inputs/ in a temporary folder,
then the exact-diagonalization script of the same example (its main(input_file)), and checks that
the largest |TN - ED| over all the compared observables is below a tolerance. For the dynamics the
tolerances are set by the expected discretization errors (Trotter step, imaginary-time step,
first-order dissipative step), not by the truncation, which is negligible for these sizes: each
tolerance is about three times the error observed when the tests were written (listed in README.md),
so that a regression is caught while round-off and DMRG noise are not. The ground-state tests
(GroundStateSearch) check the DMRG driver in the regimes where DMRG can get stuck.

Run from tests/ with `make test`, or `python3 -m unittest -v test_integration` once the examples
are built. Set MYTN_KEEP_TEST_DATA=1 to keep the temporary folders (their paths are printed).
"""
import contextlib
import importlib.util
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLES_DIR = os.path.normpath(os.path.join(TESTS_DIR, '..', 'examples'))
sys.path.insert(0, EXAMPLES_DIR)
from exact_diagonalization_tools import read_input  # noqa: E402

KEEP_DATA = bool(os.environ.get('MYTN_KEEP_TEST_DATA'))
TIMEOUT = 600   # seconds per TN program


def write_input(path, params):
    """Write an ITensor InputGroup file "input { key = value ... }"."""
    with open(path, 'w') as f:
        f.write('input\n{\n' + ''.join(f'{key} = {value}\n' for key, value in params.items()) + '}\n')


def load_exact_diagonalization(example, program):
    """Module examples/<example>/<program>_exact_diagonalization.py."""
    path = os.path.join(EXAMPLES_DIR, example, f'{program}_exact_diagonalization.py')
    spec = importlib.util.spec_from_file_location(f'{program}_exact_diagonalization', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_and_compare(example, program, input_name, **overrides):
    """Run examples/<example>/<program> and its ED script on tests/inputs/<input_name> (with the
    parameters in overrides replaced) and return {output: {column: max |TN - ED|}}."""
    params = read_input(os.path.join(TESTS_DIR, 'inputs', input_name))
    params.update(overrides)
    work = tempfile.mkdtemp(prefix=f'mytn_test_{program}_')
    try:
        write_input(os.path.join(work, 'input.txt'), params)
        binary = os.path.join(EXAMPLES_DIR, example, program)
        if not os.path.exists(binary):
            raise unittest.SkipTest(f'{binary} not built: run `make` in examples/ (or `make test` in tests/)')
        run = subprocess.run([binary, 'input.txt'], cwd=work, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                             text=True, timeout=TIMEOUT)
        if run.returncode != 0:
            raise AssertionError(f'{program} failed:\n{run.stderr[-2000:]}')
        exact_diagonalization = load_exact_diagonalization(example, program)
        with contextlib.chdir(work), contextlib.redirect_stdout(io.StringIO()):
            return exact_diagonalization.main('input.txt')
    finally:
        if KEEP_DATA:
            print(f'\n  data kept in {work}', file=sys.stderr)
        else:
            shutil.rmtree(work, ignore_errors=True)


def max_difference(comparisons):
    """Largest |TN - ED| over all outputs and columns, with its label 'output:column'."""
    return max((diff, f'{output}:{column}') for output, columns in comparisons.items() for column, diff in columns.items())


def convergence_ratio(example, program, input_name, dt):
    """max |TN - ED| with time step dt divided by the one with dt/2 (~2^order for a method of that order)."""
    coarse = max_difference(run_and_compare(example, program, input_name, dt=dt))[0]
    fine = max_difference(run_and_compare(example, program, input_name, dt=dt / 2))[0]
    return coarse / fine


class TensorNetworkVsExactDiagonalization(unittest.TestCase):
    """One test per example (and per regime of the model): TN result within tolerance of ED."""

    def assert_matches_exact_diagonalization(self, comparisons, tolerance):
        for output, columns in comparisons.items():
            self.assertIsNotNone(columns, f'{output}: no TN result to compare with')
        diff, label = max_difference(comparisons)
        print(f'\n  max |TN - ED| = {diff:.2e} ({label}), tolerance {tolerance:.1e}', file=sys.stderr, end=' ')
        self.assertLess(diff, tolerance, f'max |TN - ED| = {diff:.2e} at {label}')

    # ---- closed systems, pure states -------------------------------------------------------

    def test_rydberg_chain(self):
        """Three-site TEBD gates in three layers (make_rydberg_gates_nnn, count_gates_containing)."""
        comparisons = run_and_compare('rydberg_chain_tebd', 'rydberg_chain_tebd', 'rydberg_chain.txt')
        self.assert_matches_exact_diagonalization(comparisons, 5e-5)

    def test_ising_quench(self):
        """Spin-chain gates (make_ising_gates), entanglement entropy and full counting statistics."""
        comparisons = run_and_compare('ising_full_counting_statistics', 'ising_quench', 'ising_quench.txt')
        self.assert_matches_exact_diagonalization(comparisons, 5e-5)

    def test_cavity_dicke(self):
        """Spin-boson chain with all-to-all Dicke coupling via swap gates (make_light_matter_gates)."""
        comparisons = run_and_compare('cavity_unitary', 'cavity_unitary', 'cavity_dicke.txt')
        self.assert_matches_exact_diagonalization(comparisons, 5e-3)

    def test_cavity_tavis_cummings(self):
        """As test_cavity_dicke with the Tavis-Cummings coupling g (sigma^+ a + sigma^- a^dag)."""
        comparisons = run_and_compare('cavity_unitary', 'cavity_unitary', 'cavity_tavis.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-4)

    # ---- imaginary time ----------------------------------------------------------------------

    def test_ising_thermal(self):
        """Thermal state by imaginary-time evolution of the purified MPO, generating functions."""
        comparisons = run_and_compare('ising_full_counting_statistics', 'ising_thermal', 'ising_thermal.txt')
        self.assert_matches_exact_diagonalization(comparisons, 5e-4)

    # ---- open systems, purified density matrices ---------------------------------------------

    def test_leaky_cavity(self):
        """Lindblad dynamics: purified gates (make_purified_gates), local dissipators, folded traces."""
        comparisons = run_and_compare('leaky_cavity', 'leaky_cavity', 'leaky_cavity.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1.5e-3)

    def test_impurity_integrable(self):
        """Dephasing impurity on the central bond (make_spin_impurity_gates, make_impurity_dissipative_gates),
        compared with the exact free-fermion solution; autocorrelation via the regression theorem."""
        comparisons = run_and_compare('dissipative_impurity_ising', 'impurity_dynamics', 'impurity_integrable.txt')
        self.assert_matches_exact_diagonalization(comparisons, 5e-3)

    def test_impurity_xx_coupling(self):
        """As test_impurity_integrable with an XX coupling (non-integrable, dense Lindblad solution)."""
        comparisons = run_and_compare('dissipative_impurity_ising', 'impurity_dynamics', 'impurity_xx.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-2)

    def test_impurity_next_nearest_neighbour(self):
        """Three-site impurity gates (make_spin_impurity_nnn_gates) on N = 5 sites, where every gate
        touches an edge of its half chain. hx = 1.5 keeps the ground state gapped (unique DMRG state)."""
        comparisons = run_and_compare('dissipative_impurity_ising', 'impurity_dynamics', 'impurity_nnn.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-2)


class GroundStateSearch(unittest.TestCase):
    """find_ground_state (ground_state/dmrg.h) on short-range spin chains of 12-18 sites, through the
    example spin_chain_ground_state, against Lanczos: energy and profile j, <Z_j>, <Z_j Z_{j+1}>,
    entanglement entropy S(j, j+1). The DMRG errors are truncation errors (cutoff 1E-12), ~1e-9."""

    example = ('spin_chain_ground_state', 'spin_chain_ground_state')

    def ground_state(self, input_name, **overrides):
        return run_and_compare(*self.example, input_name, **overrides)

    def energy_error(self, comparisons):
        return comparisons['energy']['E']

    def assert_matches_exact_diagonalization(self, comparisons, energy_tolerance, profile_tolerance):
        self.assertIsNotNone(comparisons['energy'], 'no TN energy to compare with')
        self.assertIsNotNone(comparisons['profile'], 'no TN profile to compare with')
        profile_error, label = max((diff, column) for column, diff in comparisons['profile'].items())
        print(f'\n  |E_TN - E_ED| = {self.energy_error(comparisons):.2e} (tolerance {energy_tolerance:.0e}),'
              f' max profile error {profile_error:.2e} ({label}, tolerance {profile_tolerance:.0e})', file=sys.stderr, end=' ')
        self.assertLess(self.energy_error(comparisons), energy_tolerance)
        self.assertLess(profile_error, profile_tolerance, f'profile error {profile_error:.2e} at {label}')

    def test_frustrated_j1_j2_chain(self):
        """Heisenberg chain with next-nearest-neighbour J2 = 0.3 (frustrated, dimerized phase) and a
        field hz, without conserved quantities: random initial states, 2 restarts."""
        comparisons = self.ground_state('ground_state_j1j2.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-8, 1e-6)

    def test_conserved_magnetization(self):
        """XXZ chain with next-nearest-neighbour couplings (gapless) with conserved S^z: DMRG starts from
        the Neel product state (InitState overload), the noise lets the bond dimension grow."""
        comparisons = self.ground_state('ground_state_xxz_sz.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-8, 1e-6)

    def test_quasi_degenerate_ground_state(self):
        """Ising chain close to its transition, gap 8.5e-4 above the (Z2-symmetric) ground state. With
        ITensor's 2 Davidson iterations per bond DMRG stops in a mixture of the two lowest states
        (energy error ~1e-4, <Z_j> ~ 1); the admixture left here is ~1e-4 in <Z_j>."""
        comparisons = self.ground_state('ground_state_quasi_degenerate.txt')
        self.assert_matches_exact_diagonalization(comparisons, 1e-8, 2e-3)

    def test_restarts_escape_metastable_state(self):
        """Ferromagnetic Ising chain, weak transverse field hx = 0.1, longitudinal field hz = 0.002: the
        state magnetized against hz is a local minimum where a single DMRG run gets stuck (~1/3 of the
        seeds, with variance ~1e-11). Find a seed whose first run gets stuck, then check that 8 restarts
        from the same seed (the first run is the same) return the ground state."""
        stuck_seed = None
        for seed in range(1, 21):
            if self.energy_error(self.ground_state('ground_state_metastable.txt', dmrg_seed=seed)) > 1e-2:
                stuck_seed = seed
                break
        if stuck_seed is None:
            self.skipTest('no single run stuck in the metastable state for seeds 1-20')
        comparisons = self.ground_state('ground_state_metastable.txt', dmrg_seed=stuck_seed, dmrg_restarts=8)
        print(f'\n  seed {stuck_seed}: single run stuck, 8 restarts:', file=sys.stderr, end=' ')
        self.assert_matches_exact_diagonalization(comparisons, 1e-8, 1e-6)

    def test_bond_dimension_convergence(self):
        """Critical transverse-field Ising chain: the energy error is a truncation error that drops by
        orders of magnitude when the maximal bond dimension doubles (4 -> 8 -> 16), down to round-off."""
        errors = [self.energy_error(self.ground_state('ground_state_critical_ising.txt', dmrg_max_dim=D)) for D in (4, 8, 16)]
        print(f'\n  |E_TN - E_ED| for max_dim 4, 8, 16: {", ".join(f"{e:.1e}" for e in errors)}', file=sys.stderr, end=' ')
        self.assertGreater(errors[0], 100 * errors[1])
        self.assertGreater(errors[1], 100 * errors[2])
        self.assertLess(errors[2], 1e-8)


class TimeStepConvergence(unittest.TestCase):
    """The TN error is a discretization error: it decreases with the order of the time step."""

    def assert_order(self, ratio, order):
        print(f'\n  error ratio dt / (dt/2) = {ratio:.2f}, expected ~{2 ** order}', file=sys.stderr, end=' ')
        self.assertGreater(ratio, 0.8 * 2 ** order)

    def test_second_order_trotter(self):
        """make_symmetric_sweep: one symmetric sweep per step is a second-order Trotter step."""
        self.assert_order(convergence_ratio('ising_full_counting_statistics', 'ising_quench', 'ising_quench.txt', 0.02), 2)

    def test_first_order_dissipative_step(self):
        """tebd_step with dissipative gates: rho + dt D rho is first order in dt."""
        self.assert_order(convergence_ratio('dissipative_impurity_ising', 'impurity_dynamics', 'impurity_nnn.txt', 0.02), 1)

    def test_first_order_local_long_range_splitting(self):
        """cavity_unitary applies the local and the photon-matter sweeps one after the other (first order)."""
        self.assert_order(convergence_ratio('cavity_unitary', 'cavity_unitary', 'cavity_dicke.txt', 0.005), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
