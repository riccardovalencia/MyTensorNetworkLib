# Examples

Complete simulations built on MyTensorNetworkLib. Each folder contains one or more programs
`<name>.cc`, a sample input `input_<name>.txt` for each of them (small parameters, a few seconds of
run time), and the scripts `<name>_exact_diagonalization.py` (check) and `plot_<name>.py` (figures).

## Output layout

Every run writes into its own folder, named after all the parameters that change its results, with a
copy of its input file; data and figures are kept apart:

```
<folder>/data/<run>/input.txt                      input of the run: rerun it with ./<name> data/<run>/input.txt
<folder>/data/<run>/observables.txt, nj.txt, ...   TN results
<folder>/data/<run>/<file>_exact_diagonalization.txt   ED result for <file>.txt (from the ED script)
<folder>/plots/<run>/<figure>.png                  figures (from the plot script)
```

e.g. `rydberg_chain_tebd/data/rydberg_N12_M2_V2_2.00_Om_0.300_T5_dt0.05_D64/`. The folder is created by
`make_run_directory` ([../io/output.h](../io/output.h)); the ED and plot scripts find the same folder.

## Build and run

A single [Makefile](Makefile) builds all the examples (and the library, if needed):

```bash
cd examples
make LIBRARY_DIR=/path/to/itensor      # all programs (LIBRARY_DIR: folder of ITensor v3 with options.mk)
make run                               # run every program on its sample input
```

New `.cc` files in an example folder are compiled automatically. To run a single program:

```bash
cd rydberg_chain_tebd
./rydberg_chain_tebd input_rydberg_chain_tebd.txt
```

## Exact diagonalization checks

The Python scripts `<name>_exact_diagonalization.py` are meant to quickly compare the tensor-network
(TN) results of the program `<name>` with exact diagonalization (ED) on small systems. Each script reads
the same input file as the TN program, solves the same model exactly, writes its results next to the
TN files of the run and, if the TN output is present, prints the maximum difference of each observable (expected: of the order of the Trotter /
truncation errors of the TN run):

```bash
cd rydberg_chain_tebd
./rydberg_chain_tebd input_rydberg_chain_tebd.txt
python3 rydberg_chain_tebd_exact_diagonalization.py input_rydberg_chain_tebd.txt
```

`make compare` does this for all the examples. The scripts need numpy, scipy and (most of them) quimb,
and share the helpers in [exact_diagonalization_tools.py](exact_diagonalization_tools.py).
Their cost grows exponentially with the system size: keep N small (about 12 spins, or fewer for the
open systems). `impurity_dynamics_exact_diagonalization.py` uses the free-fermion solution when
`Jxx = Jzzz = 0` (any N), and integrates the Lindblad equation of the full density matrix otherwise
(N <= 8).

Each script also defines `main(input_file)`, which returns the comparisons
(`{output: {column: max |TN - ED|}}`); the integration tests in [../tests](../tests/) use it.

## Plots

Each program `<name>` has a script `plot_<name>.py` that plots its results; the shared code is in
[plot_utils.py](plot_utils.py). Site-resolved quantities (densities, magnetizations, entanglement
entropy of each cut) get a heatmap over sites and time, profiles at a few times and the time trace
of every site; the other observables are plotted against time (or field, site). Exact-diagonalization
results are optional: they are overlaid as dashed lines when present in the run folder (run the
`<name>_exact_diagonalization.py` script first), and the plots work without them. The figures are saved
in `plots/<run>/`:

```bash
cd rydberg_chain_tebd
./rydberg_chain_tebd input_rydberg_chain_tebd.txt
python3 plot_rydberg_chain_tebd.py                     # latest run in data/
python3 plot_rydberg_chain_tebd.py rydberg_N12_M2_V2_2.00_Om_0.300_T5_dt0.05_D64 --show   # a given run, open the figures
```

## Input files

All programs read their parameters from a file in the ITensor `InputGroup` format; missing entries
take the defaults written at the top of `main`. Measurement intervals (`t_measure`, `t_corr`) must be
multiples of `dt` (the programs stop with an error otherwise):

```
input
{
N = 12
dt = 0.05
state = up
}
```

The programs that compute ground states (`find_ground_state`, [../ground_state/dmrg.h](../ground_state/dmrg.h))
also read the DMRG parameters `dmrg_max_dim` (200), `dmrg_cutoff` (1E-12), `dmrg_noise` (1E-6, first
sweeps only), `dmrg_tolerance` (1E-10, relative energy change between sweeps), `dmrg_max_sweeps` (30),
`dmrg_restarts` (1: number of runs from random states, the lowest energy is kept) and `dmrg_seed`
(0: different random states at every run).

## List of examples

### `heisenberg_picture`
<psi_0| Z_j(t) |psi_0> of a spin-1/2 chain with fields along x, y, z, from a product state, in the
Schroedinger picture (evolved state) and in the Heisenberg picture (evolved MPO of Z_j, `heisenberg_step`),
with either TEBD gates (`propagator = gates`) or the first-order MPO toExpH(H, i dt) (`propagator = mpo`;
much slower, keep N <= 6). The two pictures agree up to truncation.
- Inputs: `N`, `Jx`, `Jy`, `Jz`, `hx`, `hy`, `hz`, `config`, `basis` (initial product state, see
  `make_product_state`), `site`, `T`, `dt`, `t_measure`, `max_dim`, `cut_off`, `propagator`.
- Output: t, <Z_j> in the two pictures, bond dimensions of the state and of the operator.

### `spin_chain_ground_state`
Ground state of a spin-1/2 chain with nearest- and next-nearest-neighbour couplings and fields,
H = sum_j (Jx X_j X_{j+1} + Jy Y_j Y_{j+1} + Jz Z_j Z_{j+1}) + sum_j (J2x X_j X_{j+2} + ...) + sum_j (hx X_j + hz Z_j),
with `find_ground_state` (`make_spin_chain_mpo`). With `conserve_sz = 1` the search runs in the sector
S^z = 0 of the Neel state (needs `hx = 0`, `Jx = Jy`, `J2x = J2y`, even `N`). Also the driver of the
ground-state tests ([../tests](../tests/)).
- Inputs: `N`, `Jx`, `Jy`, `Jz`, `J2x`, `J2y`, `J2z`, `hx`, `hz`, `conserve_sz`, `dmrg_*`.
- Output: energy, variance and convergence flag; energy of every restart; profile `<Z_j>`, `<Z_j Z_{j+1}>`,
  entanglement entropy of the bond (j, j+1).

### `rydberg_chain_tebd`
Closed dynamics of a 1D Rydberg chain (interactions up to next-nearest neighbours, anti-blockade
detunings) from a kink state, with 3-site TEBD gates (`make_rydberg_gates_nnn`).
Optional gaussian disorder on the atomic positions reproduces Fig. S1 of arXiv:2309.12392.
- Inputs: `N`, `M` (initial excitations), `V2`, `Omega`, `T`, `dt`, `max_dim`, `cut_off`, `sigmax`, `seed`, `t_measure`.
- Output: fidelity, half-chain entropy, bond dimension; Rydberg densities `n_j(t)`; entanglement entropy
  `S_j(t)` of every cut (j, j+1).

### `leaky_cavity`
N spin-1/2 coupled to a lossy cavity mode (open Dicke or Tavis-Cummings model), optionally with
Rydberg interactions `V`. The density matrix is purified into a bra-ket MPS; the boson reaches every
spin with swap gates.
- Inputs: `N`, `max_occ`, `h`, `g` (in units of the critical coupling), `V`, `kappa`, `theta` (initial
  polar angle of the spins, units of pi), `T`, `dt`, `cut_off`, `max_dim`, `coupling` (`dicke` or `tavis`), `t_measure`.
- Sample inputs: `input_leaky_cavity.txt` (Rydberg atoms) and `input_leaky_cavity_dicke.txt` (V = 0).
- Output: Tr(rho), first-spin magnetizations, photon number, bond dimension; profiles `<X_j>`, `<Z_j>`.
- `sweep_g.sh`: runs a sweep over `g`.

### `cavity_unitary`
Closed dynamics of the Dicke or Tavis-Cummings model (pure state), benchmarked with exact
diagonalization.
- Inputs: `N`, `max_occ`, `omega0`, `h`, `g`, `theta` (initial polar angle, units of pi), `T`, `dt`,
  `cut_off`, `max_dim`, `coupling`, `t_measure`.
- Output: fidelity, `<S^x>/N`, `<S^z>/N`, `<a^dag a>/N`, bond dimension.

### `dissipative_impurity_ising`
Ising chain with a dephasing impurity on the first site (arXiv:2404.04255),
H = -sum Z_j Z_{j+1} + Jxx sum X_j X_{j+1} + Jzzz sum Z_j Z_{j+2} + hx sum X_j, L = sqrt(gamma) Z_1.
- `impurity_dynamics`: from the ground state of H, Lindblad evolution of the purified density matrix
  up to `Tness` (profile `<X_j>`, Tr rho, bond dimension), then the autocorrelation
  `<Z_1(t) Z_1(0)>` up to `T`. Inputs: `N`, `hx`, `Jxx`, `Jzzz`, `gamma`, `Tness`, `T`, `dt`, `max_dim`,
  `cut_off`, `t_measure`, `t_corr`, `hz` (symmetry-breaking field for the ground-state search), `dmrg_*`.
- `ground_state_scan`: DMRG ground states for a range of `hx` (energy, variance, central
  magnetization). Inputs: `N`, `Jxx`, `Jzzz`, `hx_min`, `hx_max`, `dhx`, `hz` (symmetry-breaking
  longitudinal field), `dmrg_*`. With `hz = 0` the ground state is almost degenerate in the ordered
  phase and DMRG may return any combination of the two lowest states: fix `dmrg_seed` for reproducible output.

### `ising_full_counting_statistics`
Ising chain in longitudinal and transverse fields, H = -J sum (X_j X_{j+1} + hx X_j + hz Z_j), and
full counting statistics of the block magnetization (arXiv:2005.01679).
- `ising_quench`: TEBD from a product state along x (`state` = `up`, `down`, `wall`); every `t_measure` the entanglement entropy across each bond and the generating function
  G_l(theta) = <exp(i theta S^x_A)> of blocks of l = 1..N/2 sites.
  Inputs: `N`, `J`, `hx`, `hz`, `T`, `dt`, `max_dim`, `cut_off`, `state`, `t_measure`, `number_points`,
  `max_block_size`.
- `ising_thermal`: thermal state at the energy of |+x...+x> (imaginary-time evolution of the identity)
  and its generating function. Inputs: `N`, `J`, `hx`, `hz`, `dbeta`, `number_points`, `max_block_size`,
  `max_dim`, `cut_off`.
