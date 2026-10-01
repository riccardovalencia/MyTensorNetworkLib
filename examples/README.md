# Examples

Each example lives in its own folder with a `Makefile` that compiles the driver and links it against
`lib/libmytn.a` (the library is built or updated automatically by `make`).

Set `LIBRARY_DIR` to your ITensor v3 folder (the one containing `options.mk`), either by editing the
`Makefile` or on the command line:

```bash
cd examples/rydberg_chain_TEBD
make LIBRARY_DIR=/path/to/itensor
./rydberg_chain_TEBD            # default parameters
./rydberg_chain_TEBD 20 2 2.0 0.1 10 0.05 128   # N M V2 Omega T dt maxDim
./rydberg_chain_TEBD 20 2 2.0 0.1 10 0.05 128 0.01 3   # ... + position disorder sigmax, seed
```

Examples reading parameters from an input file are run as `./<app> input.txt` (an `input.txt` with
small, quick parameters is provided; missing entries take the defaults written in the source).

## List of examples

- `rydberg_chain_TEBD`: closed dynamics of a 1D Rydberg chain (interactions up to next-nearest
  neighbours, anti-blockade detunings) starting from a kink state, using 3-site TEBD gates
  (`gates_rydberg_up_to_VNNN`). Optional gaussian disorder on the atomic positions (`sigmax`, `seed`),
  reproducing Fig. S1 of https://arxiv.org/abs/2309.12392. Writes fidelity, half-chain entanglement entropy, max bond dimension
  and Rydberg densities `n_j(t)` to `data/`.

- `rydberg_in_leaky_cavity` (`./tn_rydberg_in_leaky_cavity input.txt`): N Rydberg atoms with
  nearest-neighbour interactions coupled to a single lossy cavity mode (open Dicke model with Rydberg
  interactions). The density matrix is purified into a doubled bra-ket MPS; long-range photon-matter
  gates are applied with swap gates. Writes norm, total `Sx`, `Sz`, photon number, max bond dimension,
  and local `sx_j`, `sz_j` to `data/`. `sub_tn_rydberg_in_leaky_cavity.sh` runs a sweep over the
  coupling `g`.
- `collective_light_matter_unitary` (`./tn_unitary_collective_light_matter_systems`): closed dynamics
  of the Dicke (or Tavis-Cummings, via `photon_matter_coupling`) model starting from a spin coherent
  state. Parameters are set at the top of `main`. Benchmarked against exact diagonalization.
- `collective_light_matter_dissipative` (`./tn_leaky_collective_light_matter_systems_dynamics input.txt`):
  Dicke model with photon losses, solved via the Lindblad master equation in the purified (doubled)
  space.

- `dissipative_impurity_ising`: Ising chain with a dephasing impurity on the first site,
  H = -sum_j Z_j Z_{j+1} + Jxx sum_j X_j X_{j+1} + hx sum_j X_j (optionally + Jzzz sum_j Z_j Z_{j+2}),
  L = sqrt(gamma) Z_1 (https://arxiv.org/abs/2404.04255). The density matrix is purified on 2N sites
  (bra mirrored on sites 1..N, ket on N+1..2N) and evolved from the ground state of H (DMRG).
  All programs write to `data/`:
  - `./tn_ising_model_impurity N hx Jxx gamma Tness T dt maxDim`: evolves up to `Tness` saving the state
    every 0.2 time units, then computes the autocorrelation <Z_1(t) Z_1(0)> up to time `T`.
  - `./tn_ising_model_impurity_NNN_interactions N hx Jxx Jzzz gamma Tness T dt maxDim`: same with the
    next-nearest-neighbour coupling `Jzzz` (three-site gates).
  - `./tn_ising_model_impurity_measure ...` and `./tn_ising_model_impurity_measure_NNN ...` (same
    arguments as the corresponding dynamics): post-processing of the saved states (local
    magnetizations, trace, bond dimension).
  - `./DMRG_Ising_up_NNN N hx_max Jxx Jzzz`: ground states and their energy variance and
    magnetization for hx from 1.1 to hx_max.
  - `./tn_ising_model_pure N hx Jxx T dt maxDim`: coherent TEBD of a pure state starting from the
    ground state of the same Hamiltonian (a check: the state must stay stationary).

  Note: the autocorrelation part of `tn_ising_model_impurity_measure_NNN` reads states
  (`..._psi_autocorr_t...`) that the dynamics programs do not save, so its output file stays empty.

- `ising_full_counting_statistics`: quench in the Ising chain in longitudinal and transverse fields,
  H = -J sum_j (X_j X_{j+1} + hx X_j + hz Z_j), and full counting statistics of the subsystem
  magnetization (https://arxiv.org/abs/2005.01679). The programs share files in `data/`:
  - `./TEBD_TLIC state N J hxChoice hzChoice ttotal tstep nmeas bonddim 0`: TEBD from |+x...+x>
    (state 0), |-x...-x> (1) or a domain wall (2); saves the state every `nmeas` steps.
    hx and hz are chosen by index from the lists in `get_data` (io/input.h); the last argument is
    unused (it selected local/cluster runs in the original code).
  - `./FULL_COUNTING_STATISTICS N hxChoice hzChoice tstep nmeas numberPoints maxLength 0`: generating
    function G(theta) = <exp(i theta S^x_A)> of the saved states, for blocks of 1..maxLength sites.
  - `./ENTROPY_HALF N tstep nmeas 0`, `./ENTROPY_ALL N tstep nmeas 0`: entanglement entropy (natural
    log) of the saved states across the central bond / every bond.
  - `./TLIC_compute_thermal N hx hz`: thermal state at the energy of |+x...+x>, by imaginary-time
    evolution of the identity, and its generating function.

  Typical run:
  ```bash
  ./TEBD_TLIC 0 8 1 1 3 2 0.01 50 64 0
  ./FULL_COUNTING_STATISTICS 8 1 3 0.01 50 40 4 0
  ./ENTROPY_HALF 8 0.01 50 0
  ```
