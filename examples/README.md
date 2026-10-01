# Examples

Each example lives in its own folder with a `Makefile` that compiles the driver together with the
needed sources of MyTensorNetworkLib (no copies of the library are required).

Set `LIBRARY_DIR` to your ITensor v3 folder (the one containing `options.mk`), either by editing the
`Makefile` or on the command line:

```bash
cd examples/rydberg_chain_TEBD
make LIBRARY_DIR=/path/to/itensor
./rydberg_chain_TEBD            # default parameters
./rydberg_chain_TEBD 20 2 2.0 0.1 10 0.05 128   # N M V2 Omega T dt maxDim
```

Examples reading parameters from an input file are run as `./<app> input.txt` (an `input.txt` with
small, quick parameters is provided; missing entries take the defaults written in the source).

## List of examples

- `rydberg_chain_TEBD`: closed dynamics of a 1D Rydberg chain (interactions up to next-nearest
  neighbours, anti-blockade detunings) starting from a kink state, using 3-site TEBD gates
  (`gates_rydberg_up_to_VNNN`). Writes fidelity, half-chain entanglement entropy, max bond dimension
  and Rydberg densities `n_j(t)` to `data/`.

- `rydberg_in_leaky_cavity` (`./tn_rydberg_in_leaky_cavity input.txt`): N Rydberg atoms with
  nearest-neighbour interactions coupled to a single lossy cavity mode (open Dicke model with Rydberg
  interactions). The density matrix is purified into a doubled bra-ket MPS; long-range photon-matter
  gates are applied with swap gates. Writes norm, total `Sx`, `Sz`, photon number, max bond dimension,
  and local `sx_j`, `sz_j` to `data/`.
- `collective_light_matter_unitary` (`./tn_unitary_collective_light_matter_systems`): closed dynamics
  of the Dicke (or Tavis-Cummings, via `photon_matter_coupling`) model starting from a spin coherent
  state. Parameters are set at the top of `main`. Benchmarked against exact diagonalization.
- `collective_light_matter_dissipative` (`./tn_leaky_collective_light_matter_systems_dynamics input.txt`):
  Dicke model with photon losses, solved via the Lindblad master equation in the purified (doubled)
  space.
