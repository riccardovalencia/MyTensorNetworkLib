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

## List of examples

- `rydberg_chain_TEBD`: closed dynamics of a 1D Rydberg chain (interactions up to next-nearest
  neighbours, anti-blockade detunings) starting from a kink state, using 3-site TEBD gates
  (`gates_rydberg_up_to_VNNN`). Writes fidelity, half-chain entanglement entropy, max bond dimension
  and Rydberg densities `n_j(t)` to `data/`.
