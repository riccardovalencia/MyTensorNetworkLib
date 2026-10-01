# MyTensorNetworkLib

MyTensorNetworkLib is a C++ library for tensor network methods built on top of the C++ ITensor (v3) library. It provides customized methods for one-dimensional tensor networks (Matrix Product States and Matrix Product Operators): preparation of specific initial states, equilibrium properties (DMRG), and dynamics of closed and open systems via the Time Evolving Block Decimation (TEBD) algorithm. It handles spin, bosonic and fermionic degrees of freedom.

## Prerequisites

- C++ ITensor v3 library, see the [installation guide](https://itensor.org/docs.cgi?vers=cppv3&page=install).
- A C++17-compliant compiler.

## Structure

| Module | Umbrella header | Content |
|---|---|---|
| [core](core/) | `core.h` | Model-independent tools: gate containers (`MyBondGate`, ...), `apply_gate`, entanglement entropy, generic initial states and `insert_state`, swap gates, density matrices. |
| [spin_boson](spin_boson/) | `spin_boson.h` | Spin-1/2 chains, Rydberg arrays, PXP, spin-boson / cavity QED models (Dicke, Tavis-Cummings), Lindblad dynamics on purified states, related observables. |
| [spins](spins/) | `transverse_field_ising_chain.h` | Ising chain in longitudinal and transverse fields, full counting statistics of the magnetization. |
| [bosons](bosons/) | `bosons.h` | Bosonic quantum east model: Hamiltonians, DMRG drivers, TEBD (closed and open), bosonic initial states and observables. |
| [fermions](fermions/) | `fermions.h` | Tight-binding Hamiltonians and initial states for spinful fermions. |
| [examples](examples/) | | Complete simulations using the library. |
| [legacy](legacy/) | | Old code kept for reference, not compiled. |

`mytn.h` includes every module. Each module umbrella header also includes `core.h`.

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

Include the header of the module you need and link `lib/libmytn.a` before the ITensor libraries:

```cpp
#include "spin_boson.h"

int main()
{
    auto sites = SpinHalf(10, {"ConserveQNs=", false});
    MPS psi = initial_computational_state(sites, "1100000000");     // |down down up up ...>_z
    // ...
}
```

```bash
g++ -std=c++17 -I/path/to/MyTensorNetworkLib $(ITensor CCFLAGS) main.cc -o main \
    /path/to/MyTensorNetworkLib/lib/libmytn.a $(ITensor LIBFLAGS)
```

The [examples](examples/) contain ready-to-use Makefiles that do this (and rebuild the library when needed).

## Conventions

- Sites are 1-indexed, as in ITensor.
- Spin-1/2: `|0> = |up_z>`, `|1> = |down_z>`. The Rydberg/excitation projector is `n = (1 - Z)/2 = |down_z><down_z|`.
- X, Y, Z denote Pauli matrices; S^a = sigma^a/2 (ITensor `"Sx"`, `"Sy"`, `"Sz"`).
- Gate lists returned by the TEBD builders implement one second-order Trotter step (forward sweep with dt/2, then the reversed sweep).
- Functions taking `MPS*` modify the state in place (e.g. by moving its orthogonality center).
- Open systems are simulated on the vectorized (purified) density matrix. In the impurity geometry the bra occupies sites 1..N in reversed order and the ket sites N+1..2N, so that the dissipative site sits on the central bond.
- Headers contain `using namespace std; using namespace itensor;` for backward compatibility with existing drivers.

## License

[MIT](https://choosealicense.com/licenses/mit/)
