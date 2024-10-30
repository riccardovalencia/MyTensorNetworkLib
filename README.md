# MyTensorNetworkLib

MyTensorNetworkLib is a C++ library for tensor network methods built on top of the C++ ITensor library. This library provides customized methods for working with one-dimensional tensor networks, including Matrix Product States (MPS) and Matrix Product Operators (MPO) tailored for specific purposes: initialization of specific MPS and MPO; equilibrium properties (DMRG); dynamics via Time Evolving Block Decimation (TEBD) algorithm of closed and open systems. It can handle spin, bosonic, and fermionic degrees of freedom.

## Prerequisites

C++ ITensor v3 library: You need the ITensor v3 library to be installed. See the ITensor Installation Guide for details https://itensor.org/docs.cgi?vers=cppv3&page=install

C++17: This library requires a C++17-compliant compiler.

## Installation
Clone the repository:
```bash
git clone https://github.com/riccardovalencia/MyTensorNetworkLib.git
```


## License

[MIT](https://choosealicense.com/licenses/mit/)
