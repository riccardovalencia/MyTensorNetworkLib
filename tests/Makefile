# Integration tests: tensor-network examples against exact diagonalization (see README.md).
#   make test LIBRARY_DIR=/path/to/itensor     build the library and the examples, run all the tests
#   make test TEST=test_integration.TensorNetworkVsExactDiagonalization.test_leaky_cavity
#                                              run a single test

# ITensor v3 folder (the one containing options.mk)
LIBRARY_DIR ?= /home/riccardo/itensor
TEST ?= test_integration

test:
	$(MAKE) -C ../examples build LIBRARY_DIR=$(LIBRARY_DIR)
	python3 -m unittest -v $(TEST)

.PHONY: test
