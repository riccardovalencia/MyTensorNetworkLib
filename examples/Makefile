# Build all the examples in one go:
#   make LIBRARY_DIR=/path/to/itensor          (release, default target)
#   make debug LIBRARY_DIR=/path/to/itensor    (programs <name>-g with debug symbols)
#   make run                                   (run every program on its sample input)
#
# Every examples/<folder>/<name>.cc is compiled into the program examples/<folder>/<name> and
# linked against MyTensorNetworkLib, which is (re)built in ../lib when needed. New examples are
# picked up automatically. Run a program from its folder, with an input file:
#   cd rydberg_chain_tebd && ./rydberg_chain_tebd input_rydberg_chain_tebd.txt
# Results are written to <folder>/data/.

# ITensor v3 folder (the one containing options.mk)
LIBRARY_DIR ?= /home/riccardo/itensor

# MyTensorNetworkLib root
MYTN_DIR = ..

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

SOURCES  = $(wildcard */*.cc)
PROGRAMS = $(SOURCES:.cc=)
GPROGRAMS = $(addsuffix -g,$(PROGRAMS))

MYTN_LIB  = $(MYTN_DIR)/lib/libmytn.a
MYTN_GLIB = $(MYTN_DIR)/lib/libmytn-g.a

build: $(PROGRAMS)
debug: $(GPROGRAMS)

$(MYTN_LIB): FORCE
	$(MAKE) -C $(MYTN_DIR) build LIBRARY_DIR=$(LIBRARY_DIR)

$(MYTN_GLIB): FORCE
	$(MAKE) -C $(MYTN_DIR) debug LIBRARY_DIR=$(LIBRARY_DIR)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

%-g.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

$(PROGRAMS): %: %.o $(MYTN_LIB) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $< -o $@ $(MYTN_LIB) $(LIBFLAGS)

$(GPROGRAMS): %-g: %-g.o $(MYTN_GLIB) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $< -o $@ $(MYTN_GLIB) $(LIBGFLAGS)

# run every program <folder>/<name> on <folder>/input_<name>.txt
run: build
	@for p in $(PROGRAMS); do \
		d=$$(dirname $$p); n=$$(basename $$p); \
		echo "== $$p input_$$n.txt"; \
		(cd $$d && ./$$n input_$$n.txt > /dev/null 2>&1) || { echo "FAILED: $$p"; exit 1; }; \
	done

clean:
	rm -f */*.o */*.d $(PROGRAMS) $(GPROGRAMS)

FORCE:

.PHONY: build debug run clean FORCE

-include $(wildcard */*.d)
