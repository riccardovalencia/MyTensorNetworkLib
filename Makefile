# Build MyTensorNetworkLib as a static library: lib/libmytn.a (and lib/libmytn-g.a with `make debug`).
#
# LIBRARY_DIR must point to the ITensor v3 folder containing options.mk:
#   make LIBRARY_DIR=/path/to/itensor
# The compiler and flags are taken from ITensor's options.mk, so the library is built
# consistently with ITensor itself.

LIBRARY_DIR ?= /home/riccardo/itensor

MODULES = core spin_boson spins bosons fermions

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

SOURCES  = $(wildcard $(addsuffix /*.cc,$(MODULES)))
OBJECTS  = $(patsubst %.cc,build/release/%.o,$(SOURCES))
GOBJECTS = $(patsubst %.cc,build/debug/%.o,$(SOURCES))

LIB  = lib/libmytn.a
GLIB = lib/libmytn-g.a

build: $(LIB)
debug: $(GLIB)
all: build debug

build/release/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CCCOM) -c $(CCFLAGS) -MMD -MP -o $@ $<

build/debug/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CCCOM) -c $(CCGFLAGS) -MMD -MP -o $@ $<

$(LIB): $(OBJECTS)
	@mkdir -p lib
	rm -f $@
	ar rcs $@ $^

$(GLIB): $(GOBJECTS)
	@mkdir -p lib
	rm -f $@
	ar rcs $@ $^

clean:
	rm -rf build lib

.PHONY: build debug all clean

-include $(OBJECTS:.o=.d) $(GOBJECTS:.o=.d)
