# Location of the ITensor v3 library source folder (the one containing options.mk).
# Override from the command line if needed: make LIBRARY_DIR=/path/to/itensor
LIBRARY_DIR ?= /home/riccardo/itensor

# Programs of this example (one .cc file each)
APPS = ising_quench ising_thermal

# MyTensorNetworkLib root: the static library lib/libmytn.a is (re)built there if needed
MYTN_DIR=../..

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

MYTN_LIB  = $(MYTN_DIR)/lib/libmytn.a
MYTN_GLIB = $(MYTN_DIR)/lib/libmytn-g.a

build: $(APPS)
debug: $(addsuffix -g,$(APPS))

$(MYTN_LIB): FORCE
	$(MAKE) -C $(MYTN_DIR) build LIBRARY_DIR=$(LIBRARY_DIR)

$(MYTN_GLIB): FORCE
	$(MAKE) -C $(MYTN_DIR) debug LIBRARY_DIR=$(LIBRARY_DIR)

%.o: %.cc
	$(CCCOM) -c $(CCFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

%-g.o: %.cc
	$(CCCOM) -c $(CCGFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

$(APPS): %: %.o $(MYTN_LIB) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $< -o $@ $(MYTN_LIB) $(LIBFLAGS)

$(addsuffix -g,$(APPS)): %-g: %-g.o $(MYTN_GLIB) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $< -o $@ $(MYTN_GLIB) $(LIBGFLAGS)

clean:
	rm -f *.o *.d $(APPS) $(addsuffix -g,$(APPS))

FORCE:

.PHONY: build debug clean FORCE

-include $(wildcard *.d)
