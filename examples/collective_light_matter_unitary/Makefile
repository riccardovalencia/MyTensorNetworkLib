# Location of the ITensor v3 library source folder (the one containing options.mk).
# Override from the command line if needed: make LIBRARY_DIR=/path/to/itensor
LIBRARY_DIR ?= /home/riccardo/itensor

APP=tn_unitary_collective_light_matter_systems

# MyTensorNetworkLib root: the static library lib/libmytn.a is (re)built there if needed
MYTN_DIR=../..

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

MYTN_LIB  = $(MYTN_DIR)/lib/libmytn.a
MYTN_GLIB = $(MYTN_DIR)/lib/libmytn-g.a

build: $(APP)
debug: $(APP)-g

$(MYTN_LIB): FORCE
	$(MAKE) -C $(MYTN_DIR) build LIBRARY_DIR=$(LIBRARY_DIR)

$(MYTN_GLIB): FORCE
	$(MAKE) -C $(MYTN_DIR) debug LIBRARY_DIR=$(LIBRARY_DIR)

$(APP).o: $(APP).cc
	$(CCCOM) -c $(CCFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

$(APP)-g.o: $(APP).cc
	$(CCCOM) -c $(CCGFLAGS) -I$(MYTN_DIR) -MMD -MP -o $@ $<

$(APP): $(APP).o $(MYTN_LIB) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(APP).o -o $(APP) $(MYTN_LIB) $(LIBFLAGS)

$(APP)-g: $(APP)-g.o $(MYTN_GLIB) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(APP)-g.o -o $(APP)-g $(MYTN_GLIB) $(LIBGFLAGS)

clean:
	rm -f $(APP).o $(APP)-g.o $(APP).d $(APP)-g.d $(APP) $(APP)-g

FORCE:

.PHONY: build debug clean FORCE

-include $(APP).d $(APP)-g.d
