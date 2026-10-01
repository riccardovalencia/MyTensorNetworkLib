# Location of the ITensor v3 library source folder (the one containing options.mk).
# Override from the command line if needed: make LIBRARY_DIR=/path/to/itensor
LIBRARY_DIR ?= /home/riccardo/itensor

APP=rydberg_chain_TEBD

# MyTensorNetworkLib sources used by this example
MYLIB_DIR=../../spin_boson
HEADERS=../../spin_boson.h $(wildcard $(MYLIB_DIR)/*.h)

CCFILES=$(APP).cc \
	$(MYLIB_DIR)/TEBD.cc \
	$(MYLIB_DIR)/initial_state.cc \
	$(MYLIB_DIR)/observables.cc \
	$(MYLIB_DIR)/MyClasses.cc \
	$(MYLIB_DIR)/parameters.cc \
	$(MYLIB_DIR)/general.cc \
	$(MYLIB_DIR)/state_manipulation.cc

#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<
.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	mkdir -p $(dir $@)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs $(OBJECTS) $(APP) $(APP)-g

.PHONY: build debug clean
