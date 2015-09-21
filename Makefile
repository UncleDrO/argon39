# Makefile for ARGON39
# Optimized for Mac OS X with gfortran, doxygen, dot installed

MAKE = make

FC = gfortran
#FFLAGS = 
#FFLAGS = -Ofast
FFLAGS = -fcheck=all -Wall
F90C = $(FC)
F90FLAGS = 

LDFLAGS =

INSTALL_DIR = ~/Codes/bin

program = argon39.x
headers = 
objects = argon39.o decchains.o elemnucl.o ensdf.o fission.o general.o iostuff.o modules.o neutron.o stoppow.o nr_dqtrap.o nr_dtrapzd.o nr_hunt.o
modules = abra_type.mod abun_type.mod alev_type.mod all_data.mod bbra_type.mod decch_type.mod elem_type.mod ensdf.mod inpout.mod mhunt.mod neutron.mod nucl_type.mod precision.mod
libs = 
files = 

program2 = talysprep.x
headers2 = 
objects2 = talysprep.o elemnucl.o ensdf.o general.o iostuff.o modules.o neutron.o nr_dqtrap.o nr_dtrapzd.o nr_hunt.o
modules2 = 
libs2 =
files2 =

program3 = mcnp6prep.x
headers3 = 
objects3 = mcnp6prep.o modules.o
modules3 = prepare.mod
libs3 =
files3 =

program4 = mcnp6extract.x
headers4 = 
objects4 = mcnp6extract.o modules.o
modules4 = extract.mod
libs4 =
files4 =

allprograms = $(program) $(program2) $(program3) $(program4)
allobjects = $(objects) $(objects2) $(objects3) $(objects4)
allmodules = $(modules) $(modules2) $(modules3) $(modules4)
allfiles = $(files) $(files2) $(files3) $(files4)

%.o : %.mod  # cancel m2c implicit rule

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(F90C) -c $(FFLAGS) $(F90FLAGS) $< -o $@

### BEGIN SPECIFIC TASKS ###

all: $(allprograms)

$(program): $(objects)  
	$(F90C) $(LDFLAGS) -o $@ $(objects) $(libs)

$(program2): $(objects2)  
	$(F90C) $(LDFLAGS) -o $@ $(objects2) $(libs2)

$(program3): $(objects3)  
	$(F90C) $(LDFLAGS) -o $@ $(objects3) $(libs3)

$(program4): $(objects4)  
	$(F90C) $(LDFLAGS) -o $@ $(objects4) $(libs4)

decchains.o ensdf.o elemnucl.o fission.o general.o iostuff.o neutron.o stoppow.o: modules.o
elemnucl.o: ensdf.o
decchains.o general.o stoppow.o: nr_hunt.o

TALYSprep.o: modules.o

MCNP6prep.o: modules.o

MCNP6extract.o: modules.o

.PHONY : install  # moves executable to INSTALL_DIR
install: $(allprograms)
	/bin/cp -p $(allprograms) $(INSTALL_DIR)

.PHONY : doc  # generates documentation using doxygen
doc: Doxyfile
	doxygen Doxyfile

.PHONY : opendoc  # opens HMTL documentation
opendoc: doc
	open doc-html/index.html

.PHONY : tidy  # removes objects, modules, backup files
tidy:
	rm -f $(allobjects) $(allmodules)
	rm -f *~

.PHONY : clean  # runs "tidy" + removes executable, documentation, output files
clean: tidy
	rm -f $(allprograms) $(allfiles)
	rm -rf doc-html
