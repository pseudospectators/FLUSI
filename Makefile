# Makefile for FLUSI

# see README for documentation regarding setting environment variables.

# the program file
PROG_FILE = FLUSI.f90

# the other Fortran files
FFILES = cal_nlk.f90 cal_vis.f90 FluidTimeStepper.f90 init_fields.f90 \
	create_mask.f90 params.f90 save_fields.f90 time_step.f90
OBJS := $(FFILES:%.f90=%.o)

## set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77) # sometimes FC gets defined as f77, which is bad.
FC = mpif90
endif

# GNU compiler
ifeq ($(shell $(FC) -v 2>&1 | tail -n 1 | head -c 3),gcc)
FFLAGS += -Wall # warn for unused and uninitialzied variables 
FFLAGS += -Wsurprising # warn if things might not behave as expected
PPFLAG= -cpp #preprocessor flag
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
# We are using the Intel compiler
PPFLAG= -fpp #preprocessor flag
DIFORT= -DIFORT # define the IFORT variable 
FFLAGS += -vec_report0
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
PPFLAG= -qsuffix=cpp=f90  #preprocessor flag
endif

# this seems to be the only one in use.
MPI_HEADER = mpi_duke_header.f90

PROGRAMS = main

FFT_LIB = $(FFT_ROOT)/lib
FFT_INC = $(FFT_ROOT)/include
P3DFFT_LIB = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include

# for mesocentre:
# export HDF_ROOT=/LOGINSTO/softs/hdf5/intel/1.8.8/

HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include

LDFLAGS = -L$(P3DFFT_LIB)  -lp3dfft -lfftw3 -lm  -L$(FFT_LIB)  $(HDF5_FLAGS)
LDFLAGS += -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz

FFLAGS += -I$(HDF_INC) -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(DIFORT)


# by default, compile all the programs.
all: main

# compile main program, with dependencies:
main: $(PROG_FILE) $(OBJS) mpi_header.o share_vars.o cof_p3dfft.o
	$(FC) -o $@ $^ $(LDFLAGS)

# compile modules:
mpi_header.o: $(MPI_HEADER)
	$(FC) $(FFLAGS) -c -o $@ $^ $(LDFLAGS)
share_vars.o: share_vars.f90
	$(FC) $(FFLAGS) -c -o $@ $^ $(LDFLAGS)
cof_p3dfft.o: cof_p3dfft.f90 share_vars.o mpi_header.o
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

# compile remaining objects:
%.o: %.f90 mpi_header.o share_vars.o cof_p3dfft.o
	$(FC) $(FFLAGS) -c $<  $(LDFLAGS) -o $@

clean:
	rm -f main *.o *.mod
