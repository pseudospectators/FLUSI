# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

# Non-module Fortran files to be compiled:
FFILES = cal_nlk.f90 cal_vis.f90 FluidTimeStepper.f90 init_fields.f90 \
	mask.f90 mask_fsi.f90 mask_mhd.f90 save_fields.f90 time_step.f90 \
	init_fields_mhd.f90 init_fields_fsi.f90 integrals.f90 params.f90 \
	insects.f90 postprocessing.f90 runtime_control.f90 cal_drag.f90 \
	sponge.f90 FFT_unit_test.f90
OBJS := $(FFILES:%.f90=%.o)

# Set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77)
FC = mpif90
endif

# GNU compiler
ifeq ($(shell $(FC) --version 2>&1 | tail -n 1 | head -c 3),GNU)
FFLAGS += -Wall # warn for unused and uninitialzied variables 
FFLAGS += -Wsurprising # warn if things might not behave as expected
PPFLAG= -cpp #preprocessor flag

# Debug flags for gfortran
#FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp #preprocessor flag
DIFORT= -DIFORT # define the IFORT variable 
FFLAGS += -vec_report0
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
PPFLAG= -qsuffix=cpp=f90  #preprocessor flag
endif

# This seems to be the only one in use.
MPI_HEADER = mpi_duke_header.f90

PROGRAMS = flusi mhd

# FFT_ROOT is set in envinroment.
FFT_LIB = $(FFT_ROOT)/lib
FFT_INC = $(FFT_ROOT)/include

# P3DFFT_ROOT is set in environment.
P3DFFT_LIB = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include

# HDF_ROOT is set in environment.
HDF_LIB = $(HDF_ROOT)/lib
HDF_INC = $(HDF_ROOT)/include

LDFLAGS = -L$(P3DFFT_LIB) -lp3dfft -L$(FFT_LIB) -lfftw3 -lm
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB) -lhdf5_fortran -lhdf5 -lz -ldl

FFLAGS += -I$(HDF_INC) -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(DIFORT)

# Both programs are compiled by default.
all: $(PROGRAMS)

# Compile main programs, with dependencies.
flusi: FLUSI.f90 $(OBJS) mpi_header.o share_vars.o cof_p3dfft.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
mhd: mhd.f90 $(OBJS) mpi_header.o share_vars.o cof_p3dfft.o
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran).
mpi_header.o: $(MPI_HEADER)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
share_vars.o: share_vars.f90 mpi_header.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
cof_p3dfft.o: cof_p3dfft.f90 share_vars.o mpi_header.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
%.o: %.f90 mpi_header.o share_vars.o cof_p3dfft.o
	$(FC) $(FFLAGS) -c -o $@ $<  $(LDFLAGS)

clean:
	rm -f $(PROGRAMS) *.o *.mod
