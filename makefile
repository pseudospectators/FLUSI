###########################################################################
## 			Makefile for FLUSI
###########################################################################

# see README for documentation regarding setting environment variables.

# the program file
PROG_FILE = FLUSI.f90

# the other Fortran files
FFILES = cal_nlk.f90 cal_vis.f90 FluidTimeStepper.f90 init_fields.f90 \
	obst_mask.f90 params.f90 save_dd.f90 save_fields.f90 subtr.f90 \
	time_step.f90 trextents.f90
OBJS := $(FFILES:%.f90=%.o)

PPFLAG = -cpp

## set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77)
FC = mpif90
endif

# pre-processor flag
PPFLAG= -cpp

# define variables for the intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp
DIFORT= -DIFORT
FFLAGS += -vec_report0
endif


# this seems to be the only one in use.
MPI_HEADER = mpi_duke_header.f90

PROGRAMS = main

FFT_LIB = $(FFT_ROOT)/lib
FFT_INC = $(FFT_ROOT)/include
P3DFFT_LIB = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include

LDFLAGS = -L$(P3DFFT_LIB) -lp3dfft -lfftw3 -lm  -L$(FFT_LIB)
FFLAGS += -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(DIFORT)

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
