## Makefile for flusi

<<<<<<< HEAD
# note that FFT_ROOT, P3DFFT_ROOT, and FFLAGS must be set in the shell
# (e.g. in .bashrc).

# the program file
PROG_FILE = flusi.f90
=======
MOD_FILES = share_vars.f90
SUB_FILES = FluidTimeStepper.f90 cof_p3dfft.f90 cal_vis.f90 params.f90\
time_step.f90 cal_nlk.f90 init_fields.f90 save_fields.f90 obst_mask.f90\
save_dd.f90 subtr.f90 trextents.f90 
PROG_FILE = FLUSI.f90
>>>>>>> 938f3f9e2c72ce4274f8f82fc4e17dc4dc438758

# the other Fortran files
FFILES = \
cal_nlk.f90  save_fields.f90 cal_vis.f90  FluidTimeStepper.f90 \
obst_mask.f90  share_vars.f90 cofdx.f90 params.f90 \
subtr.f90 cofdy.f90  init_fields.f90  poisson.f90 time_step.f90 \
cofdz.f90 save_dd.f90 trextents.f90
OBJS := $(FFILES:%.f90=%.o)

# this seems to be the only one in use.
MPI_HEADER = mpi_duke_header.f90

PROGRAMS=main

FC = mpif90

PPFLAG= -cpp

FFT_LIB = $(FFT_ROOT)/lib
FFT_INC = $(FFT_ROOT)/include
P3DFFT_LIB = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include

LDFLAGS = -L$(P3DFFT_LIB) -lp3dfft -lfftw3 -lm  -L$(FFT_LIB)
FFLAGS += -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(DIFORT)

# by default, compile all the programs.
all: $(PROGRAMS)

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
