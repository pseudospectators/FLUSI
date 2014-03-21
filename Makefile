# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

# Non-module Fortran files to be compiled:
FFILES = rhs.f90 vis.f90 fluid_time_step.f90 init_fields.f90 \
	mask.f90 mask_fsi.f90 mask_mhd.f90 save_fields.f90 time_step.f90 \
	init_fields_mhd.f90 init_fields_fsi.f90 integrals.f90 params.f90 \
	insects.f90 postprocessing.f90 runtime_control.f90 drag.f90 \
	sponge.f90 fft_unit_test.f90 draw_plate.f90 draw_sphere.f90 \
        kineloader.f90 set_fluid_solid_comms.f90 rigid_solid_time_stepper.f90 \
        add_channel.f90 add_cavity.f90 draw_flexible_plate.f90 \
        wings_geometry.f90 wings_motion.f90 body_motion.f90 \
        body_geometry.f90 rotation_matrices.f90 stroke_plane.f90
        

# Object and module directory:
OBJDIR=obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = vars.f90 kine.f90 cof_p3dfft.f90 solid_solver.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = src
VPATH += :src/inicond:src/inicond/hyd:src/inicond/mhd
VPATH += :src/geometry:src/geometry/hyd:src/geometry/mhd
VPATH += :src/insects:src/solid_solver

# Set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77)
FC = mpif90
endif

# GNU compiler
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
FFLAGS += -Wall # warn for unused and uninitialzied variables 
FFLAGS += -Wsurprising # warn if things might not behave as expected
FFLAGS += -pedantic -fdefault-double-8 -fdefault-real-8
FFLAGS += -Wconversion
FFLGAS += -Wunused-labels
PPFLAG= -cpp #preprocessor flag

# Debug flags for gfortran
FFLAGS += -Wuninitialized -O -fimplicit-none -fbounds-check -g -ggdb
endif

# Intel compiler
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp #preprocessor flag
DIFORT= -DIFORT # define the IFORT variable 
FFLAGS += -module $(OBJDIR) # specify directory for modules.
FFLAGS += -vec_report0
endif

#IBM compiler
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
FFLAGS += -qmoddir=$(OBJDIR)
FFLAGS += -I$(OBJDIR)
PPFLAG= -qsuffix=cpp=f90  #preprocessor flag
endif

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
LDFLAGS += -llapack
FFLAGS += -I$(HDF_INC) -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(DIFORT)


# Both programs are compiled by default.
all: directories $(PROGRAMS) 

# Compile main programs, with dependencies.
flusi: flusi.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)
mhd: mhd.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/vars.o: vars.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/kine.o: kine.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/cof_p3dfft.o: cof_p3dfft.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/solid_solver.o: solid_solver.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
