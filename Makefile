# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

# Non-module Fortran files to be compiled:
FFILES = rhs.f90 vis.f90 fluid_time_step.f90 init_fields.f90 \
	mask.f90 mask_fsi.f90 mask_mhd.f90 save_fields.f90 time_step.f90 \
	init_fields_mhd.f90 init_fields_fsi.f90 integrals.f90 params.f90 \
	runtime_control.f90 drag.f90 \
	sponge.f90 fft_unit_test.f90 draw_plate.f90 draw_sphere.f90 \
  rotation_matrices.f90 add_channel.f90 add_cavity.f90 init_scalar.f90 dry_run.f90 \
  noncircular_cylinder.f90 draw_flexible_plate.f90 \
	bin2hdf.f90 compare_key.f90 compare_timeseries.f90 convert_abs_vorticity.f90 \
	convert_velocity.f90 convert_vorticity.f90 copy_hdf_file.f90 energy_post.f90 \
	extract_subset.f90 field_analysis.f90 hdf2bin.f90 keyvalues.f90 \
	magnitude_post.f90 mean_2d.f90 post_helicity.f90 postprocessing.f90 \
	post_spectrum.f90 pressure_to_Qcriterion.f90 set_hdf5_attribute.f90 \
	simple_field_operation.f90 time_avg_HDF5.f90 tke_mean.f90 \
	turbulence_analysis.f90 upsample.f90 stl2dist.f90 dist2chi.f90 force_decomposition.f90 \
	extend_domain.f90 basic_file_routines.f90 fractal_trees.f90 convert_to_wing_system.f90 \
	pressure_force.f90 flusi_hdf5_interface.f90 runtime_backuping.f90 io_test.f90 transpose_test.f90


# Object and module directory:
OBJDIR=obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

# Files that create modules:
MFILES = vars.f90 helpers.f90 cof_p3dfft.f90 solid_solver.f90 \
	interpolation.f90 basic_operators.f90 insects.f90 turbulent_inlet.f90 \
	slicing.f90 ghostpoints.f90 passive_scalar.f90 ini_files_parser.f90 stlreader.f90 \
	ini_files_parser_mpi.f90 hdf5_wrapper.f90 wavelet_library.f90
MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

# Source code directories (colon-separated):
VPATH = src
VPATH += :src/inicond:src/inicond/hyd:src/inicond/mhd:src/inicond/scalar
VPATH += :src/geometry:src/geometry/hyd:src/geometry/mhd:src/wavelets
VPATH += :src/insects:src/solid_solver:src/postprocessing:src/file_io

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
#FFLAGS += -Wall # warn for unused and uninitialzied variables
#FFLAGS += -Werror # warnings are errors
#FFLAGS += -pedantic
#FFLAGS += -O3
PPFLAG= -cpp #preprocessor flag
# Debug flags for gfortran:
FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -g -ggdb
FFLAGS += -O3
FFLAGS += -ffpe-trap=zero,overflow,underflow -ffree-line-length-none
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
DIFORT=-WF,-DTURING # this defines the TURING with the IBM compiler
PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
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
$(OBJDIR)/cof_p3dfft.o: cof_p3dfft.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/insects.o: insects.f90 $(OBJDIR)/vars.o \
	body_geometry.f90 body_motion.f90 rigid_solid_time_stepper.f90 wings_geometry.f90 \
	wings_motion.f90 stroke_plane.f90 type_definitions.f90\
  kineloader.f90 $(OBJDIR)/helpers.o $(OBJDIR)/ini_files_parser_mpi.o $(OBJDIR)/stlreader.o $(OBJDIR)/interpolation.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/solid_solver.o: solid_solver.f90 $(OBJDIR)/vars.o  $(OBJDIR)/interpolation.o $(OBJDIR)/basic_operators.o $(OBJDIR)/insects.o $(OBJDIR)/helpers.o \
	mouvement.f90 integrate_position.f90 init_beam.f90 save_beam.f90 BeamForces.f90 plate_geometry.f90 aux.f90 \
	prescribed_beam.f90 solid_solver_wrapper.f90 $(OBJDIR)/ghostpoints.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ghostpoints.o: ghostpoints.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/interpolation.o: interpolation.f90 $(OBJDIR)/vars.o $(OBJDIR)/basic_operators.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/basic_operators.o: basic_operators.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/turbulent_inlet.o: turbulent_inlet.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o $(OBJDIR)/basic_operators.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/slicing.o: slicing.f90 $(OBJDIR)/vars.o $(OBJDIR)/hdf5_wrapper.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/helpers.o: helpers.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser_mpi.o: ini_files_parser_mpi.f90 $(OBJDIR)/vars.o $(OBJDIR)/ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser.o: ini_files_parser.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/stlreader.o: stlreader.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/hdf5_wrapper.o: hdf5_wrapper.f90 $(OBJDIR)/vars.o $(OBJDIR)/helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/passive_scalar.o: passive_scalar.f90 $(OBJDIR)/vars.o $(OBJDIR)/basic_operators.o $(OBJDIR)/ghostpoints.o $(OBJDIR)/helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/wavelet_library.o: wavelet_library.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod a.out

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
