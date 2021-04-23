# Makefile for fsi and mhd codes. See README for necessary environment
# variables.

#--------------------------------------------------------------
# Non-module Fortran files to be compiled:
#--------------------------------------------------------------
FFILES = rhs.f90 vis.f90 fluid_time_step.f90 init_fields.f90 \
	mask.f90 mask_fsi.f90 mask_mhd.f90 save_fields.f90 time_step.f90 \
	init_fields_mhd.f90 init_fields_fsi.f90 integrals.f90 params.f90 \
	runtime_control.f90 drag.f90 uCT_assemble_HDF5.f90 \
	sponge.f90 fft_unit_test.f90 draw_plate.f90 draw_sphere.f90 \
	rotation_matrices.f90 add_channel.f90 add_cavity.f90 init_scalar.f90 dry_run.f90 \
	noncircular_cylinder.f90 draw_flexible_plate.f90 \
	runtime_backuping.f90 io_test.f90 POD.f90 post_force.f90

ifndef NOHDF5
# Case WITH HDF5 (all machines except earth simulator)
FFILES += postprocessing.f90 \
	bin2hdf.f90 compare_key.f90 compare_timeseries.f90 \
	convert_velocity.f90 convert_vorticity.f90 copy_hdf_file.f90 energy_post.f90 \
	extract_subset.f90 field_analysis.f90 hdf2bin.f90 keyvalues.f90 \
	magnitude_post.f90 mean_2d.f90 post_helicity.f90 \
	post_spectrum.f90 pressure_to_Qcriterion.f90 set_hdf5_attribute.f90 \
	simple_field_operation.f90 time_avg_HDF5.f90 tke_mean.f90 \
	turbulence_analysis.f90 upsample.f90 stl2dist.f90 dist2chi.f90 force_decomposition.f90 \
	zero_padd.f90 convert_to_wing_system.f90 convert_pressure.f90 \
	pressure_force.f90 flusi_hdf5_interface.f90 post_CVE.f90 transpose_test.f90 \
	post_grad.f90 pointcloud2mask.f90 crop.f90 flexible_wing_mask.f90 post_div.f90

else
# Case WITHOUT Hdf
FFILES += flusi_nohdf5_interface.f90 postprocessing_nohdf5.f90
endif
# Object and module directory:
OBJDIR=obj
OBJS := $(FFILES:%.f90=$(OBJDIR)/%.o)

#--------------------------------------------------------------
# Files that create modules:
#--------------------------------------------------------------

MFILES = vars.f90 module_helpers.f90 cof_p3dfft.f90 solid_solver.f90 flexible_solver.f90\
	interpolation.f90 basic_operators.f90 module_insects.f90 turbulent_inlet.f90 \
	ghostpoints.f90 passive_scalar.f90 ini_files_parser.f90 \
	ini_files_parser_mpi.f90 wavelet_library.f90 module_insects_integration_flusi_wabbit.f90 \
	module_timing.f90

ifndef NOHDF5
MFILES += hdf5_wrapper.f90 slicing.f90 stlreader.f90
else
MFILES += slicing_nohdf5.f90
endif

MOBJS := $(MFILES:%.f90=$(OBJDIR)/%.o)

#--------------------------------------------------------------
# Source code directories (colon-separated):
#--------------------------------------------------------------
VPATH = src
VPATH += :src/inicond:src/inicond/hyd:src/inicond/mhd:src/inicond/scalar
VPATH += :src/geometry:src/geometry/hyd:src/geometry/mhd:src/wavelets
VPATH += :src/insects:src/solid_solver:src/postprocessing:src/file_io
VPATH += :src/flexible_solver

# Set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77)
FC = mpif90
endif
ifeq ($(FC),sxf90)
FC = sxmpif90
endif

#-------------------------------------------------------------------------------
# SX compiler
#-------------------------------------------------------------------------------
ifeq ($(FC),sxmpif90)
FFLAGS += -I$(OBJDIR)
FFLAGS += -R0 -P stack -C hopt -pi nest=2 -pi exp="periodize_coordinate,wave_z,wave_y" expin="src/vars.f90,src/cof_p3dfft.f90" -f2003
FFLAGS += -Wf,"-L fmtlist,summary,transform -pvctl fullmsg"
FFLAGS += -gmalloc
#FFLAGS += -ftrace
PPFLAG =
# the DNOHDF5 flag disables all HDF5 compilation (if present)
ifdef NOHDF5
PRE_FLAGS=-DNOHDF5
endif

# Note that shell $(FC) makes an error on FC system and should be bypassed
else

#-------------------------------------------------------------------------------
# GNU compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) --version 2>&1 | head -n 1 | head -c 3),GNU)
# Specify directory for compiled modules:
FFLAGS += -J$(OBJDIR) # specify directory for modules.
#FFLAGS += -Wall # warn for unused and uninitialzied variables
#FFLAGS += -Werror # warnings are errors
#FFLAGS += -pedantic
PPFLAG= -cpp #preprocessor flag
# Debug flags for gfortran:
#FFLAGS += -Wuninitialized -fimplicit-none -fbounds-check -g -ggdb -fbacktrace
FFLAGS += -O3
#FFLAGS += -ffpe-trap=zero,overflow,underflow -ffree-line-length-none -fbacktrace
ifdef NOHDF5
PRE_FLAGS=-DNOHDF5
endif
endif

#-------------------------------------------------------------------------------
# Intel compiler
#-------------------------------------------------------------------------------
ifort:=$(shell $(FC) --version | head -c 5)
ifeq ($(ifort),ifort)
PPFLAG= -fpp # preprocessor flag
PRE_FLAGS= -DIFORT # define the IFORT variable
ifdef NOHDF5
PRE_FLAGS += -DNOHDF5
endif
FFLAGS += -module $(OBJDIR) # specify directory for modules.
#FFLAGS += -vec_report0
# debug flags
# FFLAGS+= -g -debug all -traceback -check all -CB
endif

#-------------------------------------------------------------------------------
# IBM compiler
#-------------------------------------------------------------------------------
ifeq ($(shell $(FC) -qversion 2>&1 | head -c 3),IBM)
FFLAGS += -qmoddir=$(OBJDIR)
FFLAGS += -I$(OBJDIR)
PRE_FLAGS = -WF,-DTURING
# this defines the TURING with the IBM compiler
ifdef NOHDF5
PRE_FLAGS:=$(PRE_FLAGS),-DNOHDF5
endif
PPFLAG=-qsuffix=cpp=f90  #preprocessor flag
endif


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

# Common build flags
LDFLAGS = -L$(P3DFFT_LIB) -lp3dfft -L$(FFT_LIB) -lfftw3
ifndef NOHDF5
LDFLAGS += $(HDF5_FLAGS) -L$(HDF_LIB)
LDFLAGS += -lhdf5_fortran -lhdf5
endif
ifdef MKLROOT
# Use MKL lapack
LDFLAGS += -I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
else
# Use custom LAPACK installation
LDFLAGS += -llapack -lblas
endif


# SX build flags
ifeq ($(FC),sxmpif90)
LDFLAGS += -L${ASL_ROOT} -laslf90 -lasl -lblas
else
LDFLAGS += -lz -ldl
endif
LDFLAGS += -lm

# Other common compile flags
FFLAGS += -I$(P3DFFT_INC) -I$(FFT_INC) $(PPFLAG) $(PRE_FLAGS)

# HDF5 compile flags
ifndef NOHDF5
FFLAGS += -I$(HDF_INC)
endif



# Both programs are compiled by default.
all: directories $(PROGRAMS)

# Compile main programs, with dependencies.
flusi: flusi.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

mhd: mhd.f90 $(MOBJS) $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Compile modules (module dependency must be specified by hand in
# Fortran). Objects are specified in MOBJS (module objects).
$(OBJDIR)/vars.o: vars.f90 $(OBJDIR)/module_timing.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_timing.o: module_timing.f90
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/cof_p3dfft.o: cof_p3dfft.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_insects.o: module_insects.f90 active_grid_winglets.f90 $(OBJDIR)/module_insects_integration_flusi_wabbit.o \
	body_geometry.f90 body_motion.f90 rigid_solid_time_stepper.f90 wings_geometry.f90 \
	wings_motion.f90 stroke_plane.f90 pointcloud.f90 fractal_trees.f90 insect_init_clean.f90 \
	kineloader.f90
	$(FC) -Isrc/insects/ $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_insects_integration_flusi_wabbit.o: module_insects_integration_flusi_wabbit.f90 $(OBJDIR)/vars.o $(OBJDIR)/module_helpers.o $(OBJDIR)/ini_files_parser_mpi.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/solid_solver.o: solid_solver.f90 $(OBJDIR)/vars.o  $(OBJDIR)/interpolation.o $(OBJDIR)/basic_operators.o $(OBJDIR)/module_insects.o $(OBJDIR)/module_helpers.o \
	mouvement.f90 integrate_position.f90 init_beam.f90 save_beam.f90 BeamForces.f90 plate_geometry.f90 aux.f90 \
	prescribed_beam.f90 solid_solver_wrapper.f90 $(OBJDIR)/ghostpoints.o
	$(FC) -Isrc/solid_solver/ $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/flexible_solver.o: flexible_solver.f90 $(OBJDIR)/vars.o  $(OBJDIR)/interpolation.o $(OBJDIR)/basic_operators.o $(OBJDIR)/module_helpers.o \
	flexible_tri_mask.f90 internal_force.f90 internal_force_derivative.f90 \
	init_wing.f90 flexible_solver_wrapper.f90
	$(FC) -Isrc/flexible_solver/ $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ghostpoints.o: ghostpoints.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/interpolation.o: interpolation.f90 $(OBJDIR)/vars.o $(OBJDIR)/basic_operators.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/basic_operators.o: basic_operators.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o $(OBJDIR)/ghostpoints.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/turbulent_inlet.o: turbulent_inlet.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o $(OBJDIR)/basic_operators.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/slicing.o: slicing.f90 $(OBJDIR)/vars.o $(OBJDIR)/hdf5_wrapper.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/slicing_nohdf5.o: slicing_nohdf5.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/module_helpers.o: module_helpers.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser_mpi.o: ini_files_parser_mpi.f90 $(OBJDIR)/vars.o $(OBJDIR)/ini_files_parser.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/ini_files_parser.o: ini_files_parser.f90 $(OBJDIR)/vars.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/stlreader.o: stlreader.f90 $(OBJDIR)/vars.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/hdf5_wrapper.o: hdf5_wrapper.f90 $(OBJDIR)/vars.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/passive_scalar.o: passive_scalar.f90 $(OBJDIR)/vars.o $(OBJDIR)/basic_operators.o $(OBJDIR)/ghostpoints.o $(OBJDIR)/module_helpers.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/params.o: params.f90 $(OBJDIR)/vars.o $(OBJDIR)/ini_files_parser.o $(OBJDIR)/module_insects.o $(OBJDIR)/solid_solver.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/runtime_backuping.o: runtime_backuping.f90 $(OBJDIR)/vars.o $(OBJDIR)/solid_solver.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/io_test.o: io_test.f90 $(OBJDIR)/vars.o $(OBJDIR)/runtime_backuping.o
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)
$(OBJDIR)/wavelet_library.o: wavelet_library.f90 $(OBJDIR)/vars.o $(OBJDIR)/cof_p3dfft.o coherent_vortex_extraction.f90 FWT3_PO.f90 \
	IWT3_PO.f90
	$(FC) -Isrc/wavelets/ $(FFLAGS) -c -o $@ $< $(LDFLAGS)

# Compile remaining objects from Fortran files.
$(OBJDIR)/%.o: %.f90 $(MOBJS)
	$(FC) $(FFLAGS) -c -o $@ $< $(LDFLAGS)

clean:
	rm -rf $(PROGRAMS) $(OBJDIR)/*.o $(OBJDIR)/*.mod *.mod *.o a.out *.L

tidy:
	rm -rf $(OBJDIR)/*.o $(OBJDIR)/*.mod *.mod *.o a.out

.PHONY: test

test:
	@cd tests/; ./unittest.sh

# If the object directory doesn't exist, create it.
.PHONY: directories

directories: ${OBJDIR}

${OBJDIR}:
	mkdir -p ${OBJDIR}
