## Makefile for the flusi post-processing tools.

## set the default compiler if it's not already set, make sure it's not F77.
ifndef FC
FC = mpif90
endif
ifeq ($(FC),f77) # sometimes FC gets defined as f77, which is bad.
FC = mpif90
endif

PROGRAMS = convert_mpiio2vtk convert_mpiio2vtk_ALL convert_mpiio2binary

all: $(PROGRAMS)

convert_mpiio2vtk: convert_vtk.f90 
	$(FC) $(FFLAGS) $^ -o $@

convert_mpiio2vtk_ALL: convert_vtk_ALL.f90
	$(FC) $(FFLAGS) $^ -o $@

convert_mpiio2binary: convert_mpiio.f90
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f *.o *.mod $(PROGRAMS)
