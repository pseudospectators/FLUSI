## Makefile for HIVE
## Build on Babel, Octopus or Duke

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),mesocentre)

FF = mpif90
FFLAGS = 
MPI_HEADER = mpi_duke_header.f90

endif

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),sugiton)

FF = mpif90
FFLAGS = -g
MPI_HEADER = mpi_duke_header.f90

endif

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),duke)

FF = mpiifort
FFLAGS = -O3 -xW -fpp -r8
MPI_HEADER = mpi_duke_header.f90

endif


help:
	@echo Usage: make {babel/sugiton/duke}
mesocentre:
	$(MAKE) build CONF=mesocentre
sugiton:
	$(MAKE) build CONF=sugiton
duke:
	$(MAKE) build CONF=duke
build: 
	$(FF) $(FFLAGS) $(MPI_HEADER) convert_vtk.f90 -o convert_mpiio2vtk.out 
	$(FF) $(FFLAGS) $(MPI_HEADER) convert_vtk_ALL.f90 -o convert_mpiio2vtk_ALL.out 
	$(FF) $(FFLAGS) $(MPI_HEADER) convert_mpiio.f90 -o convert_mpiio2binary.out
clean:
	rm -f *.out *.o *.mod
