## Makefile for HIVE
## Build on Babel, Octopus or Duke

MOD_FILES = share_vars.f90
SUB_FILES = FluidTimeStepper.f90 cof_p3dfft.f90 cofdx.f90 cofdy.f90 cofdz.f90 cal_vis.f90 poisson.f90 params.f90\
time_step.f90 cal_nlk.f90 init_fields.f90 save_fields.f90 obst_mask.f90\
save_dd.f90 subtr.f90 trextents.f90 ## obst_pos.f90 move_mask.f90 channel.f90 
PROG_FILE = flusi.f90

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),mesocentre)

FF = mpif90
P3DFFT_ROOT = /home/ethomas/p3dfft.2.3.2_fixed
P3DFFT_LIB = p3dfft
P3DFFT_LOC = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include
FFT_LOC = /home/ethomas/fftw/lib
FFT_INC = /home/ethomas/fftw/include
FFLAGS = 
MPI_HEADER = mpi_duke_header.f90

endif
#------------------------------------------------------------------------------------------------
ifeq ($(CONF),babel)

FF = mpixlf90_r
P3DFFT_ROOT = /homegpfs/rech/sch/rsch509/DMITRY/P3DFFT/p3dfft.2.3.2
P3DFFT_LIB = p3dfft
P3DFFT_LOC = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include
FFT_LOC = /bglocal/pub/fftw/3.1.2/lib
FFT_INC = /bglocal/pub/fftw/3.1.2/include
FFLAGPROF =  
FFLAGS = -WF,"-P" -O4 -qarch=450d -qnostrict -qsave -qtune=450 -qcache=auto -qrealsize=8 -d $(FFLAGPROF)
MPI_HEADER = mpi_babel_header.f90

endif

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),sugiton)

FF = mpif90
P3DFFT_ROOT = /home/tommy/Documents/Research/src/p3dfft.2.3.2_fixed
P3DFFT_LIB = p3dfft
P3DFFT_LOC = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include
FFT_LOC = /home/tommy/fftw/lib
FFT_INC = /home/tommy/fftw/include
FFLAGS = -g
MPI_HEADER = mpi_duke_header.f90

endif

#------------------------------------------------------------------------------------------------
ifeq ($(CONF),duke)

FF = mpiifort
P3DFFT_ROOT = /home/dkolom/P3DFFT/p3dfft.2.3.2
P3DFFT_LIB = p3dfft
P3DFFT_LOC = $(P3DFFT_ROOT)/lib
P3DFFT_INC = $(P3DFFT_ROOT)/include
FFT_LOC = /Softs/fftw/lib
FFT_INC = /Softs/fftw/include
FFLAGS = -O3 -xW -fpp -r8
MPI_HEADER = mpi_duke_header.f90

endif


help:
	@echo Usage: make {babel/sugiton/duke}
babel:
	$(MAKE) build CONF=babel
mesocentre:
	$(MAKE) build CONF=mesocentre
sugiton:
	$(MAKE) build CONF=sugiton
duke:
	$(MAKE) build CONF=duke
build: 
	$(FF) $(FFLAGS) -o main $(MPI_HEADER) $(MOD_FILES) $(SUB_FILES) $(PROG_FILE) -I$(P3DFFT_INC) -I$(FFT_INC) -L$(P3DFFT_LOC) -l$(P3DFFT_LIB) -L$(FFT_LOC) -lfftw3 -lm	
clean:
	rm -f main *.o *.mod
