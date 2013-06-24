program MHD3D
  use mpi_header
  use mhd_vars
  implicit none

  integer                :: mpicode
  character (len=80)     :: infile

  ! Arrays needed for simulation  
  real(kind=pr),dimension(:,:,:),allocatable :: workvis  
  ! ub and ubk, and nlk are 4-dimensional arrays, with the last index
  ! indicating the field.  The first three fields are for the
  ! velocity, the last three are for the magnetic field.
  complex(kind=pr),dimension(:,:,:,:),allocatable :: ubk
  real(kind=pr),dimension(:,:,:,:),allocatable :: ub
  complex(kind=pr),dimension(:,:,:,:,:),allocatable :: nlk
  ! wj is also 4-dimensional array, with the last index indicating
  ! the field.  The first three fields are the vorticity, the last
  ! three are the current density field.
  real(kind=pr),dimension(:,:,:,:),allocatable :: wj
  ! work is a 3-dimensional array which is used for what? unused? FIXME
  real(kind=pr),dimension(:,:,:),allocatable :: work

  ! Set method information in vars module:
  method="mhd" ! We are doing fluid-structure intergrep actions
  nf=6 ! There are three velocity fields, 3 magnetic fields
  
  ! Initialize MPI, get size and rank
  call MPI_INIT(mpicode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpicode)

  ! FIXME/TODO: initialize time integrals to zero

  if (mpirank == 0) write(*,'(A)') 'Starting MHD3D'

  ! Read input parameters
  call get_command_argument(1,infile) ! infile from command line
  if (mpirank == 0) write(*,'(A)') 'Reading parameters from'//infile
  call get_params(infile)

  ! Initialize FFT
  call fft_initialize

  !  Set up output directory
  call system('mkdir -p fields')
  
  ! FIXME: modify file names, etc, for mhd
  ! Overwrite drag_data file? only if we're not resuming a backup!
  if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then 
    open  (14, file = 'drag_data', status = 'replace')
    close (14)
  endif

  ! Allocate memory:
  ! FIXME: make sure that this is all the right dimension
  allocate(ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
  allocate(ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf))
  allocate(wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf))     
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf,0:1))
  allocate(workvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))) ! unused?

  ! Step forward in time
  if (mpirank == 0)  write(*,'(A)') 'Info: Starting time iterations.'
  call MPI_barrier(MPI_COMM_world,mpicode)
  call time_step(ub,ubk,nlk,wj,work,workvis) ! Actual time-stepping function

  ! Output information on where the algorithm spent the most time.
  if (mpirank == 0) write(*,'(A)') 'Finished computation.'
  
  deallocate(workvis)
  deallocate(ubk)
  deallocate(nlk)
  deallocate(ub)
  deallocate(wj)
  deallocate(work) ! unused?

  call fft_free 
  call MPI_FINALIZE(mpicode)
end program MHD3D
