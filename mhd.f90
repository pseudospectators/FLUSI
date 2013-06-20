program MHD3d
  use mpi_header
  use mhd_vars
  implicit none

  integer                :: mpicode
  character (len=80)     :: infile

  ! Set method information in vars module:
  method="mhd" ! We are doing fluid-structure interactions
  nf=6 ! There are three velocity fields, 3 magnetic fields

  ! Initialize MPI, get size and rank
  call MPI_INIT (mpicode)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)

  ! FIXME/TODO: initialize time integrals to zero

  if (mpirank == 0) write(*,'(A)') 'Starting MHD3D'

  ! Read input parameters
  call get_command_argument(1,infile) ! infile from command line
  if (mpirank == 0) write(*,'(A)') 'Reading parameters from'//infile
  call get_params (infile)

  ! Initialize FFT
  call fft_initialize 

  !  Set up output directory
  call system('mkdir -p fields')
  
  ! Overwrite drag_data file? only if we're not resuming a backup!
  if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then 
    open  (14, file = 'drag_data', status = 'replace')
    close (14)
  endif

  ! FIXME: allocate memory here?

  ! Step forward in time
  if (mpirank == 0) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '*** info: Starting time iterations...'
     write(*,'(A)') '--------------------------------------'
  endif
  call MPI_barrier (MPI_COMM_world, mpicode)
  call time_step() ! Actual time-stepping function

  ! Output information on where the algorithm spent the most time.
  if (mpirank == 0) write(*,'(A)') 'Finished computation.'

  call fft_free 
  call MPI_FINALIZE(mpicode)
end program MHD3D
