
program FLUSI
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  integer                :: mpicode
  integer, dimension (3) :: ifield
  real (kind=pr)         :: t1,t2
  character (len=80)     :: infile

  !---- Initialize MPI, get size and rank, create mpi data types
  call MPI_INIT (mpicode)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)
  

  if (mpirank == 0) then
    write(*,'(A)') '--------------------------------------'
    write(*,'(A)') '		FLUSI'
    write(*,'(A)') '--------------------------------------'
    write(*,'("Running on ",i3," CPUs")') mpisize
    call system('mkdir fields')
  endif

  mpiinteger 	= MPI_INTEGER
  mpireal 	= MPI_DOUBLE_PRECISION
  mpicomplex 	= MPI_DOUBLE_COMPLEX

  !---- Read input parameters
  if (mpirank == 0) then
    write(*,'(A)') '*** info: Reading input data...'
  endif

  call get_command_argument(1, infile)
  call get_params (infile)

  !---- Step forward in time
  if (mpirank == 0) then
    write(*,'(A)') '--------------------------------------'
    write(*,'(A)') '*** info: Starting time iterations...'
    write(*,'(A)') '--------------------------------------'
  endif
  
  t1 = MPI_wtime()
  call time_step ()
  t2 = MPI_wtime() - t1
  
  if (mpirank ==0) then
  write(*,'("$$$ info: total elapsed time time_step=",es12.4, " on ",i2," CPUs")') t2, mpisize
  endif

  !---- Deallocate memory
  if (mpirank == 0) then
    write(*,'(A)') '--------------------------------------'
    write(*,'(A)') '!!! warning: Finalizing computation...'
    write(*,'(A)') '--------------------------------------'
  endif
  
  call fft_free 
  


  !---- Finalize MPI
  call MPI_FINALIZE (mpicode)

end program FLUSI

