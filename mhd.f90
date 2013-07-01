program MHD3D
  use mpi_header
  use mhd_vars
  implicit none

  integer :: mpicode
  character(len=80) :: infile

  ! Arrays needed for simulation  

  ! 4D array (nf 3D arrays) for the implicit time-stepping routines.
  real(kind=pr),dimension(:,:,:,:),allocatable :: explin 

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

  ! work is a 3-dimensional array. FIXME: what is it used in?
  real(kind=pr),dimension(:,:,:),allocatable :: work
  
  ! Initialize MPI, get size and rank
  call MPI_INIT(mpicode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpicode)

  ! Set method information in vars module:
  method="mhd" ! We are doing fluid-structure intergrep actions
  nf=2 ! We are evolving two fields: u and B.
  nd=3*nf ! Each field has three dimensions, for six total.
  allocate(lin(nf)) ! Set up the linear term

  ! FIXME/TODO: initialize time integrals to zero

  if (mpirank == 0) write(*,'(A)') 'Starting MHD3D'

  tab = char(9) ! set horizontal tab character (please kill me now)

  ! Read input parameters

  call get_command_argument(1,infile) ! infile from command line
  if (mpirank == 0) write(*,'(A)') 'Reading parameters from'//infile
  call get_params(infile)

  ! Initialize FFT
  call fft_initialize

  !  Set up output directory
  call system('mkdir -p fields') ! NB: this does not work on turing.
  
  ! Overwrite integral output file? only if we're not resuming a
  ! backup!
  if(mpirank == 0 .and. inicond(1:8).ne."backup::") then
     ! ekvt
     open(14,file='ekvt',status='replace')
10   format(19A)
     write(14,10) "#time",tab,"Ekin",tab,"Ekinx",tab,"Ekiny",tab,"Ekinz"
     close(14)

     ! ebvt
     open(14,file='ebvt',status='replace')
11   format(19A)
     write(14,11) "#time",tab,"Emag",tab,"Emagx",tab,"Emagy",tab,"Emagz"
     close(14)

     ! jvt
     open(14,file='jvt',status='replace')
12   format(15A)
     write(14,12) "#time",tab,"meanjx",tab,"meanjy",tab,"meanjz",tab,&
          "jmax",tab,"jxmax",tab,"jymax",tab,"jzmax"
     close(14) 

     ! dvt
     open(14,file='dvt',status='replace')    
13   format(5A)
     write(14,13) "#time",tab,"divu",tab,"divb"
     close(14)
  endif

  ! Allocate memory:
  allocate(ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd))
  allocate(ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))     
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1))
  allocate(explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))) ! FIXME: unused?

  ! Step forward in time
  if (mpirank == 0)  write(*,'(A)') 'Info: Starting time iterations.'
  call time_step(ub,ubk,nlk,wj,work,explin)

  ! Output information on where the algorithm spent the most time.
  if (mpirank == 0) write(*,'(A)') 'Finished computation.'
  
  deallocate(lin)
  deallocate(explin)
  deallocate(ubk)
  deallocate(nlk)
  deallocate(ub)
  deallocate(wj)
  deallocate(work)

  call fft_free 
  call MPI_FINALIZE(mpicode)
end program MHD3D
