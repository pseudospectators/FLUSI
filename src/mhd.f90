program mhd
  use mpi
  use p3dfft_wrapper
  use vars
  use module_insects !TODO: MAKE MHD INDEPENDENT OF THIS
  use solid_model!TODO: MAKE MHD INDEPENDENT OF THIS
  use flexible_model!TODO: MAKE MHD INDEPENDENT OF THIS
  use penalization ! mask array etc
  implicit none

  integer :: mpicode
  character(len=strlen) :: infile
  real(kind=pr) :: time,dt0,dt1 ! FIXME: move to vars.
  integer :: n0=0,n1=1
  integer :: it
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

  real(kind=pr),dimension(:,:,:,:),allocatable :: work
  ! the following arrays are dummies and never allocated:
  real(kind=pr),dimension(:,:,:),allocatable :: press
  real(kind=pr),dimension(:,:,:,:),allocatable :: scalars
  real(kind=pr),dimension(:,:,:,:,:),allocatable :: scalars_rhs

  ! complex work array, currently unused in the MHD case (not allocated)
  complex(kind=pr),dimension(:,:,:,:),allocatable :: workc

  ! this is a hack and will be removed later:
  type(diptera) :: dummy_insect
  type(flexible_wing),dimension(1:nWings) :: dummy_wings
  type(solid), dimension(1:nBeams) :: dummy_beams

  ! Initialize MPI, get size and rank
  call MPI_INIT(mpicode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpicode)
  if (mpirank==0) root=.true.

  ! MHD does not use ghost points
  ng=0
  ! MHD does not use scalars
  n_scalars=0

  ! Set method information in vars module:
  method="mhd" ! We are doing fluid-structure intergrep actions
  nf=2 ! We are evolving two fields: u and B.
  nd=3*nf ! Each field has three dimensions, for six total.
  neq=nd ! the vector of unknowns has 6 dimensions (ux,uy,uz,bx,by,bz)
  nrhs=2 ! number of right hand side registers

  nrw = 1 ! number of real work arrays in work
  ncw = 0 ! number of complex work arrays in workc

  allocate(lin(nf)) ! Set up the linear term

  ! this helps displaying the walltime (in the global var time_total)
  time_total=MPI_wtime()

  if (mpirank == 0) write(*,'(A)') 'Starting MHD3D'

  tab = char(9) ! set horizontal tab character (please kill me now)

  ! Read input parameters

  call get_command_argument(1,infile) ! infile from command line
  if (mpirank == 0) write(*,'(A)') 'Reading parameters from'//infile
  call get_params(infile,dummy_insect,.true.)

  call fft_initialize

  ! Overwrite integral output file? only if we're not resuming a
  ! backup!
  if(mpirank == 0 .and. inicond(1:8).ne."backup::") then
     ! Formatting things for FORTRAN text output, including tabs.
     ! Please, please kill me now.
5    format(5A) ! 5 outputs, including tabs
9    format(9A) ! 9 outputs, including tabs
15   format(15A) ! 15 outputs, including tabs

     ! ek.t: 9 outputs, including tabs
     open(14,file='ek.t',status='replace')
     write(14,9) "#time",tab,"Ekin",tab,"Ekinx",tab,"Ekiny",tab,"Ekinz"
     close(14)

     ! eb.t: 9 outputs, including tabs
     open(14,file='eb.t',status='replace')
     write(14,9) "#time",tab,"Emag",tab,"Emagx",tab,"Emagy",tab,"Emagz"
     close(14)

     ! j.t: 15 outputs, including tabs
     open(14,file='j.t',status='replace')
     write(14,15) "#time",tab,"meanjx",tab,"meanjy",tab,"meanjz",tab,&
          "jmax",tab,"jxmax",tab,"jymax",tab,"jzmax"
     close(14)

     ! diss.t: 5 outputs, including tabs
     open(14,file='diss.t',status='replace')
     write(14,5) "#time",tab,"disskin",tab,"dissmag"
     close(14)

     ! d.t: 5 outputs, including tabs
     open(14,file='d.t',status='replace')
     write(14,5) "#time",tab,"divu",tab,"divb"
     close(14)

     ! this file contains, time, iteration#, time step and performance
     open(14,file='timestep.t',status='replace')
     close(14)
  endif

  ! initialize runtime control file
  if (mpirank==0) call initialize_runtime_control_file()

  ! Allocate memory:
  allocate(ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  allocate(ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1))
  allocate(explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw))
  allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))

  ! HACK HACK On newer intel compilers with array bounds checks (i.e. ifort -CB)
  ! passing unallocated arrays to suborutines causes errors (although these
  ! arrays are of course unused). So for the intel ifort compiler, we allocate
  ! always one scalar.
#ifdef IFORT
  n_scalars = 1
  if(mpirank==0) write(*,*) "scalar module is in use: allocate additional memory (IFORT EXCEPTION)"
  allocate(scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars))
  allocate(scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1))
#endif

  ! Check if at least FFT works okay
  call fft_unit_test(ub(:,:,:,1),ubk(:,:,:,1))

  time=0.0
  it=0
  dt0=1.0d0
  dt1=2.0d0

  ! Initialize vorticity or read values from a backup file
  if (mpirank == 0) write(*,*) "Set up initial conditions:"
  call init_fields(time,it,dt0,dt1,n0,n1,ub,ubk,nlk,wj,explin,work,workc,press,&
       scalars,scalars_rhs,dummy_insect,dummy_beams,dummy_wings)

  if (mpirank == 0) write(*,*) "Create mask variables:"
  call create_mask_mhd
  if (mpirank == 0) write(*,*) "update_us:"
  call update_us(ub)

  if (mpirank == 0) write(*,*) "Start time-stepping:"
  call time_step(time,dt0,dt1,n0,n1,it,ub,ubk,nlk,wj,work,workc,explin,&
       press,scalars,scalars_rhs,infile,dummy_insect,dummy_beams,dummy_wings)
  if (mpirank == 0) write(*,'(A)') 'Finished computation.'

  deallocate(ubk)
  deallocate(ub)
  deallocate(wj)
  deallocate(nlk)
  deallocate(explin)
  deallocate(mask)
  deallocate(us)
  deallocate(lin)
  deallocate(work)
  deallocate(workc)
  deallocate(ra_table,rb_table)

  call fft_free
  call MPI_FINALIZE(mpicode)
end program mhd
