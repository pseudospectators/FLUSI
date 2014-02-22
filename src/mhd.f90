program mhd
  use mpi
  use p3dfft_wrapper
  use mhd_vars
  implicit none

  integer :: mpicode
  character(len=80) :: infile
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

  call set_fluid_solid_communicators()
  
  ! FIXME/TODO: initialize time integrals to zero

  if (mpirank == 0) write(*,'(A)') 'Starting MHD3D'

  tab = char(9) ! set horizontal tab character (please kill me now)

  ! Read input parameters

  call get_command_argument(1,infile) ! infile from command line
  if (mpirank == 0) write(*,'(A)') 'Reading parameters from'//infile
  call get_params(infile)
  
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

  ! Allocate memory:
  allocate(ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd))
!  call alloccomplexnd(ubk)
  call allocrealnd(ub)
  call allocrealnd(wj)
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1))
  allocate(explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
  call allocreal(work)
  call allocreal(mask)
  call allocrealnd(us)

  ! Check if at least FFT works okay
  call fft_unit_test(work,ubk(:,:,:,1))

  time=0.0
  it=0
  dt0=1.0d0
  dt1=2.0d0

  ! Initialize vorticity or read values from a backup file
  if (mpirank == 0) write(*,*) "Set up initial conditions:"
  call init_fields(n1,time,it,dt0,dt1,ubk,nlk,wj,explin)
  n0=1 - n1 !important to do this now in case we're retaking a backp
  
  if (mpirank == 0) write(*,*) "Create mask variables:"
  call create_mask_mhd
  if (mpirank == 0) write(*,*) "update_us:"
  call update_us(ub)
  
  if (mpirank == 0) write(*,*) "Start time-stepping:"
  call time_step(ub,ubk,nlk,wj,work,explin,infile,time,dt0,dt1,n0,n1,it)

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

  call fft_free 
  call MPI_FINALIZE(mpicode)
end program mhd
