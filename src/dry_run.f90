! Dry run mode - just generate the mask function every tsave time steps
! and write it to disk. Only the minimum of memory is allocated and only
! create_mask is called from here.
! Note this is a very reduced version from start_simulation()
! in flusi.f90
subroutine dry_run()
  use fsi_vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  use penalization ! mask array etc
  use kine ! kinematics from file (Dmitry, 14 Nov 2013)
  implicit none
  real(kind=pr)          :: time,memory,mem_field
  integer                :: it
  character(len=strlen)  :: infile
  character(len=6) :: name
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the solid model beams:
  type(solid), dimension(1:nBeams) :: beams


  ! Set method information in vars module.
  method="fsi" ! We are doing fluid-structure interactions
  nf=1    ! We are evolving one field (that means 1 integrating factor)
  nd=3*nf ! The one field has three components.
  neq=nd  ! number of equations, can be higher than 3 if using passive scalar


  if (root) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '  FLUSI--dry run'
     write(*,'(A)') '--------------------------------------'
     write(*,'("Running on ",i5," CPUs")') mpisize
     write(*,'(A)') '--------------------------------------'
  endif

  !-----------------------------------------------------------------------------
  ! Read input parameters
  !-----------------------------------------------------------------------------
  allocate(lin(nf)) ! Set up the linear term
  if (root) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(2,infile)
  ! read all parameters from that file
  call get_params(infile,Insect)

  !-----------------------------------------------------------------------------
  ! Initialize FFT (this also defines local array bounds for real and cmplx arrays)
  !-----------------------------------------------------------------------------
  call fft_initialize

  !-----------------------------------------------------------------------------
  ! Allocate memory:
  !-----------------------------------------------------------------------------
  ! size (in bytes) of one field
  mem_field = dble(nx)*dble(ny)*dble(nz)*8.d0
  memory = 0.0d0

  ! mask function (defines the geometry)
  allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  memory = memory + mem_field

  ! mask color function (distinguishes between different parts of the mask)
  allocate(mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  memory = memory + mem_field/4.d0

  ! solid body velocities
  allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  memory = memory + dble(nd)*mem_field

  !-----------------------------------------------------------------------------
  ! show memory consumption for information
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,'("Allocated ",i1," real and ",i1," complex work arrays")') nrw,ncw
    write(*,'("FLUSI allocated ",f7.1,"MB (",f5.1,"GB) of memory in total")')&
    memory/(1.0d6),memory/(1.0d9)
    write(*,'("which is ",f7.1,"MB (",f4.1,"GB) per CPU")') &
    memory/(1.0d6)/dble(mpisize),memory/(1.0d9)/dble(mpisize)
    write(*,'(80("-"))')
  endif

  !-----------------------------------------------------------------------------
  ! initalize some insect stuff, if used
  !-----------------------------------------------------------------------------
  ! Load kinematics from file (Dmitry, 14 Nov 2013)
  if (iMask=="Insect") then
    if (Insect%KineFromFile/="no") then
      if (mpirank==0) write(*,*) "Initializing kinematics loader..."
      call load_kine_init(mpirank)
    endif
    ! If required, initialize rigid solid dynamics solver
    if (Insect%BodyMotion=="free_flight") then
      call rigid_solid_init(0.d0,Insect)
    endif
  endif


  if (tsave == 0.d0) then
    if(mpirank==0) write(*,*) "Warning, tsave NOT set assuming 0.05d0!!!"
    tsave = 0.05d0
  endif

  !*****************************************************************************
  ! Step forward in time
  !*****************************************************************************
  time = 0.d0
  it = 0
  do while (time<=tmax)

    call create_mask( time,Insect,beams )

    write(name,'(i6.6)') floor(time*1000.d0)

    if(mpirank==0) then
      write(*,'("Dry run: Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A)') &
      time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
    endif

    call save_field_hdf5(time,'./mask_'//name,mask*eps)
    if (isaveSolidVelocity == 1) then
      call save_field_hdf5(time,'./usx_'//name,us(:,:,:,1))
      call save_field_hdf5(time,'./usy_'//name,us(:,:,:,2))
      call save_field_hdf5(time,'./usz_'//name,us(:,:,:,3))
    endif

    it = it+1
    time = dble(it)*tsave
  enddo


  !-----------------------------------------------------------------------------
  ! Deallocate memory
  !-----------------------------------------------------------------------------
  deallocate(us, mask, mask_color)

  if (iMask=="Insect") then
    ! Clean kinematics (Dmitry, 14 Nov 2013)
    if (Insect%KineFromFile/="no") call load_kine_clean
    ! Clean insect (the globally stored arrays for Fourier coeffs etc..)
    call insect_clean(Insect)
  endif

  ! release other memory
  call fft_free
end subroutine dry_run
