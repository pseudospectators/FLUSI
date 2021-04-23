! Dry run mode - just generate the mask function every tsave time steps
! and write it to disk. Only the minimum of memory is allocated and only
! create_mask is called from here.
! Note this is a very reduced version from start_simulation()
! in flusi.f90
subroutine dry_run()
  use vars
  use p3dfft_wrapper
  use solid_model
  use module_insects
  use turbulent_inlet_module
  use penalization ! mask array etc
  implicit none
  real(kind=pr)          :: time,memory,mem_field
  integer                :: it
  character(len=strlen)  :: infile, mode
  character(len=6) :: name
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the solid model beams:
  type(solid), dimension(1:nBeams) :: beams
  logical :: exists


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
     write(*,'(A)') 'Usage: ./flusi --dry-run PARAMS.ini [MODE]'
     write(*,'(A)') 'where mode can be'
     write(*,'(A)') ' --kinematics   specify insect mask parameters in command line'
     write(*,'(A)') ' --post         reconstruct mask to an existing simulation, reading data from *.t files'
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
  call get_params(infile,Insect,.true.)

  ! is the position of body and wings given by the command line?
  call get_command_argument(3,mode)
  if (mode == "--kinematics") then
    ! for insects:
    ! the flag --kinematics can be used to construct a single mask function with
    ! position and angles (=12 parameters) given by the command line call
    if (root) then
      write(*,*) "parameters are given by command line call"
      write(*,*) "note the mask has NO velocity field! (todo: implement that)"
      write(*,*) "./flusi --dry-run PARAMS.ini --kinematics x y z psi beta gamma &
      &phi_l alpha_l theta_l phi_r alpha_r theta_r eta"
    endif
    ! the following parameters are overwritten (thus not read from ini file)
    ! they trigger that the routines themselves read the kinematics from the
    ! command line arguments.
    Insect%BodyMotion = "command-line"
    Insect%FlappingMotion_left = "command-line-left"
    Insect%FlappingMotion_right = "command-line-right"
    if (Insect%second_wing_pair) then
      Insect%FlappingMotion_left2 = "command-line-left"
      Insect%FlappingMotion_right2 = "command-line-right"
    endif
    ! since the kinematics do not change in time (we have a single snapshot), we
    ! save only one file
    tmax = 0.d0
  endif

  !-----------------------------------------------------------------------------
  ! Initialize domain decomposition, but not FFT (we dont need thats)
  !-----------------------------------------------------------------------------
  call decomposition_initialize()

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

    call insect_init( 0.d0, infile, Insect, .false., "", (/xl,yl,zl/), nu, dx, periodic=periodic)

    ! If required, initialize rigid solid dynamics solver. Note that if the --post flag
    ! is set, the insect state is read from file, so we skip the initialization .
    if (Insect%BodyMotion=="free_flight" .and. mode/="--post") then
      call rigid_solid_init(0.d0, Insect, .false., "")
    endif

    ! tell insect module not to write to kinematics.dry-run.t file, so content of
    ! original file is not touched
    Insect%kinematics_file = "kinematics.dry-run.t"

    if (root) then
      open  (14,file=Insect%kinematics_file, status='replace')
      if (Insect%second_wing_pair) then
        write (14,'(44(A15,1x))') "%          time","xc_body_g","yc_body","zc_body",&
        "psi","beta","gamma","eta_stroke",&
        "alpha_l","phi_l","theta_l",&
        "alpha_r","phi_r","theta_r",&
        "rot_l_x","rot_l_y","rot_l_z",&
        "rot_r_x","rot_r_y","rot_r_z",&
        "rot_dt_l_x","rot_dt_l_y","rot_dt_l_z",&
        "rot_dt_r_x","rot_dt_r_y","rot_dt_r_z",&
        "alpha_l2","phi_l2","theta_l2",&
        "alpha_r2","phi_r2","theta_r2",&
        "rot_l2_x","rot_l2_y","rot_l2_z",&
        "rot_r2_x","rot_r2_y","rot_r2_z",&
        "rot_dt_l2_x","rot_dt_l2_y","rot_dt_l2_z",&
        "rot_dt_r2_x","rot_dt_r2_y","rot_dt_r2_z"
      else
        write (14,'(26(A15,1x))') "%          time","xc_body_g","yc_body","zc_body",&
        "psi","beta","gamma","eta_stroke",&
        "alpha_l","phi_l","theta_l",&
        "alpha_r","phi_r","theta_r",&
        "rot_l_x","rot_l_y","rot_l_z",&
        "rot_r_x","rot_r_y","rot_r_z",&
        "rot_dt_l_x","rot_dt_l_y","rot_dt_l_z",&
        "rot_dt_r_x","rot_dt_r_y","rot_dt_r_z"
      endif
      close (14)
    endif
  endif


  if (tsave == 0.d0) then
    if (mpirank==0) write(*,*) "Warning, tsave NOT set assuming 0.05d0!!!"
    tsave = 0.05d0
  endif

  ! read in turbulent inlet fields
  if (use_turbulent_inlet=="yes") then
    call init_turbulent_inlet ( )
  endif

  !*****************************************************************************
  ! Step forward in time
  !*****************************************************************************
  ! tstart is read in parameter file and default value 0.0
  ! this way you can have better control over dry-runs
  time = tstart
  it = 0
  do while (time<=tmax)

    ! if the motion is free_flight, dry run can still be used e.g. to integrate
    ! a constant velocity which is set in the ini file. of course, a dry run cannot
    ! take the FSI coupling into account.
    ! sometimes, one wants to run a dry-run as postprocessing to an existing simulation
    ! (to save HDD space and erase the mask), and if that run is free_flight than the code
    ! should read the insect state from the rigidsolidsolver.t file (instead of simply
    ! integrating the rigid solid time stepper.)
    if (Insect%BodyMotion == "free_flight") then
      inquire( file='rigidsolidsolver.t', exist=exists )
      if (exists) then
        ! in this case, we find the rigidsolidsolver file, and we assume the user
        ! wants to reconstruct the mask from that file.
        if (root) write(*,*) "DRY-RUN found rigidsolidsolver.t file and use that for mask generation!"
        call read_insect_STATE_from_file(time, Insect, 'rigidsolidsolver.t', verbose=.true.)
      else
        ! use rigid solid solver to integrate the body motion state; this is useful only
        ! if a constant velocity is set
        call rigid_solid_time_step(time, tsave, tsave, it, Insect, (/0.0_pr, 0.0_pr, 0.0_pr/), (/0.0_pr, 0.0_pr, 0.0_pr/))
      endif
    endif


    ! create the mask
    call create_mask(time, Insect, beams)

    if (iMask=="Insect".and.modulo(it,itdrag)==0) then
        call write_kinematics_file( time, Insect )
    endif

    ! Save data
    write(name,'(i6.6)') floor(time*1000.d0)

    if(mpirank==0) then
      write(*,'("Dry run: Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A)') &
      time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
    endif

    call save_field_hdf5(time, 'mask_'//name, mask)
    if (isaveSolidVelocity == 1) then
      call save_field_hdf5(time,'usx_'//name,us(:,:,:,1))
      call save_field_hdf5(time,'usy_'//name,us(:,:,:,2))
      call save_field_hdf5(time,'usz_'//name,us(:,:,:,3))
    endif

    it = it+1
    time = tstart + dble(it)*tsave
  enddo


  !-----------------------------------------------------------------------------
  ! Deallocate memory
  !-----------------------------------------------------------------------------
  deallocate(us, mask, mask_color)

  if (iMask=="Insect") then
    ! Clean insect
    call insect_clean(Insect)
  endif

  ! release other memory
  call fft_free
end subroutine dry_run

! The same dry_run subroutine as above but modified for flexible wing model
subroutine dry_run_flexible_wing()
  use vars
  use p3dfft_wrapper
  use solid_model
  use flexible_model
  use penalization ! mask array etc
  use stl_file_reader
  use module_insects
  !use helpers
  !use module_ini_files_parser_mpi
  implicit none
  real(kind=pr)          :: t1,t2
  real(kind=pr)          :: time,memory,mem_field
  integer                :: it, i, j
  character(len=strlen)  :: infile, mode
  character(len=6) :: name
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the solid model beams:
  type(solid), dimension(1:nBeams) :: beams
  ! this is the wings we're using (object oriented)
  type(flexible_wing),dimension(1:nWings) :: Wings
  logical :: exists


  ! Set method information in vars module.
  method="fsi" ! We are doing PSEUDO fluid-structure interactions
  nf=1    ! We are evolving one field (that means 1 integrating factor)
  nd=3*nf ! The one field has three components.
  neq=nd  ! number of equations, can be higher than 3 if using passive scalar

  if (root) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '  FLUSI--dry run flexible wing'
     write(*,'(A)') '--------------------------------------'
     write(*,'("Running on ",i5," CPUs")') mpisize
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') 'Usage: ./flusi --dry-run-flexible-wing PARAMS.ini [MODE]'
     write(*,'(A)') 'where mode can be'
     write(*,'(A)') ' --fixed   wings are fixed at the wingbases'
     write(*,'(A)') ' --kinematics   specify wings mask parameters in command line'
     write(*,'(A)') '--------------------------------------'
  endif


  !-----------------------------------------------------------------------------
  ! Read input parameters and mesh data
  !-----------------------------------------------------------------------------
  allocate(lin(nf)) ! Set up the linear term
  if (root) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(2,infile)

  write(*,*) infile

  ! read all parameters from that file
  call get_params(infile,Insect,.true.)

  ! is the position of body and wings given by the command line?
  call get_command_argument(3,mode)

  !-----------------------------------------------------------------------------
  ! Checking
  !-----------------------------------------------------------------------------
  if (iMask/="Flexible_wing") then
    call abort(476659, "dry-run-flexible-wing is used only for flexible wing model, &
    change iMask into Flexible_wing to continue")
  endif



  !-----------------------------------------------------------------------------
  ! Initialize domain decomposition, but not FFT (we dont need thats)
  !-----------------------------------------------------------------------------
  call decomposition_initialize()

  !-----------------------------------------------------------------------------
  ! Allocate memory:
  !-----------------------------------------------------------------------------
  ! size (in bytes) of one field
  mem_field = dble(nx)*dble(ny)*dble(nz)*8.d0
  memory = 0.0d0

  ! mask function (defines the geometry)
  allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  memory = memory + mem_field

  ! mask function (defines the geometry)
  !allocate(unsigned_distance(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  !memory = memory + mem_field

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
  ! initalize wings
  !-----------------------------------------------------------------------------
  call init_wings( infile, Wings, dx)

  if (tsave == 0.d0) then
    if(mpirank==0) write(*,*) "Warning, tsave NOT set assuming 0.05d0!!!"
    tsave = 0.05d0
  endif

  !*****************************************************************************
  ! Step forward in time
  !*****************************************************************************
  ! tstart is read in parameter file and default value 0.0
  ! this way you can have better control over dry-runs
  time = tstart
  it = 0

  ! create the startup mask function
  call create_mask_fsi(time,Insect,beams,wings)

  ! Save data
  write(name,'(i6.6)') floor(time*1000.d0)

  call save_field_hdf5(time,'mask_'//name,mask)
  !call save_field_hdf5(time,'unsigned_distance_'//name,unsigned_distance)
  if (isaveSolidVelocity == 1) then
      call save_field_hdf5(time,'usx_'//name,us(:,:,:,1))
      call save_field_hdf5(time,'usy_'//name,us(:,:,:,2))
      call save_field_hdf5(time,'usz_'//name,us(:,:,:,3))
  endif

  do while (time<tmax)
      t1 = MPI_wtime()

      do i=1,nWings
          ! call flexible_wing_motions ( time, wings(i) )

          !
          call flexible_solid_time_step(time, tsave, tsave, it, wings(i))

      enddo
      ! create the mask
      call create_mask_fsi(time,Insect,beams,wings)


      it = it+1
      time = tstart + dble(it)*tsave

      ! Save data
      write(name,'(i6.6)') floor(time*1000.d0)

      if(mpirank==0 .and. modulo(it,20)==0) then
          write(*,'("Dry run flexible wing: Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A)') &
          time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
      endif

      call save_field_hdf5(time,'mask_'//name,mask)
      !call save_field_hdf5(time,'unsigned_distance_'//name,unsigned_distance)
      if (isaveSolidVelocity == 1) then
          call save_field_hdf5(time,'usx_'//name,us(:,:,:,1))
          call save_field_hdf5(time,'usy_'//name,us(:,:,:,2))
          call save_field_hdf5(time,'usz_'//name,us(:,:,:,3))
      endif

      if (root) then
          call SaveWingData( time, wings )
      endif

      call toc("DryRun time step", MPI_wtime()-t1)
  enddo


  !-----------------------------------------------------------------------------
  ! Deallocate memory
  !-----------------------------------------------------------------------------
  deallocate(us, mask, mask_color)

  ! release other memory
  call fft_free
end subroutine dry_run_flexible_wing
