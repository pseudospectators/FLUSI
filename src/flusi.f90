program FLUSI
  use mpi
  use vars
  use solid_model
  use module_insects

  implicit none
  integer                :: mpicode
  character (len=strlen) :: infile
  type(diptera)::dummyinsect

  ! Initialize MPI, get size and rank
  call MPI_INIT (mpicode)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)

  if (mpirank==0) root=.true.

  ! this is a fallback:
  ! the reason is that many routines take vector fields as input (either with neq or nd components)
  ! but some postprocessing tools do not set these variables. they are overwritten in almost all
  ! cases by the main programs.
  neq = 3
  nd = 3

  ! get filename of PARAMS file from command line
  call get_command_argument(1,infile)

  if ( index(infile,'.ini') .ne. 0) then
    !-------------------------------------------------------------------------
    ! the file is an *.ini file -> we run a normal simulation
    !-------------------------------------------------------------------------
    call Start_Simulation()

  elseif ((infile=="--postprocess").or.(infile=="-p").or.(infile=="-h").or.(infile=="--help")) then
    !-------------------------------------------------------------------------
    ! the first argument tells us that we're postprocessing
    !-------------------------------------------------------------------------
    call postprocessing()

  elseif ( infile == "--dry-run" ) then
    !-------------------------------------------------------------------------
    ! dry run that only generates and dumps the mask function, without
    ! allocating or computing the fluid.
    !-------------------------------------------------------------------------
    call dry_run()

  elseif ( infile == "--io-test" ) then
    call io_test()

  elseif ( infile=="--solid" .or. infile=="--solid-time-convergence" ) then
    !-------------------------------------------------------------------------
    ! run solid model only
    !-------------------------------------------------------------------------
    method="fsi" ! We are doing fluid-structure interactions
    nf=1 ! We are evolving one field.
    nd=3*nf ! The one field has three components.
    allocate(lin(1)) ! Set up the linear term
    ! read params file from SECOND input parameter
    call get_command_argument(2,infile)
    call get_params(infile,dummyinsect,.true.)
    ! fetch postprocessing argument again:
    call get_command_argument(1,infile)

    if (infile=="--solid") then
      call OnlySolidSimulation()
    elseif (infile=="--solid-time-convergence") then
      call SolidModelConvergenceTest()
    endif

  else
    if (mpirank==0) write(*,*) "nothing to do; the argument " // &
    trim(adjustl(infile)) // " is unkown.."
  endif

  ! normal exit
  if (root) then
    open (15, file='return', status='replace')
    write(15,'(i1)') 0
    close(15)
  endif

  call MPI_FINALIZE(mpicode)
  call exit(0)
end program FLUSI




  subroutine Start_Simulation()
    use mpi
    use vars
    use p3dfft_wrapper
    use solid_model
    use module_insects
    use slicing
    use turbulent_inlet_module
    use penalization ! mask array etc
    implicit none
    real(kind=pr)          :: t1,t2
    real(kind=pr)          :: time,dt0,dt1,memory, mem_field
    integer                :: n0=0,n1=1,it
    integer                :: mpicode
    character (len=strlen)     :: infile
    ! Arrays needed for simulation
    real(kind=pr),dimension(:,:,:,:),allocatable :: explin
    real(kind=pr),dimension(:,:,:,:),allocatable :: u,vort
    complex(kind=pr),dimension(:,:,:,:,:),allocatable :: nlk
    complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
    ! real valued work arrays (there will be "nrw" of them)
    real(kind=pr),dimension(:,:,:,:),allocatable :: work
    ! complex work array, used for sponge and/or passive scalar
    complex(kind=pr),dimension(:,:,:,:),allocatable :: workc
    ! pressure array, with ghost points
    real(kind=pr),dimension(:,:,:),allocatable :: press
    real(kind=pr),dimension(:,:,:,:),allocatable :: scalars
    real(kind=pr),dimension(:,:,:,:,:),allocatable :: scalars_rhs
    ! this is the insect we're using (object oriented)
    type(diptera) :: Insect
    ! this is the solid model beams:
    type(solid), dimension(1:nBeams) :: beams

    ! Set method information in vars module.
    method="fsi" ! We are doing fluid-structure interactions
    nf=1    ! We are evolving one field (that means 1 integrating factor)
    nd=3*nf ! The one field has three components.
    neq=nd  ! number of equations, can be higher than 3 if using passive scalar
    nrw=1   ! number of real valued work arrays
    ncw=1   ! number of complex values work arrays (decide that later)
    nrhs=2  ! number of right-hand side registers

    ! initialize timing variables
    time_fft=0.d0; time_ifft=0.d0; time_vis=0.d0; time_mask=0.d0; time_nlk2=0.d0
    time_vor=0.d0; time_curl=0.d0; time_p=0.d0; time_nlk=0.d0; time_fluid=0.d0
    time_bckp=0.d0; time_save=0.d0; time_total=MPI_wtime(); time_u=0.d0; time_sponge=0.d0
    time_scalar=0.d0
    time_solid=0.d0; time_drag=0.d0; time_surf=0.d0; time_LAPACK=0.d0
    time_hdf5=0.d0; time_integrals=0.d0; time_rhs=0.d0; time_nlk_scalar=0.d0
    tslices=0.d0

    if (root) then
      write(*,'(A)') '--------------------------------------'
      write(*,'(A)') '  FLUSI'
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
    call get_command_argument(1,infile)
    ! read all parameters from that file
    call get_params(infile,Insect,.true.)

    !-----------------------------------------------------------------------------
    ! ghost points. only the "active" FSI part, i.e. with flexible obstacles,
    ! currently needs the ghost point system for interpolating the pressure on the
    ! surface
    !-----------------------------------------------------------------------------
    if (use_solid_model=="yes") then
      if (interp=='linear') ng=1 ! one ghost point
      if (interp=='delta')  ng=3 ! three ghost points
    else
      ! we dont need ghosts when not solving the solid model
      ng=0 ! zero ghost points
    endif
    ! for new passive scalars (with FD discretization) we do also need ghosts
    if (use_passive_scalar==1) then
      ! no less than 3 ghost points:
      ng = maxval( (/ng,3/) )
    endif

    if (root) write(*,'("Set up ng=",i1," ghost points")') ng

    ! we need more memory for RK4:
    if (iTimeMethodFluid=="RK4") nrhs=5
    if (root) write(*,'("Using nrhs=",i1," right hand side registers")') nrhs

    ! the number of entries in the state vector depends on the equation we solve
    if (equation=='navier-stokes') neq = nd ! 3 fields (velocity)
    if (equation=='artificial-compressibility') neq = 4 ! 3 velocities + pressure

    !-----------------------------------------------------------------------------
    ! Initialize FFT (this also defines local array bounds for real and cmplx arrays)
    !-----------------------------------------------------------------------------
    ! Initialize p3dfft
    call fft_initialize()
    ! Setup communicators used for ghost point update
    call setup_cart_groups()

    !-----------------------------------------------------------------------------
    ! Allocate memory:
    !-----------------------------------------------------------------------------
    ! size (in bytes) of one field
    mem_field = dble(nx)*dble(ny)*dble(nz)*8.d0
    ! memory reserved by p3dffft:
    memory = dble(nx)*dble(ny)*dble(nz)*1.6d-5*1000.d0*1000.d0

    ! integrating factors
    allocate(explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
    memory = memory + dble(nf)*mem_field

    ! velocity in Fourier space
    allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
    memory = memory + dble(neq)*mem_field

    ! right hand side of navier-stokes
    allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1))
    memory = memory + dble(nrhs)*dble(neq)*mem_field

    ! velocity in physical space
    allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
    memory = memory + dble(nd)*mem_field

    ! vorticity in physical space (TODO: remove this, add it to work)
    allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
    memory = memory + dble(nd)*mem_field

    ! mask function (defines the geometry)
    allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    mask=0.d0
    memory = memory + mem_field

    ! mask color function (distinguishes between different parts of the mask)
    allocate(mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    mask_color=0
    memory = memory + mem_field/4.d0

    ! solid body velocities
    allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
    us=0.d0
    memory = memory + dble(nd)*mem_field

    ! real valued work array(s)
    ! allocate one work array
    nrw = 1
    if (equation=="artificial-compressibility") nrw = 4
    allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw))
    memory = memory + dble(nrw)*mem_field

    ! pressure array. this is with ghost points for interpolation
    if (use_solid_model=="yes") then
      allocate(press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
      if(mpirank==0) write(*,*) "press array is allocated"
      memory = memory + mem_field
    else
      allocate(press(0:1,0:1,0:1))
    endif

    ! vorticity sponge, work array that is used for sponge and/or passive scalar
    if (iVorticitySponge=="yes") then
      ! three complex work arrays
      ncw = 3
    else
      ! one complex work array, if using scalar
      if (use_passive_scalar==1) ncw = 1
    endif
    if (equation=="artificial-compressibility") ncw = 4
    allocate (workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) )

    ! reserve additional space for scalars?
    if (use_passive_scalar==1) then
      if(mpirank==0) write(*,*) "scalar module is in use: allocate additional memory"
      allocate(scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars))
      allocate(scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1))
      memory = memory + dble((1+nrhs)*n_scalars)*mem_field

      ! this logical "activates" the scalar. if, for example, a NaN in the scalar occurs,
      ! it is set to false and the scalar is skipped, since the fluid can still be okay
      compute_scalar = .true.

    else
      ! HACK HACK On newer intel compilers with array bounds checks (i.e. ifort -CB)
      ! passing unallocated arrays to suborutines causes errors (although these
      ! arrays are of course unused). So for the intel ifort compiler, we allocate
      ! always one scalar.
#ifdef IFORT
      n_scalars = 1
      if(mpirank==0) write(*,*) "scalar module is in use: allocate additional memory (IFORT EXCEPTION)"
      allocate(scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars))
      allocate(scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1))
      memory = memory + dble((1+nrhs)*n_scalars)*mem_field
#endif
    endif

    ! for time averaging
    if (time_avg=="yes") then
      if(mpirank==0) write(*,*) "averaging module is in use: allocate additional memory"

      if (vel_avg=="yes") then
        allocate(uk_avg(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
        memory = memory + dble(3)*mem_field
      endif
      if (ekin_avg=="yes") then
        allocate(e_avg(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
        memory = memory + dble(1)*mem_field
      endif
      if (enstrophy_avg=="yes") then
        allocate(Z_avg(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
        memory = memory + dble(1)*mem_field
      endif
    endif

    ! read in turbulent inlet fields
    if (use_turbulent_inlet=="yes") then
      call init_turbulent_inlet ( )
    endif

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
    ! check if at least FFT works okay
    !-----------------------------------------------------------------------------
    call fft_unit_test(work(:,:,:,1),uk(:,:,:,1))

    !-----------------------------------------------------------------------------
    ! Initial condition
    !-----------------------------------------------------------------------------
    call init_fields(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,explin,work,workc,&
    press,scalars,scalars_rhs,Insect,beams)

    !-----------------------------------------------------------------------------
    ! Initialize time series output files, if not resuming a backup
    !-----------------------------------------------------------------------------
    if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then
      if (time == 0.d0) then
        ! the inicond "infile" reads the time from the hdf5 files and is often used
        ! if runtime backuping failed for some reason. therefore in that case, do not
        ! delete existing time series files.
        call initialize_time_series_files()
      else
        ! do nothing
        write(*,*) "As initial condition set t/=0.d0, I do NOT reset the time series files."
      endif
    endif

    ! initialize runtime control file
    if (mpirank==0) call initialize_runtime_control_file()

    ! Print domain decomposition
    ! call print_domain_decomposition()

    if (use_slicing=="yes") then
      call slice_init(time)
    endif

    !*****************************************************************************
    ! Step forward in time
    !*****************************************************************************
    if (mpirank==0) then
      write(*,'("Penalization parameter C_eta=",es12.4," and K_eta=",es12.4)') eps, &
      sqrt(nu*eps)/dx
    endif

    t1 = MPI_wtime()
    call time_step(time,dt0,dt1,n0,n1,it,u,uk,nlk,vort,work,workc,explin,&
    press,scalars,scalars_rhs,infile,Insect,beams )
    t2 = MPI_wtime() - t1

    !-----------------------------------------------------------------------------
    ! Deallocate memory
    !-----------------------------------------------------------------------------
    ! All processes should reach this point before the arrays are deallocated
    call MPI_barrier(MPI_COMM_WORLD,mpicode)
    deallocate(lin)
    deallocate(explin)
    deallocate(vort,work,workc)
    deallocate(u,uk,nlk)
    deallocate(us)
    deallocate(mask)
    deallocate(mask_color)
    deallocate(ra_table,rb_table)
    if (allocated(press))  deallocate(press)
    if (allocated(uk_old))  deallocate(uk_old)
    if (allocated(nlk_tmp))  deallocate(nlk_tmp)

    if (allocated(uk_avg))  deallocate(uk_avg)
    if (allocated(e_avg)) deallocate(e_avg)
    if (allocated(scalars)) deallocate(scalars)
    if (allocated(scalars_rhs)) deallocate(scalars_rhs)

    if (iMask=="Insect") then
      ! Clean insect (the globally stored arrays for Fourier coeffs etc..)
      call insect_clean(Insect)
    endif

    if (use_slicing=="yes") then
      call slice_free
    endif


    ! release other memory
    call fft_free
    !-------------------------
    ! Show the breakdown of timing information
    !-------------------------
    call show_timings(t2)
  end subroutine Start_Simulation




  ! Output information on where the algorithm spent the most time.
  subroutine show_timings(t2)
    use vars
    use module_helpers
    implicit none
    real (kind=pr) :: t2, t3

    3 format(80("-"))
    8 format(es12.4," (",f5.1,"%) :: ",A)
    t3 = t2
    t2 = mpisum(t2)

    time_fft = mpisum(time_fft)
    time_ifft = mpisum(time_ifft)
    time_mask = mpisum(time_mask)
     time_fluid= mpisum(time_fluid)
     time_integrals =mpisum(time_integrals)
     time_save=mpisum(time_save)
     time_bckp=mpisum(time_bckp)
     tslices=mpisum(tslices)
     time_hdf5=mpisum(time_hdf5)
     time_vis=mpisum(time_vis)
     time_solid=mpisum(time_solid)
     time_rhs=mpisum(time_rhs)
     time_surf=mpisum(time_surf)
     time_p=mpisum(time_p)
     time_nlk2=mpisum(time_nlk2)
     time_nlk_scalar=mpisum(time_nlk_scalar)
     time_scalar=mpisum(time_scalar)
     time_u=mpisum(time_u)
     time_vor=mpisum(time_vor)
     time_sponge=mpisum(time_sponge)
     time_curl=mpisum(time_curl)

if (mpirank/=0) return

    write(*,3)
    write(*,'("*** Timings")')
    write(*,3)
    write(*,'("of the total time ",es12.4,", FLUSI spend ",es12.4," (",f5.1,"%) on FFTS")') &
    t2, (time_fft)+(time_ifft),100.d0*((time_fft)+(time_ifft))/t2
    write(*,3)
    write(*,'("time stepping (top level tasks)")')

    write(*,8) (time_mask), 100.d0*(time_mask)/t2, "create_mask"
    write(*,8) (time_fluid), 100.d0*(time_fluid)/t2, "fluid time stepping"
    write(*,8) (time_integrals), 100.d0*(time_integrals)/t2, "integrals"
    write(*,8) (time_save), 100.d0*(time_save)/t2, "save fields"
    write(*,8) (time_bckp), 100.d0*(time_bckp)/t2, "backuping"
    write(*,8) (tslices), 100.d0*(tslices)/t2, "slicing"
    write(*,3)
    write(*,'("save fields:")')
    write(*,8) (time_hdf5), 100.d0*(time_hdf5)/t2, "hdf5 disk dumping"
    write(*,3)

    write(*,'("Fluid time stepping:")')
    write(*,8) (time_vis),100.d0*(time_vis/t2),"cal_vis"
    write(*,8) (time_rhs),100.d0*(time_rhs/t2),"cal_nlk"
    write(*,8) (time_solid),100.d0*(time_solid/t2),"solid  model"
    write(*,8) (time_surf),100.d0*(time_surf/t2),"surface interpolation"
    write(*,3)

    write(*,'("Fluid right hand side:")')
    write(*,8) (time_nlk2),100.d0*(time_nlk2/t2),"cal_nlk_fsi"
    write(*,8) (time_p),100.d0*(time_p/t2),"pressure"
    write(*,8) (time_nlk_scalar),100.d0*(time_nlk_scalar/t2),"scalar rhs"
    write(*,8) (time_scalar),100.d0*(time_scalar/t2),"passive scalar"
    write(*,3)

    write(*,'("cal_nlk_fsi:")')
    write(*,8) (time_u),100.d0*(time_u)/t2,"velocity"
    write(*,8) (time_vor),100.d0*(time_vor)/t2,"vorticity"
    write(*,8) (time_sponge),100.d0*(time_sponge)/t2,"sponge"
    write(*,8) (time_curl),100.d0*(time_curl)/t2,"nonlinear term"
    write(*,3)
    write(*,'("Integral walltime ",es12.4," (",i7," CPUh)")') t2, nint( t2*dble(mpisize)/3600.d0 )
    write(*,'("Actual   walltime ",es12.4," (",i7," CPUh)")') t3, nint( t3*dble(mpisize)/3600.d0 )
    write(*,3)
    write(*,'(A)') 'Finalizing computation....'
    write(*,3)
  end subroutine show_timings



  subroutine initialize_time_series_files()
    use vars
    implicit none

    ! For insect wing/body forces
    if (iMask=='Insect') then
      open  (14,file='forces.t',status='replace')
      write (14,'(15(A15,1x))') "%          time","Forcex","Forcey","Forcez",&
      "Forcex_unst","Forcey_unst","Forcez_unst",&
      "Momentx","Momenty","Momentz",&
      "Momentx_unst","Momenty_unst","Momentz_unst",&
      "Aero_Power", "Inert power"
      close (14)
      open  (14,file='forces_part1.t',status='replace')
      write (14,'(15(A15,1x))') "%          time","Forcex","Forcey","Forcez",&
      "Forcex_unst","Forcey_unst","Forcez_unst",&
      "Momentx","Momenty","Momentz",&
      "Momentx_unst","Momenty_unst","Momentz_unst",&
      "Aero_Power", "Inert power"
      close (14)
      open  (14,file='forces_part2.t',status='replace')
      write (14,'(15(A15,1x))') "%          time","Forcex","Forcey","Forcez",&
      "Forcex_unst","Forcey_unst","Forcez_unst",&
      "Momentx","Momenty","Momentz",&
      "Momentx_unst","Momenty_unst","Momentz_unst",&
      "Aero_Power", "Inert power"
      close (14)
      open  (14,file='forces_part3.t',status='replace')
      write (14,'(15(A15,1x))') "%          time","Forcex","Forcey","Forcez",&
      "Forcex_unst","Forcey_unst","Forcez_unst",&
      "Momentx","Momenty","Momentz",&
      "Momentx_unst","Momenty_unst","Momentz_unst",&
      "Aero_Power", "Inert power"
      close (14)
      open  (14,file='kinematics.t',status='replace')
      write (14,'(26(A15,1x))') "%          time","xc_body_g","yc_body_g","zc_body_g",&
      "psi","beta","gamma","eta_stroke",&
      "alpha_l","phi_l","theta_l",&
      "alpha_r","phi_r","theta_r",&
      "rot_l_w_x","rot_l_w_y","rot_l_w_z",&
      "rot_r_w_x","rot_r_w_y","rot_r_w_z",&
      "rot_dt_l_w_x","rot_dt_l_w_y","rot_dt_l_w_z",&
      "rot_dt_r_w_x","rot_dt_r_w_y","rot_dt_r_w_z"
      close (14)
      open  (14,file='muscle.t',status='replace')
      close (14)
      open  (14,file='insect_state.t',status='replace')
      close (14)
      ! If this is not an insect
    else
      open  (14,file='forces.t',status='replace')
      write (14,'(13(A15,1x))') "%          time","Forcex","Forcey","Forcez",&
      "Forcex_unst","Forcey_unst","Forcez_unst",&
      "Momentx","Momenty","Momentz",&
      "Momentx_unst","Momenty_unst","Momentz_unst"
      close (14)
    endif

    open  (14,file='divu.t',status='replace')
    write (14,'(3(A15,1x))') "%          time","max_div","max_div_fluid"
    close (14)

    open  (14,file='ekin.t',status='replace')
    write (14,'(2(A15,1x))') "%          time","E_kin_tot"
    close (14)

    open  (14,file='u_residual.t',status='replace')
    write (14,'(7(A15,1x))') "%          time","u_res_color0","u_res_color1",&
    "u_res_color2","u_res_color3","u_res_color4","u_res_color5"
    close (14)

    open  (14,file='energy.t',status='replace')
    write (14,'(21(A15,1x))') "%          time",&
    "E_kin_f","E_kin_x_f","E_kin_y_f","E_kin_z_f",&
    "diss_f","diss_x_f","diss_y_f","diss_z_f",&
    "E_kin_tot","E_kin_tot_x","E_kin_tot_y","E_kin_tot_z",&
    "diss_tot","diss_tot_x","diss_tot_y","diss_tot_z",&
    "Power-penal","Power-penalx","Power-penaly","Power-penalz"
    close (14)

    ! this file contains, time, iteration#, time step and performance
    open  (14,file='timestep.t',status='replace')
    write (14,'(5(A15,1x))') "%            it","time","dt","avg sec/step", "sec/step"
    close (14)

    open  (14,file='dt.t',status='replace')
    write (14,'(5(A15,1x))') "%        time","dt","CFL","viscous", "penalization"
    close (14)

    open  (14,file='meanflow.t',status='replace')
    write (14,'(4(A15,1x))') "%          time","mean_ux","mean_uy","mean_uz"
    close (14)

    open  (14,file='mask_volume.t',status='replace')
    write (14,'(5(A15,1x))') "%          time","volume","mask*usx","mask*usy","mask*usz"
    close (14)

    call init_empty_file('iterations.t')

    open  (14,file='rigidsolidsolver.t',status='replace')
    write (14,'(14(A15,1x))') "%          time","x","y","z","vx","vy","vz",&
    "eps0","eps1","eps2","eps3","rot_body_b_x","rot_body_b_y","rot_body_b_z"
    close (14)


    if (use_passive_scalar==1) call init_empty_file('scalar1.t')
  end subroutine




  subroutine print_domain_decomposition(fname)
    use vars
    use mpi
    implicit none
    character(len=*)::fname
    integer :: mpicode
return
    ! open  (14,file=fname,status='replace')

    if (root) then
      write(14,'(A)') '--------------------------------------'
      write(14,'(A)') '*** Domain decomposition:'
      write(14,'(A)') '--------------------------------------'
    endif
    call MPI_barrier (MPI_COMM_world, mpicode)
    write (14,'("mpirank=",i5," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
    &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
    mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
    call MPI_barrier (MPI_COMM_world, mpicode)
    close(14)
  end subroutine
