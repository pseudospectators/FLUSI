program FLUSI
  use vars
  use solid_model
  use insect_module
  implicit none
  integer                :: mpicode
  character (len=strlen) :: infile
  type(diptera)::dummyinsect

  ! Initialize MPI, get size and rank
  call MPI_INIT (mpicode)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode) 
  
  if (mpirank==0) root=.true.
  
  ! get filename of PARAMS file from command line
  call get_command_argument(1,infile)
  
  if ( index(infile,'.ini') .ne. 0) then  
      !-------------------------------------------------------------------------
      ! the file is an *.ini file -> we run a normal simulation 
      !-------------------------------------------------------------------------
      call Start_Simulation()    
      
  elseif ( infile == "--postprocess" .or. infile == "--post") then 
      !-------------------------------------------------------------------------
      ! the first argument tells us that we're postprocessing 
      !-------------------------------------------------------------------------
      call postprocessing()
      
  elseif ( infile == "--solid" ) then
      !-------------------------------------------------------------------------
      ! run solid model only
      !-------------------------------------------------------------------------
      method="fsi" ! We are doing fluid-structure interactions
      call get_command_argument(2,infile)
      call get_params(infile,dummyinsect)
      call OnlySolidSimulation()
      
  else
      if (mpirank==0) write(*,*) "nothing to do..."      
  endif

  
  call MPI_FINALIZE(mpicode)
  call exit(0)
end program FLUSI


subroutine Start_Simulation()
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  use kine ! kinematics from file (Dmitry, 14 Nov 2013)
  implicit none
  real(kind=pr)          :: t1,t2
  real(kind=pr)          :: memory, mem_field
  type(timetype) :: time
  character (len=strlen) :: infile
  
  ! mask color function
  integer(kind=2),dimension (:,:,:),allocatable,save :: mask_color
  ! mask containing the obstacle
  real(kind=pr),dimension(:,:,:),allocatable :: mask
  ! solution vector u
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  ! velocity-field inside solid 
  real(kind=pr),dimension(:,:,:,:),allocatable :: us
  ! work arrays
  real(kind=pr),dimension(:,:,:,:),allocatable :: work
  ! the right hand side
  real(kind=pr),dimension(:,:,:,:,:),allocatable :: nlk
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the solid model beams:
  type(solid), dimension(1:nBeams) :: beams
  
  ! Set method information in vars module.
  t1 = MPI_wtime()
  neq=4  ! number of equations, can be higher than 3 if using passive scalar
  nrw=1  ! number of real valued work arrays
  
  ! initialize timing variables
  time_mask=0.d0
  time_fluid=0.d0
  time_bckp=0.d0; time_save=0.d0; time_total=MPI_wtime()
  time_insect_body=0.d0; 
  time_insect_wings=0.d0; time_insect_vel=0.d0; time_scalar=0.d0
  time_solid=0.d0; time_drag=0.d0; time_surf=0.d0; time_LAPACK=0.d0
  time_hdf5=0.d0; time_integrals=0.d0; time_rhs=0.d0; time_nlk_scalar=0.d0
  
  if (root) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '  FLUSI-mistral'
     write(*,'(A)') '--------------------------------------'
     write(*,'("Running on ",i5," CPUs")') mpisize
     write(*,'(A)') '--------------------------------------'
  endif

  
  !-----------------------------------------------------------------------------
  ! Read input parameters
  !-----------------------------------------------------------------------------
  if (root) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(1,infile)
  ! read all parameters from that file
  call get_params(infile,Insect)
  
  !-----------------------------------------------------------------------------
  ! decide about memory
  !-----------------------------------------------------------------------------
  if (iTimeMethodFluid=="RK2") then
    nrhs=2  ! number of registers for right hand side vectors
  elseif (iTimeMethodFluid=="RK4") then
    nrhs=5  ! number of registers for right hand side vectors
  elseif (iTimeMethodFluid=="AB2") then
    nrhs=2  ! number of registers for right hand side vectors
    time%n0=1
    time%n1=2
  elseif (iTimeMethodFluid=="FSI_RK2_semiimplicit") then
    nrhs = 2
  elseif (iTimeMethodFluid=="FSI_RK4_semiimplicit") then
    nrhs = 5
  else
    if (root) write(*,*) "flusi.f90 :: error: iTimeMethodFluid is unknown"
    if (root) write(*,*) iTimeMethodFluid
    call abort()
  endif
  
  !-----------------------------------------------------------------------------
  ! ghost points. the number of ghosts we need depends on the Discretization
  ! and the interpolation method for FSI pressure
  !-----------------------------------------------------------------------------
  ng = 0
  if (method=="centered_2nd") ng = 1
  if (method=="centered_4th") ng = 3
  if (use_solid_model=="yes") ng = 3
  if (root) write(*,'("Set up ng=",i1," ghost points")') ng
  if (root) write(*,'("Discretization: ",A)') trim(method)
  
  !-----------------------------------------------------------------------------
  ! Initialize FFT (this also defines local array bounds for arrays)
  !-----------------------------------------------------------------------------  
  call fft_initialize 
  
  !-----------------------------------------------------------------------------
  ! Initialize time series output files, if not resuming a backup
  !-----------------------------------------------------------------------------
  if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then 
    call initialize_time_series_files()
  endif
  
  ! initialize runtime control file
  if (mpirank==0) call initialize_runtime_control_file()
 
  ! Print domain decomposition
  call print_domain_decomposition()

  !-----------------------------------------------------------------------------
  ! Allocate memory:
  !-----------------------------------------------------------------------------
!   ! reserve additional space for scalars?
!   if (use_passive_scalar==1) then
!     ! add n_scalars to number of equations. Note the difference between neq and nd
!     neq = neq + n_scalars 
!     ! this logical "activates" the scalar. if, for example, a NaN in the scalar occurs,
!     ! it is set to false and the scalar is skipped, since the fluid can still be okay
!     compute_scalar = .true.
!   endif

  ! size (in bytes) of one field  
  mem_field = dble(nx)*dble(ny)*dble(nz)*8.d0
  memory = 0.d0
  
  ! right hand side of navier-stokes (possibly with passive scalar)
  allocate(nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs))
  memory = memory + dble(neq*nrhs)*mem_field
  
  ! velocity in physical space (WITHOUT passive scalar)
  allocate(u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq))
  memory = memory + dble(neq)*mem_field
  
  ! mask function (defines the geometry)
  allocate(mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))  
  memory = memory + mem_field
  
  ! mask color function (distinguishes between different parts of the mask)
  allocate(mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
  memory = memory + mem_field/4.d0
  
  ! solid body velocities
  allocate(us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)) 
  memory = memory + dble(neq)*mem_field
  
  ! allocate one work array
  allocate(work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw))
  memory = memory + dble(nrw)*mem_field

  !-----------------------------------------------------------------------------
  ! show memory consumption for information
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,'("Allocated ",i1," real and ",i1," complex work arrays")') nrw,0
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
!   if (iMask=="Insect") then
!     if (Insect%KineFromFile/="no") then
!       if (mpirank==0) write(*,*) "Initializing kinematics loader..."
!       call load_kine_init(mpirank)
!     endif
!     ! If required, initialize rigid solid dynamics solver
!     ! and set idynamics flag on or off
!     call rigid_solid_init(SolidDyn%idynamics,Insect)
!   endif
  
  !-----------------------------------------------------------------------------
  ! Initial condition
  !-----------------------------------------------------------------------------
  call init_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  
  
  !*****************************************************************************
  ! Step forward in time
  !*****************************************************************************
  call time_step(time,u,nlk,work,mask,mask_color,us,Insect,beams,infile)
  
  !-----------------------------------------------------------------------------
  ! Deallocate memory
  !-----------------------------------------------------------------------------
  deallocate(work)
  deallocate(u,nlk)
  deallocate(us)
  deallocate(mask)
  deallocate(mask_color)
  deallocate(ra_table,rb_table)
  
  
!   if (iMask=="Insect") then
!     ! Clean kinematics (Dmitry, 14 Nov 2013)
!     if (Insect%KineFromFile/="no") call load_kine_clean
!     ! Clean insect (the globally stored arrays for Fourier coeffs etc..)
!     call insect_clean(Insect)
!   endif
  
  ! write empty success file
  if (root) call init_empty_file("success")
  
  ! release other memory
  call fft_free 
  
  !-------------------------
  ! Show the breakdown of timing information
  !-------------------------
  t2 = MPI_wtime() - t1
  if (root) call show_timings(t2)
end subroutine Start_Simulation




! Output information on where the algorithm spent the most time.
subroutine show_timings(t2)
  use vars
  implicit none
  real (kind=pr) :: t2

3 format(80("-"))
8 format(es12.4," (",f5.1,"%) :: ",A)


  write(*,3)
  write(*,'("*** Timings")')
  write(*,3)
  write(*,'("time stepping (top level tasks)")')
  
  write(*,8) time_fluid, 100.d0*time_fluid/t2, "fluid time stepping"
  write(*,8) time_integrals, 100.d0*time_integrals/t2, "integrals"
  write(*,8) time_save, 100.d0*time_save/t2, "save fields"
  write(*,8) time_bckp, 100.d0*time_bckp/t2, "backuping"
  write(*,3)
  write(*,'("Create Mask:")')
  write(*,8) time_insect_body, 100.d0*time_insect_body/t2, "insect::body"
  write(*,8) time_insect_wings,100.d0*time_insect_wings/t2,"insect::wings"
  write(*,8) time_insect_vel,100.d0*time_insect_vel/t2,"insect::roration"
  write(*,3)
  write(*,'("save fields:")')
  write(*,8) time_hdf5, 100.d0*time_hdf5/t2, "hdf5 disk dumping"
  write(*,3)
  
  write(*,'("Fluid time stepping:")')
  write(*,8) time_mask, 100.d0*time_mask/t2, "create_mask"
  write(*,8) time_rhs,100.d0*time_rhs/t2,"cal_nlk"
  write(*,8) time_solid,100.d0*time_solid/t2,"solid  model"
  write(*,8) time_surf,100.d0*time_surf/t2,"surface interpolation"
  write(*,3)
  
  write(*,'("Total walltime ",es12.4," (",i7," CPUh)")') t2, nint( t2*dble(mpisize)/3600.d0 )
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
    write (14,'(26(A15,1x))') "%          time","xc_body","yc_body","zc_body",&
                      "psi","beta","gamma","eta_stroke",&
                      "alpha_l","phi_l","theta_l",&
                      "alpha_r","phi_r","theta_r",&
                      "rot_l_x","rot_l_y","rot_l_z",&
                      "rot_r_x","rot_r_y","rot_r_z",&
                      "rot_dt_l_x","rot_dt_l_y","rot_dt_l_z",&
                      "rot_dt_r_x","rot_dt_r_y","rot_dt_r_z"
                      
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
  
  open  (14,file='energy.t',status='replace')
  write (14,'(18(A15,1x))') "%          time",&
                     "E_kin_f","E_kin_x_f","E_kin_y_f","E_kin_z_f",&
                     "diss_f","diss_x_f","diss_y_f","diss_z_f",&
                     "E_kin_tot","E_kin_tot_x","E_kin_tot_y","E_kin_tot_z",&
                     "diss_tot","diss_tot_x","diss_tot_y","diss_tot_z",&
                     "flux_penal"
                      
  close (14)
  
  ! this file contains, time, iteration#, time step and performance
  open  (14,file='timestep.t',status='replace')
  write (14,'(4(A15,1x))') "%            it","time","dt","avg sec/step", "sec/step"
  close (14)    
  
  open  (14,file='meanflow.t',status='replace')
  write (14,'(4(A15,1x))') "%          time","mean_ux","mean_uy","mean_uz"
  close (14)  

  call init_empty_file('iterations.t')
  call init_empty_file('mask_volume.t')
  if (use_passive_scalar==1) call init_empty_file('scalar.t')
end subroutine




subroutine print_domain_decomposition()
  use vars
  use mpi
  implicit none
  integer :: mpicode
  
  if (root) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '*** Domain decomposition:'
     write(*,'(A)') '--------------------------------------'
  endif
  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("rank=",i5," (",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
       mpirank, ga(1),gb(1), ga(2),gb(2),ga(3),gb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)
end subroutine 
