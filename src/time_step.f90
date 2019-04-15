subroutine time_step(time,dt0,dt1,n0,n1,it,u,uk,nlk,vort,work,workc,explin,&
  press,scalars,scalars_rhs,params_file,Insect,beams,wings)
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use flexible_model
  use module_insects
  use slicing

  implicit none

  integer :: inter
  integer,intent(inout) :: n0,n1,it
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  integer :: it_start
  real(kind=pr),intent(inout) :: time,dt0,dt1
  real(kind=pr) :: t1,t2,t3,t4
  character(len=strlen)  :: command ! for runtime control
  character(len=strlen),intent(in)  :: params_file ! for runtime control
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  ! pressure array. this is with ghost points for interpolation
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)

  type(flexible_wing),dimension(1:nWings), intent(inout) :: Wings
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  logical :: continue_timestepping
  t1 = MPI_wtime()

  ! first backup is in "truntime" hours
  truntimenext = truntime
  continue_timestepping = .true.
  it_start = it

  ! After init, output integral quantities. (note we can overwrite only
  ! nlk(:,:,:,:,n0) when retaking a backup) (if not resuming a backup)
  if (index(inicond,'backup::')==0) then
    if (root) write(*,*) "Initial output of integral quantities...."
    call write_integrals(time,uk,u,vort,nlk(:,:,:,:,n0),work,scalars,Insect,beams,Wings)
  endif


  !-----------------------------------------------------------------------------
  ! Loop over time steps
  !-----------------------------------------------------------------------------
  if (root) write(*,*) "Start time-stepping...."
  do while ((time<tmax) .and. (it<=nt) .and. (continue_timestepping))
    t4 = MPI_wtime()
    dt0 = dt1

    !---------------------------------------------------------------------------
    ! advance fluid/B-field in time
    !---------------------------------------------------------------------------
    ! note: the array "vort" is a real work array and has neither input nor
    ! output values after fluid time stepper
    call fluidtimestep(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,workc,explin,&
         press,scalars,scalars_rhs,Insect,beams,Wings)

    !---------------------------------------------------------------------------
    ! time step done: advance iteration + time
    !---------------------------------------------------------------------------
    inter=n1 ; n1=n0 ; n0=inter
    ! Advance in time so that uk contains the evolved field at time 'time+dt1'
    time = time + dt1
    it = it + 1

    !---------------------------------------------------------------------------
    ! Output of INTEGRALS after every tintegral time units or itdrag time steps
    !---------------------------------------------------------------------------
    if (time_for_output(time, dt1, it, tintegral, itdrag, tmax, 0.d0)) then
      call write_integrals(time,uk,u,vort,nlk(:,:,:,:,n0),work,scalars,Insect,beams,Wings)
    endif

    !---------------------------------------------------------------------------
    ! Save solid model data, if using it
    !---------------------------------------------------------------------------
    if (use_solid_model=="yes" .and. root .and. modulo(it,itbeam)==0) then
      call SaveBeamData( time, beams )
    endif

    !---------------------------------------------------------------------------
    ! Save flexible wing model data, if using it
    !---------------------------------------------------------------------------
    if (use_flexible_wing_model=="yes" .and. root) then
      call SaveWingData( time, wings )
    endif

    !---------------------------------------------------------------------------
    ! Output FIELDS DATA (after tsave time units, but not before tsave_first)
    !---------------------------------------------------------------------------
    if (time_for_output(time, dt1, it, tsave, 99999999, tmax, tsave_first)) then
      ! Note: we can safely delete nlk(:,:,:,1:neq,n0). for RK2 it never matters,
      ! and for AB2 this is the one to be overwritten in the next step. This frees
      ! 3 complex arrays, which are then used in Dump_Runtime_Backup.
      call save_fields(time,it,uk,u,vort,nlk(:,:,:,:,n0),work,workc,press,scalars,scalars_rhs,Insect,beams,wings)
      call are_we_there_yet(time,t1,dt1)
    endif

    ! Backup if that's specified in the PARAMS.ini file. We try to do
    ! backups every "truntime" hours (precise to one time step)
    if (idobackup==1 .and. truntimenext<(MPI_wtime()-time_total)/3600.d0) then
      call dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,&
      work(:,:,:,1),scalars,scalars_rhs,Insect,beams,wings)
      truntimenext = truntimenext+truntime
    endif

    !---------------------------------------------------------------------------
    ! save slices to hard disk
    !---------------------------------------------------------------------------
    if ((method=="fsi").and.(use_slicing=="yes")) then
      if (time_for_output(time, dt1, it, tslice, itslice, 9.9d9, tslice_first)) then
        call save_slices( time, u )
      endif
    endif

    !---------------------------------------------------------------------------
    ! Output how much time remains
    !---------------------------------------------------------------------------
    if ((modulo(it,10)==0).or.(it==20)) then
      call are_we_there_yet(time,t1,dt1)
    endif

    !-------------------------------------------------
    ! escape from loop if walltime is about to be exceeded
    !-------------------------------------------------
    if (wtimemax < (MPI_wtime()-time_total)/3600.d0) then
      if (root) write(*,*) "Out of walltime!"
      continue_timestepping=.false.
    endif

    !-----------------------------------------------
    ! Runtime remote control (every 10 time steps)
    !-----------------------------------------------
    if (modulo(it,10)==0) then
      ! fetch command from file
      call runtime_control_command( command )
      ! execute it
      select case ( command )
      case ("reload_params")
        if (root) write (*,*) "runtime control: Reloading PARAMS file.."
        ! read all parameters from the params.ini file
        call get_params(params_file,Insect,.true.)
        ! overwrite control file
        if (root) call initialize_runtime_control_file()
      case ("save_stop")
        if (root) write (*,*) "runtime control: Safely stopping..."
        continue_timestepping = .false. ! this will stop the time loop
        ! overwrite control file
        if (root) call initialize_runtime_control_file()
      end select
    endif

    if (root) call save_time_stepping_info(time, dt1, it, t4)

    ! timing module: sum the time for all time steps
    call toc("MAIN (complete time stepping loop)", MPI_wtime()-t4)
  enddo


  !-----------------------------------------------------------------------------
  ! save final backup so we can resume where we left
  !-----------------------------------------------------------------------------
  if(idobackup==1) then
    if(root) write (*,*) "final backup..."
    call dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,&
    work(:,:,:,1),scalars,scalars_rhs,Insect,beams,wings)
  endif

  if(root) write(*,'("Done time stepping; did nt=",i7," steps")') it-it_start
end subroutine time_step


!-------------------------------------------------------------------------------
! dump information about performance, time, time step and iteration
! number to disk. this helps estimating how many time steps will be
! required when passing to higher resolution. This routine is called after every
! time step
!-------------------------------------------------------------------------------
subroutine save_time_stepping_info(time, dt, it, wtime_start)
  use vars
  implicit none

  integer,intent(in) :: it ! current time step number
  real(kind=pr),intent(in) :: time, dt
  real(kind=pr),intent(in) :: wtime_start

  open  (14,file='timesteps_info.t',status='unknown',position='append')
  write (14,'(2(g15.8,1x),i9,1x,g15.8,1x,i6)') time, MPI_wtime()-wtime_start, it, dt, mpisize
  close (14)

end subroutine save_time_stepping_info


!-------------------------------------------------------------------------------
! Output how much time remains in the simulation.
!-------------------------------------------------------------------------------
subroutine are_we_there_yet(time, wtime_tstart, dt1)
  use vars
  implicit none

  real(kind=pr),intent(inout) :: time, wtime_tstart, dt1
  real(kind=pr):: time_left, t2

  if (mpirank==0) then
    ! compute the time elapsed since the time stepping started
    t2 = MPI_wtime() - wtime_tstart

    ! compute (estimated) remaining time
    time_left = (tmax-time) * ( t2 / (time-tstart) )

    ! print information
    write(*,'("time left: ",i2,"d ",i2,"h ",i2,"m ",i2,"s wtime=",f4.1,"h dt=",es10.2,"s t=",g10.2)') &
    floor(time_left/(24.d0*3600.d0))   ,&
    floor(mod(time_left,24.d0*3600.d0)/3600.d0),&
    floor(mod(time_left,3600.d0)/60.d0),&
    floor(mod(mod(time_left,3600.d0),60.d0)),&
    (MPI_wtime()-time_total)/3600.d0,&
    dt1,time
  endif
end subroutine are_we_there_yet
