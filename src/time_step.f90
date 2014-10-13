subroutine time_step(time,u,nlk,work,mask,mask_color,us,Insect,beams,params_file)
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none
  
  
  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  character(len=strlen),intent(in)  :: params_file ! for runtime control  
  
  logical :: continue_timestepping
  integer :: inter
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  real(kind=pr) :: t1,t2,t3,t4
  character(len=strlen)  :: command ! for runtime control
 
  ! first backup is in "truntime" hours
  truntimenext = truntime 
  continue_timestepping = .true.
  time%it_start = time%it

  ! After init, output integral quantities
  if (index(inicond,'backup::')==0 .and. dry_run_without_fluid/="yes") then
    if (root) write(*,*) "Initial output of integral quantities...."
    call write_integrals(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  endif


  if (root) write(*,*) "Start time-stepping...."
  
  ! Loop over time steps
  t1=MPI_wtime()
  do while ((time%time<tmax).and.(time%it<=nt).and.(continue_timestepping))
     t4=MPI_wtime()
     !-------------------------------------------------
     ! advance fluid/B-field in time
     !-------------------------------------------------
     if(dry_run_without_fluid/="yes") then
       ! make a time step of the fluid
       call fluidtimestep(time,u,nlk,work,mask,mask_color,us,Insect,beams)
     endif
     
     !-----------------------------------------------
     ! time step done: advance iteration + time
     !-----------------------------------------------
     ! change registers for AB2:
     inter=time%n1 ; time%n1=time%n0 ; time%n0=inter
     ! Advance in time so that uk contains the evolved field at time 'time+dt1'
     time%time = time%time + time%dt_new
     time%dt_old=time%dt_new
     time%it = time%it + 1
     
     !-------------------------------------------------
     ! Output of INTEGRALS after every tintegral time 
     ! units or itdrag time steps
     !-------------------------------------------------
     if ((modulo(time%time,tintegral) <= time%dt_new .or. &
          modulo(time%it,itdrag)==0) .and. dry_run_without_fluid/="yes" ) then
       call write_integrals(time,u,nlk,work,mask,mask_color,us,Insect,beams)
     endif
    
     !-------------------------------------------------
     ! Save solid model data, if using it
     !-------------------------------------------------
     if (use_solid_model=="yes" .and. mpirank==0 .and. modulo(time%it,itbeam)==0) then
       call SaveBeamData( time%time, beams )
     endif
    
     !-------------------------------------------------
     ! Output FIELDS+BACKUPING (after tsave)
     !-------------------------------------------------
     if (((modulo(time%time,tsave)<time%dt_new).and.(time%time>=tsave_first)).or.(time%time==tmax)) then
        call are_we_there_yet(time%it,time%it_start,time%time,t2,t1,time%dt_new)
        call save_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)   
     endif

     ! Backup if that's specified in the PARAMS.ini file. We try to do 
     ! backups every "truntime" hours (precise to one time step)
     if(idobackup==1 .and. truntimenext<(MPI_wtime()-time_total)/3600.d0) then
         call dump_runtime_backup(time,nbackup,u,Insect,beams)
         truntimenext = truntimenext+truntime
     endif
     
     !-----------------------------------------------
     ! Output how much time remains
     !-----------------------------------------------
     if ((modulo(time%it,300)==0).or.(time%it==20)) then
        call are_we_there_yet(time%it,time%it_start,time%time,t2,t1,time%dt_new)
     endif
     
     !-------------------------------------------------
     ! escape from loop if walltime is about to be exceeded
     !-------------------------------------------------     
     if(wtimemax < (MPI_wtime()-time_total)/3600.d0) then
        if (root) write(*,*) "Out of walltime!"
        continue_timestepping=.false.   
     endif
     
     !-----------------------------------------------
     ! Runtime remote control (every 10 time steps)
     !-----------------------------------------------
     if (modulo(time%it,10)==0) then
        ! fetch command from file
        call runtime_control_command( command )
        ! execute it
        select case ( command )
        case ("reload_params")
          if (root) write (*,*) "runtime control: Reloading PARAMS file.."
          ! read all parameters from the params.ini file
          call get_params(params_file,Insect)           
          ! overwrite control file
          if (root) call initialize_runtime_control_file()
        case ("save_stop")
          if (root) write (*,*) "runtime control: Safely stopping..."
          continue_timestepping = .false. ! this will stop the time loop
          ! overwrite control file
          if (root) call initialize_runtime_control_file()
        end select
     endif 
     
     if(root) call save_time_stepping_info (time%it,time%it_start,time%time,t2,t1,time%dt_new,t4)
  enddo

  !-----------------------------------------------------------------------------
  ! save final backup so we can resume where we left 
  !-----------------------------------------------------------------------------
  if(idobackup==1) then
    if(root) write (*,*) "final backup..."
    call dump_runtime_backup(time,nbackup,u,Insect,beams)
  endif

  if(root) write(*,'("Done time stepping; did nt=",i5," steps")') time%it-time%it_start
end subroutine time_step


!-------------------------------------------------------------------------------
! dump information about performance, time, time step and iteration
! number to disk. this helps estimating how many time steps will be 
! required when passing to higher resolution. This routine is called after every
! time step
!-------------------------------------------------------------------------------
subroutine save_time_stepping_info(it,it_start,time,t2,t1,dt1,t3)
  use vars
  implicit none

  real(kind=pr),intent(inout) :: time,t2,t1,dt1,t3  
  integer,intent(inout) :: it,it_start
  
  ! t2 is time [sec] per time step, averaged
  t2 = (MPI_wtime() - t1) / dble(it-it_start)
  ! t3 is the time per time step, not averaged (on input, t3 is start point of 
  ! time step
  t3 = MPI_wtime() - t3
  
  open  (14,file='timestep.t',status='unknown',position='append')
  write (14,'(i15,1x,es15.8,1x,3(es15.8,1x))') it,time,dt1,t2,t3
  close (14)

end subroutine save_time_stepping_info


!-------------------------------------------------------------------------------
! Output how much time remains in the simulation.
!-------------------------------------------------------------------------------
subroutine are_we_there_yet(it,it_start,time,t2,t1,dt1)
  use vars
  implicit none

  real(kind=pr),intent(inout) :: time,t2,t1,dt1
  integer,intent(inout) :: it,it_start
  real(kind=pr):: time_left

  ! This is done every 300 time steps, but it may happen that this is
  ! too seldom or too often. in future versions, maybe we try doing it
  ! once in an hour or so. We also output a first estimate after 20
  ! time steps.
  if(root) then  
     t2= MPI_wtime() - t1
     time_left=(((tmax-time)/dt1)*(t2/dble(it-it_start)))
     write(*,'("time left: ",i3,"d ",i2,"h ",i2,"m ",i2,"s wtime=",f4.1,"h dt=",es10.2,"s t=",es10.2)') &
          floor(time_left/(24.d0*3600.d0))   ,&
          floor(mod(time_left,24.*3600.d0)/3600.d0),&
          floor(mod(time_left,3600.d0)/60.d0),&
          floor(mod(mod(time_left,3600.d0),60.d0)),&
          (MPI_wtime()-time_total)/3600.d0,&
          dt1,time
  endif
end subroutine are_we_there_yet
