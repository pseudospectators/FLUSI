subroutine time_step(u,uk,nlk,vort,work,explin,params_file,time,dt0,dt1,n0,n1,it)
  use mpi
  use vars
  use fsi_vars
  use p3dfft_wrapper
  implicit none
  
  integer :: inter
  integer,intent(inout) :: n0,n1,it
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  integer :: it_start
  real(kind=pr),intent(inout) :: time,dt0,dt1 
  real(kind=pr) :: t1,t2
  character (len=80)  :: command ! for runtime control
  character (len=80),intent(in)  :: params_file ! for runtime control  
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  logical :: continue_timestepping
  
  continue_timestepping = .true.
  it_start=it
  
  ! save initial conditions 
  call save_fields_new(time,uk,u,vort,nlk(:,:,:,:,n0),work)    
  
  ! initialize runtime control file
  if (mpirank == 0) call initialize_runtime_control_file()
   
  ! After init, output integral quantities. (note we can overwrite only 
  ! nlk(:,:,:,:,n0) when retaking a backup)
  if (mpirank == 0) write(*,*) "Initial output of integral quantities...."
  call write_integrals(time,uk,u,vort,nlk(:,:,:,:,n0),work)


  if (mpirank == 0) write(*,*) "Start time-stepping...."
  
  ! Loop over time steps
  t1=MPI_wtime()
  do while ((time<=tmax) .and. (it<=nt) .and. (continue_timestepping) )
     dt0=dt1
     !-------------------------------------------------
     ! If the mask is time-dependend,we create it here
     !-------------------------------------------------
     if(iMoving == 1 .and. iPenalization == 1) call create_mask(time)
     
     !-------------------------------------------------
     ! advance fluid/B-field in time
     !-------------------------------------------------
     if(dry_run_without_fluid/="yes") then
       call fluidtimestep(time,dt0,dt1,n0,n1,u,uk,nlk,vort,work,explin,it)
     endif

     !-------------------------------------------------
     ! Compute hydrodynamic forces at time level n (FSI only)
     ! NOTE:    is done every itdrag time steps. this condition will be changed
     !          in future versions if free-flight (ie solving eq of motion) is
     !          required.
     !-------------------------------------------------
     if ( method=="fsi" ) then
        ! compute unst corrections in every time step
        if (unst_corrections ==1) then
          call cal_unst_corrections ( time, dt0 )
        endif
        ! compute drag only if required
        if ((modulo(it,itdrag)==0).or.(SolidDyn%idynamics/=0)) then
        if (compute_forces==1) then
          call cal_drag ( time, u ) ! note u is OLD time level
        endif
        endif
        ! note dt0 is OLD time step t(n)-t(n-1)
        ! advance in time ODEs that describe rigid solids
        if (SolidDyn%idynamics==1) call rigid_solid_time_step(time,dt0,dt1,it)
     endif     
     
     !-----------------------------------------------
     ! time step done: advance iteration + time
     !-----------------------------------------------
     inter=n1 ; n1=n0 ; n0=inter
     ! Advance in time so that uk contains the evolved field at time 'time+dt1'
     time=time + dt1
     it=it + 1
     
     !-------------------------------------------------
     ! Output of INTEGRALS after every tintegral time 
     ! units or itdrag time steps
     !-------------------------------------------------
     if ((modulo(time,tintegral) <= dt1).or.(modulo(it,itdrag) == 0)) then
       call write_integrals(time,uk,u,vort,nlk(:,:,:,:,n0),work)
     endif
    
     !-------------------------------------------------
     ! Output FIELDS (after tsave)
     !-------------------------------------------------
     if (((modulo(time,tsave)<=dt1).and.(it>2)).or.(time==tmax)) then
        call are_we_there_yet(it,it_start,time,t2,t1,dt1)
        ! Note: we can safely delete nlk(:,:,:,1:nd,n0). for RK2 it
        ! never matters,and for AB2 this is the one to be overwritten
        ! in the next step.  This frees 3 complex arrays, which are
        ! then used in Dump_Runtime_Backup.
        call save_fields_new(time,uk,u,vort,nlk(:,:,:,:,n0),work)       
        
        ! Backup if that's specified in the PARAMS.ini file
        if(iDoBackup == 1) then
           call dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
        endif
     endif
     
     !-----------------------------------------------
     ! Output how much time remains
     !-----------------------------------------------
     if ((modulo(it,300)==0).or.(it==20)) then
       call are_we_there_yet(it,it_start,time,t2,t1,dt1)
     endif
     
     !-----------------------------------------------
     ! Runtime remote control (every 10 time steps)
     !-----------------------------------------------
     if ( modulo(it,10) == 0 ) then
        ! fetch command from file
        call runtime_control_command( command )
        ! execute it
        select case ( command )
        case ("reload_params")
          if (mpirank==0) write (*,*) "runtime control: Reloading PARAMS file.."
          ! read all parameters from the params.ini file
          call get_params(params_file)           
          ! overwrite control file
          if (mpirank == 0) call initialize_runtime_control_file()
        case ("save_stop")
          if (mpirank==0) write (*,*) "runtime control: Safely stopping..."
          call dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
          continue_timestepping = .false. ! this will stop the time loop
          ! overwrite control file
          if (mpirank == 0) call initialize_runtime_control_file()
        end select
     endif 
     
     if(mpirank==0) call save_time_stepping_info (it,it_start,time,t2,t1,dt1)
  end do

  if(mpirank==0) then
     write(*,'("Finished time stepping; did it=",i5," time steps")') it
  endif
end subroutine time_step


!-------------------------------------------------------------------------------
! dump information about performance, time, time step and iteration
! number to disk. this helps estimating how many time steps will be 
! required when passing to higher resolution.
!-------------------------------------------------------------------------------
subroutine save_time_stepping_info(it,it_start,time,t2,t1,dt1)
  use vars
  implicit none

  real(kind=pr),intent(inout) :: time,t2,t1,dt1
  integer,intent(inout) :: it,it_start
  
  if (mpirank == 0) then
  ! t2 is time [sec] per time step
  t2 = (MPI_wtime() - t1) / dble(it-it_start)
  
  open  (14,file='timestep.t',status='unknown',position='append')
  write (14,'(e12.5,A,i7.7,A,es12.5,A,es12.5)') time,tab,it,tab,dt1,tab,t2
  close (14)
  endif
 
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

  ! this is done every 300 time steps, but it may happen that this is too seldom
  ! or too oftern. in future versions, maybe we try doing it once in an hour or
  ! so. we also output a first estimate after 20 time steps
  if(mpirank == 0) then  
     t2= MPI_wtime() - t1
     time_left=(((tmax-time)/dt1)*(t2/dble(it-it_start)))
     write(*,'("time left: ",i3,"d ",i2,"h ",i2,"m ",i2,"s wtime=",f4.1," dt=",es10.2,"s t=",es10.2)') &
          floor(time_left/(24.d0*3600.d0))   ,&
          floor(mod(time_left,24.*3600.d0)/3600.d0),&
          floor(mod(time_left,3600.d0)/60.d0),&
          floor(mod(mod(time_left,3600.d0),60.d0)),&
          (MPI_wtime()-time_total)/3600.d0,&
          dt1,time
  endif
end subroutine are_we_there_yet
