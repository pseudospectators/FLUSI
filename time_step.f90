subroutine time_step(u,uk,nlk,vort,work,workvis)
  use mpi_header
  use fsi_vars
  implicit none
  
  integer :: inter,it
  integer :: n0=0,n1=1
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  integer :: it_start
  real(kind=pr) :: time,dt0,dt1,t1,t2
  
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent (inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent (inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent (inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent (inout) :: workvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  if (mpirank == 0) write(*,'(A)') 'Info: Starting time iterations.'

  time=0.0

  dt0=1.0d0 ! useful to trigger cal_vis
  dt1=2.d0  ! just add a comment to test branching...
     
  ! Initialize vorticity or read values from a backup file
  ! FIXME: move to FLUSI.f90 or mhd.f90?
  call init_fields(n1,time,it,dt0,dt1,uk,nlk,vort,workvis)
  n0=1 - n1
  it_start=it ! FIXME: when moving the above, it would be nice to move this as well. 

  ! Loop over time steps
  t1=MPI_wtime()
  do while((time<=tmax) .and. (it<=nt))
     dt0=dt1

     ! If the mask is time-dependend,we create it here
     if(iMoving == 1) call Create_Mask(time)

     ! Do a fluid time step
     call FluidTimeStep(time,dt0,dt1,n0,n1,u,uk,nlk,vort,work,workvis,it)

     ! Switch time levels
     inter=n1 ; n1=n0 ; n0=inter
     ! Advance in time so that uk contains the evolved field at time
     ! 'time+dt1'
     time=time + dt1
     it=it + 1

     ! Output of integrals after itdrag
     if(mpirank == 0 .and. modulo(it,itdrag) == 0) call write_integrals(time)
     
     ! Output how much time remains
     if(mpirank == 0) call are_we_there_yet(it,it_start,time,t2,t1,dt1)

     ! Output(after tsave)
     if(modulo(time - tstart,tsave) <= dt1) then
        ! Note: we can safely delete nlk(:,:,:,1:nf,n0). for RK2 it
        ! never matters,and for AB2 this is the one to be overwritten
        ! in the next step.  This frees 3 complex arrays.
        ! FIXME: why is the above comment important?
        ! FIX: occasionally, I put comments that took me a while to understand
        ! in the code. 
        call save_fields_new(time,uk,u,vort,nlk(:,:,:,:,n0),work)
        
        ! Backup if that's specified in the PARAMS.ini file
        if(iDoBackup == 1) then
           call Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
        endif
     endif
  end do

  if(mpirank==0) then
     write(*,'("Finished time stepping; did it=",i5," time steps")') it
  endif
end subroutine time_step


! Output how much time remains in the simulation.
subroutine are_we_there_yet(it,it_start,time,t2,t1,dt1)
  use vars
  implicit none

  real(kind=pr),intent(inout) :: time,t2,t1,dt1
  integer,intent(inout) :: it,it_start
  real(kind=pr):: time_left

  ! this is done every 300 time steps, but it may happen that this is too seldom
  ! or too oftern. in future versions, maybe we try doing it once in an hour or
  ! so.
  if(modulo(it,300) == 0) then     
     t2= MPI_wtime() - t1
     time_left=(((tmax-time)/dt1)*(t2/dble(it-it_start)))
     write(*,'("time left: ",i3,"d ",i2,"h ",i2,"m ",i2,"s dt=",es7.1)') &
          floor(time_left/(24.d0*3600.d0))   ,&
          floor(mod(time_left,24.*3600.d0)/3600.d0),&
          floor(mod(time_left,3600.d0)/60.d0),&
          floor(mod(mod(time_left,3600.d0),60.d0)),dt1
  endif
end subroutine are_we_there_yet


! Wrapper for writing integral quantities to file
subroutine write_integrals(time)
  use vars
  implicit none

  real(kind=pr), intent(in) :: time

  select case(method(1:3))
  case("fsi") 
     call write_integrals_fsi(time)
  case("mhd") 
     call write_integrals_mhd(time)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in write_integrals"
     call abort
  end select
end subroutine write_integrals


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time)
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  
  open(14,file='drag_data',status='unknown',position='append')
  write(14,'(7(es12.4,1x))')  time,GlobIntegrals%E_kin,&
       GlobIntegrals%Dissip, GlobIntegrals%Force(1),&
       GlobIntegrals%Force(2),GlobIntegrals%Force(3),&
       GlobIntegrals%Volume
  close(14)
end subroutine write_integrals_fsi


! mhd version of writing integral quantities to disk
subroutine write_integrals_mhd(time)
  use mhd_vars
  implicit none

  real(kind=pr), intent(in) :: time
 
  !FIXME: do things here?
end subroutine write_integrals_mhd
