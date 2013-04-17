subroutine time_step 
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  integer :: inter, it, inicond1,ix,iy,iz,vis_tmp
  integer :: n0 = 0, n1 = 1
  integer :: nbackup = 0  ! 0 - backup to file runtime_backup0, 1 - to runtime_backup1, 2 - no backup 
  integer :: mpicode
  real (kind=pr) :: time, v_rms, v_rms_loc, dt0, dt1,t1,t2
  real (kind=pr), 	dimension (:,:,:), allocatable :: workvis  
  real (kind=pr), 	dimension (:,:,:,:), allocatable :: u, vort
  complex (kind=pr), 	dimension (:,:,:,:), allocatable :: uk
  complex (kind=pr), 	dimension (:,:,:,:,:), allocatable :: nlk  
  real (kind=pr), 	dimension (:,:,:), allocatable :: work
  type(Integrals) :: GlobIntegrals
  GlobIntegrals%E_kin  = 0.d0
  GlobIntegrals%Dissip = 0.d0
  GlobIntegrals%Force  = 0.d0
  
  dt0 = 1.0d0 ! useful to trigger cal_vis

  !---------------------------------------------------------------
  !-- Allocate memory
  !---------------------------------------------------------------
  allocate ( workvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)) )
  allocate ( uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3) )		! velocity in fourier space
  allocate ( nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1) )
  allocate ( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) ) 		! velocity in phy space
  allocate ( vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) ) 		! vorticity in phy space
  allocate ( work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  
  !---------------------------------------------------------------
  ! Create obstacle mask
  !---------------------------------------------------------------
  if (iPenalization>0) then 
    allocate ( mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
    if (iMoving==1) then
    allocate ( us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )
    endif
    call obst_mask()
  endif
  
  !---------------------------------------------------------------
  ! Initialize vorticity or read values from a backup file
  !---------------------------------------------------------------
  call init_fields (n1, time, it,  dt0, dt1, uk, nlk, vort, workvis )
  n0 = 1 - n1
  
  !---------------------------------------------------------
  !     LOOP OVER TIME STEPS
  !---------------------------------------------------------
  t2=0.d0
  
  
  t1=MPI_wtime()
  do while ((time<=tmax).and.(it<=nt))
     dt0 = dt1
   
     !-----------------------------
     !     do a fluid time step
     !-----------------------------
     if (imask == 555) then
      call create_mask ( time )  
     elseif (imask == 666) then
      call Draw_Jerry( time )
     endif
     
     call FluidTimeStep ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, it, GlobIntegrals )
     

     !--Switch time levels
     inter = n1 ; n1 = n0 ; n0 = inter
     !--Advance in time: at this point, uk contains the velocity field at time 'time' 
     time = time + dt1
     it = it + 1

     !---------------------------------------------------------------------------------
     !--				Output (after tdrag)
     !---------------------------------------------------------------------------------
     if ( modulo (it,itdrag) == 0 ) then
       if (mpirank==0) then
	open  (14, file = 'drag_data', status = 'unknown', access = 'append')
	write (14, '(7(es12.4,1x))')  time, GlobIntegrals%E_kin, &
	GlobIntegrals%Dissip,  GlobIntegrals%Force(1),&
	GlobIntegrals%Force(2),GlobIntegrals%Force(3),&
	GlobIntegrals%Volume
	close (14)
	endif
     endif
       
     if ( modulo (it,300) == 0 ) then     
	if (mpirank==0) then
	t2= MPI_wtime() - t1
	write(*,'("time left: ",f7.2,"min dt=",es7.1)') (((tmax-time)/dt1)*(t2/dble(it)))/60.d0 , dt1
	endif
     endif
     
     !---------------------------------------------------------------------------------
     !--				Output (after tsave)
     !---------------------------------------------------------------------------------
     if( modulo (time - tstart, tsave) <= dt1 ) then
	! note: we can safely delete nlk(:,:,:,1:3,n0). for RK2 it never matters, 
	! and for AB2 this is the one to be overwritten in the next step.
	! this frees 3 complex arrays
	call save_fields_new ( time, dt1, uk, u, vort, nlk(:,:,:,:,n0), work)
	! Backup if that's specified in the PARAMS.ini file
	if (iDoBackup==1) call Dump_Runtime_Backup ( time, dt0, dt1, n1,it,  nbackup, uk, nlk, workvis)
     endif

  end do




  if (mpirank==0) then
  v_rms = dsqrt ( sum( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)  )
  write(*,'("u_rms=",es15.8," max=",es15.8)') v_rms, maxval( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
  write(*,'("did it=",i5," time steps")') it
  endif


  !-- Deallocate memory
  deallocate ( workvis )
  deallocate ( vort, work )
  deallocate ( u, uk, nlk)
  
  if ( iPenalization == 1 ) then
     deallocate ( mask )
     if (iMoving == 1) deallocate (us)
  endif
!---------------------------------------------------------
!     END LOOP OVER TIME STEPS
!---------------------------------------------------------
end subroutine time_step

