subroutine time_step (tdrag, ifield, tsave, tstart, nt)
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  integer, intent (in) :: nt
  integer, dimension (3), intent (in) :: ifield
  integer :: inter, it, inicond1,ix,iy,iz,vis_tmp
  integer :: n0 = 0, n1 = 1
  integer :: nbackup = 0  ! 0 - backup to file runtime_backup0, 1 - to runtime_backup1, 2 - no backup 
  integer :: mpicode
  real (kind=pr), intent (in) :: tdrag, tsave, tstart
  real (kind=pr) :: time, v_rms, v_rms_loc, dt0, dt1,t1,t2
  real (kind=pr), 	dimension (:,:,:), allocatable :: workvis  
  real (kind=pr), 	dimension (:,:,:,:), allocatable :: u, vort
  complex (kind=pr), 	dimension (:,:,:,:), allocatable :: uk
  complex (kind=pr), 	dimension (:,:,:,:,:), allocatable :: nlk  
  real (kind=pr), 	dimension (:,:,:), allocatable :: work
  

  interface                                                                
     function vis (dt)                                                     
       use share_vars                                                      
       real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)) :: vis
       real (kind=pr), intent (in) :: dt                                     
     end function vis                                                        
  end interface

  dt0 = 1.0 ! useful to trigger cal_vis

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
  if (iobst>0) then 
    allocate ( mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
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
     endif
     
     call FluidTimeStep ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, it )
     

     
     !--Switch time levels
     inter = n1 ; n1 = n0 ; n0 = inter
     !--Advance in time: at this point, uk contains the velocity field at time 'time' 
     time = time + dt1
     it = it + 1

     !---------------------------------------------------------------------------------
     !--				Output (after tdrag)
     !---------------------------------------------------------------------------------
!      if ( modulo (it,1000) == 0 ) then
! 
!      endif
     
     !---------------------------------------------------------------------------------
     !--				Output (after tsave)
     !---------------------------------------------------------------------------------
     if( modulo (time - tstart, tsave) <= dt1 ) then
	! note: we can safely delete nlk(:,:,:,1:3,n0). for RK2 it never matters, 
	! and for AB2 this is the one to be overwritten in the next step.
	! this frees 3 complex arrays
	call save_fields_new ( time, dt1, uk, u, vort, nlk(:,:,:,:,n0), work)
! 	call Dump_Runtime_Backup ( time, dt0, dt1, n1,it,  nbackup, uk, nlk, workvis)
	t2= MPI_wtime() - t1
	if (mpirank==0) then
	write(*,'("time left: ",f7.2,"min dt=",es7.1)') (((tmax-time)/dt1)*(t2/dble(it)))/60.d0 , dt1
	endif
     endif

  end do




  v_rms = dsqrt ( sum( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)  )
  write(*,'("u_rms=",es15.8," max=",es15.8)') v_rms, maxval( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )


  !-- Deallocate memory
  deallocate ( workvis )
  deallocate ( vort, work )
  deallocate ( u, uk, nlk)
  
  if ( iobst > 0 ) then
     deallocate ( mask )
  endif
!---------------------------------------------------------
!     END LOOP OVER TIME STEPS
!---------------------------------------------------------
end subroutine time_step

