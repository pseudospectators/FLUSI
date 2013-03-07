subroutine cal_nlk (dt1, nlk, uk, u , vort, work )
  ! -------------------------------------------------------------------------------------------
  ! 		VERSION 5 / mars / 2013
  ! - no allocating, all work arrays passed as arguments
  ! - this version 7 real work arrays; I' afraid this is the minimum.
  ! - cofdx..cofdz and poisson are replaced by more elegant loops to minimize
  !   both comput. time and memory consumption.
  !
  ! ----------------------------
  ! This routine computes the RHS of Navier-stokes eqn in fourier space.
  ! -------------------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (in) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (out):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), intent (out) :: dt1
  real (kind=pr) :: u_max_w, divu, divumax, t1=0.d0,t2=0.d0,t3=0.d0,t4=0.d0,t5=0.d0,t6=0.d0
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  real (kind=pr), dimension (3) :: u_max, u_loc
  integer :: mpicode,ix,iy,iz
  complex (kind=pr) :: qk

  

  !-----------------------------------------------
  !--Calculate ux and uy in physical space
  !-----------------------------------------------
  call cofitxyz (uk(:,:,:,1), u(:,:,:,1))
  call cofitxyz (uk(:,:,:,2), u(:,:,:,2))
  call cofitxyz (uk(:,:,:,3), u(:,:,:,3))

  !-----------------------------------------------
  !-- compute vorticity
  !-----------------------------------------------
  do iy = ca(3), cb(3)  		! ky : 0..ny/2-1 ,then, -ny/2..-1     
    ky = scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
    do ix = ca(2), cb(2)		! kx : 0..nx/2
	kx = scalex*dble(ix)                
	do iz = ca(1),cb(1)		! kz : 0..nz/2-1 ,then, -nz/2..-1           
	  kz = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
	  nlk(iz,ix,iy,1) = dcmplx(0d0,1d0)*( ky*uk(iz,ix,iy,3) - kz*uk(iz,ix,iy,2) )
	  nlk(iz,ix,iy,2) = dcmplx(0d0,1d0)*( kz*uk(iz,ix,iy,1) - kx*uk(iz,ix,iy,3) )
	  nlk(iz,ix,iy,3) = dcmplx(0d0,1d0)*( kx*uk(iz,ix,iy,2) - ky*uk(iz,ix,iy,1) )
	enddo
    enddo
  enddo
  call cofitxyz (nlk(:,:,:,1), vort(:,:,:,1)) ! Transform it to physical space
  call cofitxyz (nlk(:,:,:,2), vort(:,:,:,2)) ! Transform it to physical space
  call cofitxyz (nlk(:,:,:,3), vort(:,:,:,3)) ! Transform it to physical space
 
      
  !-------------------------------------------------------------
  !-- Calculate omega x u (cross-product)
  !-- and transform the result into Fourier space 
  !-------------------------------------------------------------
  if ( iobst > 0 ) then
     ! -- x component
     work = u(:,:,:,2) * vort(:,:,:,3) - u(:,:,:,3) * vort(:,:,:,2) - u(:,:,:,1)*mask
     call coftxyz (work, nlk(:,:,:,1))
     ! -- y component
     work = u(:,:,:,3) * vort(:,:,:,1) - u(:,:,:,1) * vort(:,:,:,3) - u(:,:,:,2)*mask
     call coftxyz (work, nlk(:,:,:,2))
     ! -- z component
     work = u(:,:,:,1) * vort(:,:,:,2) - u(:,:,:,2) * vort(:,:,:,1) - u(:,:,:,3)*mask
     call coftxyz (work, nlk(:,:,:,3))
  else
     ! -- x component
     work = u(:,:,:,2) * vort(:,:,:,3) - u(:,:,:,3) * vort(:,:,:,2)
     call coftxyz (work, nlk(:,:,:,1))
     ! -- y component
     work = u(:,:,:,3) * vort(:,:,:,1) - u(:,:,:,1) * vort(:,:,:,3)
     call coftxyz (work, nlk(:,:,:,2))
     ! -- z component
     work = u(:,:,:,1) * vort(:,:,:,2) - u(:,:,:,2) * vort(:,:,:,1)
     call coftxyz (work, nlk(:,:,:,3))  
  endif  
  

  ! add this place, we will one day compute drag/lift forces, if required. they are nothing but the sum of the penalty term

  
  !-------------------------------------------------------------
  !-- add pressure, new version
  !-- p = (i*kx*sxk + i*ky*syk + i*kz*szk) / k**2
  !-- note: we use rotational formulation: p is NOT the physical pressure
  !------------------------------------------------------------- 
  do iy = ca(3), cb(3)  		! ky : 0..ny/2-1 ,then, -ny/2..-1     
    ky = scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
    ky2 = ky*ky
    do ix = ca(2), cb(2)		! kx : 0..nx/2
	kx = scalex*dble(ix)                
	kx2 = kx*kx
	do iz = ca(1),cb(1)		! kz : 0..nz/2-1 ,then, -nz/2..-1           
	  kz      = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
	  kz2     = kz*kz
	  k_abs_2 = kx2+ky2+kz2
	  if (abs(k_abs_2) .ne. 0.0) then  
	  qk    = (kx*nlk(iz,ix,iy,1) + ky*nlk(iz,ix,iy,2) + kz*nlk(iz,ix,iy,3)) / k_abs_2
	  nlk(iz,ix,iy,1) = nlk(iz,ix,iy,1) - kx*qk 	! add gradient of pressure
	  nlk(iz,ix,iy,2) = nlk(iz,ix,iy,2) - ky*qk
	  nlk(iz,ix,iy,3) = nlk(iz,ix,iy,3) - kz*qk
	  endif
	enddo
    enddo
  enddo  
  

  !-------------------------------------------------------------
  !--Calculate maximum velocity + TIME STEP
  !-------------------------------------------------------------
  u_loc(1) = maxval(abs(u(:,:,:,1)))  ! local maximum of x-velocity magnitude
  u_loc(2) = maxval(abs(u(:,:,:,2)))  ! local maximum of y-velocity magnitude
  u_loc(3) = maxval(abs(u(:,:,:,3)))  ! local maximum of z-velocity magnitude
  
  call MPI_REDUCE (u_loc, u_max, 3, mpireal, MPI_MAX, 0, MPI_COMM_WORLD, mpicode)  ! max at 0th process

  !--Adjust time step at 0th process
  if ( mpirank == 0 ) then
     u_max_w = max ( u_max(1)/dx, u_max(2)/dy, u_max(3)/dz )
     dt1 = cfl / u_max_w        
     call truncate(dt1,dt1) ! round time step to one digit: saves time because no need to remocompute cal_vis
     if (iobst > 0 .and. dt1 >= eps) dt1 = 0.99*eps ! time step is smaller than eps (Dmitry, 26 feb 08)  
  endif 

  !-- Broadcast time step to all processes
  call MPI_BCAST (dt1, 1, mpireal, 0, MPI_COMM_WORLD, mpicode)
  
end subroutine cal_nlk

subroutine truncate(a,b)
  use share_vars
  implicit none
  real(kind=pr) :: a,b
  character (len=7) :: str

  write (str,'(es7.1)') a
  read (str,*) b
  ! write (*,'(2(es15.8,1x))') a,b


end subroutine


  !-------------------------------------------------------------
  ! compute divergence at the old time step
  !-------------------------------------------------------------

!   work2k = vxk
!   call cofdx (work2k)
!   work1k = vyk
!   call cofdy (work1k)
!   work1k = work1k + work2k
!   work2k = vzk
!   call cofdz (work2k)
!   work1k = work1k + work2k
!   ! work1k now contains the divergence
!   call cofitxyz (work1k, work)
!   ! work now contains the divergence in physical space
!   divu = maxval(abs(work))
!   call MPI_REDUCE (divu, divumax, 1, mpireal, MPI_MAX, 0, MPI_COMM_WORLD, mpicode)  ! max at 0th process
!   if ( mpirank == 0) then
!    write(*,'("div(u)=",es12.4)') divumax
!   endif