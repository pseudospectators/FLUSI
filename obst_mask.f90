subroutine create_mask (time)
!---------------------------------------------------------------
!     Sets up obstacle mask (with optional boundary 
!     smoothing), size is the object's diameter or length
!---------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real(kind=pr), intent(in) :: time
  integer :: irad, iradmax, ix, iy, iz, ix0, iy0, iz0, kx, ky, kz
  integer :: ixmax, ixmin, iymax, iymin, izmax, izmin
  real (kind=pr) :: x, y, z, xs,ys,zs, alpha,x00,y00,z00, L,B,H, y_tmp, z_tmp

  alpha = 30.d0*pi/180.d0*sin(time*2.d0*pi)
  
  x00=0.5d0*xl
  y00=0.25d0*yl
  z00=0.5d0*zl
  
  L= 1.d0!0.5d0*yl
  H= 3.d0*dz
  B= 0.5d0*xl
  
  mask = 0.d0
  
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
	x = xl*dble(ix)/dble(nx) - x00
	y = yl*dble(iy)/dble(ny) - y00
	z = zl*dble(iz)/dble(nz) - z00
	
	! transformed
	xs = x
	ys = dcos(alpha)*y - dsin(alpha)*z
	zs = dsin(alpha)*y + dcos(alpha)*z

	if (dabs(zs)<H-dz) then
	  z_tmp = 1.d0
	elseif ( (abs(zs)>= H-dz).and.(abs(zs)<H+dz) ) then
	  z_tmp = 1.d0 - (abs(zs) - (H-dz) ) / (2.d0*dz)
	else
	  z_tmp = 0.d0	
	endif
	
	if ( (ys>= dz).and.(ys<L-dz) ) then
	  y_tmp = 1.d0
	elseif ( (ys>=L-dz).and.(ys<=L+dz) ) then
	  y_tmp = 1.d0 - (ys - (L-dz) ) / (2.d0*dz)
	elseif ((ys>=-dz ).and.(ys<=dz)) then 
	  y_tmp = (ys+dz)/(2.d0*dz)
	else
	  y_tmp = 0.d0
	endif
	
	if ( (xs>=-0.5d0*B).and.(xs<=0.5d0*B) ) then
	mask(ix,iy,iz) = z_tmp*y_tmp/eps
	endif
	
      enddo
    enddo
  enddo 

  
  
end subroutine


subroutine obst_mask 
!---------------------------------------------------------------
!     Sets up obstacle mask (with optional boundary 
!     smoothing), size is the object's diameter or length
!---------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  integer :: irad, iradmax, ix, iy, iz, ix0, iy0, iz0, kx, ky, kz
  integer :: ixmax, ixmin, iymax, iymin, izmax, izmin
  real (kind=pr) :: x, y, z


  mask = 0.d0

!---------------------------------------------------------------
! Cubic obstacle
!---------------------------------------------------------------
  if (imask == 1) then

     ixmin = int( ( x0 - size/2.0 ) * real(nx)/xl )
     ixmax = int( ( x0 + size/2.0 ) * real(nx)/xl )
     iymin = int( ( y0 - size/2.0 ) * real(ny)/yl )
     iymax = int( ( y0 + size/2.0 ) * real(ny)/yl )
     izmin = int( ( z0 - size/2.0 ) * real(nz)/zl )
     izmax = int( ( z0 + size/2.0 ) * real(nz)/zl )

     do iz = ra(3), rb(3)
       do iy = ra(2), rb(2)
         do ix = ra(1), rb(1)
           if ( (ix>=ixmin).and.(ix<=ixmax) .and. (iy>=iymin).and.(iy<=iymax) .and. (iz>=izmin).and.(iz<=izmax) ) then
             mask (ix, iy, iz) = 1.0
           endif
         enddo
       enddo
     enddo

  endif

!---------------------------------------------------------------
! Spherical obstacle
!---------------------------------------------------------------
  if (imask == 2) then

     ix0 = int( x0 * real(nx)/xl )
     iy0 = int( y0 * real(ny)/yl )
     iz0 = int( z0 * real(nz)/zl )
     iradmax = int( size/2.0 * real(nx)/xl )

     do iz = ra(3), rb(3)
       do iy = ra(2), rb(2)
         do ix = ra(1), rb(1)
              irad = int( sqrt( real( (ix-ix0)**2 + (iy-iy0)**2 + (iz-iz0)**2 ) ) )
              if ( irad <= iradmax ) then
                 mask (ix, iy, iz) = 1.0
              endif
         enddo
       enddo
     enddo

  endif


!---------------------------------------------------------------
! Cylinder
!---------------------------------------------------------------
  if (imask == 4) then
     ix0 = int( x0 * real(nx)/xl )
     iy0 = int( y0 * real(ny)/yl )
     iradmax = int( size/2.0 * real(nx)/xl )
     do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
           irad = int( sqrt( real( (ix-ix0)**2 + (iy-iy0)**2 ) ) )
           if ( irad <= iradmax ) then
              mask (ix, iy, :) = 1.0
           endif
        enddo
     enddo
  endif


!---------------------------------------------------------------
! dipole-wall first version
!---------------------------------------------------------------
 
  if (imask == 111) then

      do iz = ra(3), rb(3)
       do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
         x = dble(ix)*dx
         y = dble(iy)*dy
         z = dble(iz)*dz
         if ( x > xl-12.d0*dx  ) then
            mask (ix, iy, iz) = 1.0
         endif
        enddo
       enddo
      enddo
  endif






  mask = mask/eps
  
end subroutine obst_mask
