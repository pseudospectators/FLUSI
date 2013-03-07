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
         x = xl*dble(ix)/dble(nx)
         y = yl*dble(iy)/dble(ny)
         z = zl*dble(iz)/dble(nz)
         if ( x > xl-12.d0*xl/dble(nx)  ) then
            mask (ix, iy, iz) = 1.0
         endif
        enddo
       enddo
      enddo
  endif






  mask = mask/eps
  
end subroutine obst_mask
