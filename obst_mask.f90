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
  real (kind=pr) :: t_star, R, alpha_t, un, alpha_max, t1
  real (kind=pr) :: x, y, z, xs,ys,zs, alpha,L,B,H, y_tmp, z_tmp, tmp,N, eps_inv
  t1 = MPI_wtime()
  eps_inv=1.d0/eps
  
  alpha_max = 30.d0*pi/180.d0
  alpha   =           alpha_max*dsin(time*2.d0*pi)
  alpha_t = (2.d0*pi)*alpha_max*dcos(time*2.d0*pi)
  
  y0 = 1.d0
  z0 = 1.d0
  
  L = 1.d0
  H = 0.0625

  ! initialize fields
  mask = 0.d0
  us = 0.d0
   
  N=3.0 ! smoothing coefficient
  
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
	y = dble(iy)*dy - y0
	z = dble(iz)*dz - z0
	
	! transformed
	ys =  dcos(alpha)*y + dsin(alpha)*z
	zs = -dsin(alpha)*y + dcos(alpha)*z
	
	if ( (ys>=0.d0) .and. (ys<=L) )  then
	  
	  call SmoothStep (tmp, abs(zs), H, N*max(dx,dy,dz))
	  mask(:,iy,iz) = tmp*eps_inv
	  
	  if (dabs(zs)<=2.d0*H) then
	    R = dsqrt( y**2 + z**2  )
	    un = R*alpha_t
	    us(:,iy,iz,2) = -dsin(alpha)*un
	    us(:,iy,iz,3) = +dcos(alpha)*un
	  endif
	  
	endif
    enddo
  enddo 

  ! draw the endpoints
  do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
	  y = dble(iy)*dy - y0
	  z = dble(iz)*dz - z0
	  
	  ! origin
	  if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
	    R = dsqrt( y**2 + z**2 )
	    call SmoothStep (tmp, R, H, N*max(dx,dy,dz))
	    if (mask(1,iy,iz)<tmp) then
	      mask(:,iy,iz) = tmp*eps_inv
	      
	      R = dsqrt( y**2 + z**2  )
	      un = R*alpha_t
	      us(:,iy,iz,2) = -dsin(alpha)*un
	      us(:,iy,iz,3) = +dcos(alpha)*un
	    endif
	  endif
	  
	  
	  
	  y = dble(iy)*dy - (y0 + L*dcos(alpha))
	  z = dble(iz)*dz - (z0 + L*dsin(alpha))
	  
	  if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
	    R = dsqrt( y**2 + z**2 )
	    call SmoothStep (tmp, R, H, N*max(dx,dy,dz))
	    if (mask(1,iy,iz)<tmp) then
	      mask(:,iy,iz) = tmp*eps_inv

	      y = dble(iy)*dy - y0
	      z = dble(iz)*dz - z0
	      R = dsqrt( y**2 + z**2  )
	      un = R*alpha_t
	      us(:,iy,iz,2) = -dsin(alpha)*un
	      us(:,iy,iz,3) = +dcos(alpha)*un
	    endif
	  endif
	  
      enddo
  enddo  
    
!     write (*,*) maxval(mask)
  
  time_mask = time_mask + MPI_wtime() - t1
  
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
  real (kind=pr) :: x, y, z, t1
  t1 = MPI_wtime()

  mask = 0.d0

!---------------------------------------------------------------
! Cubic obstacle
!---------------------------------------------------------------
  if (imask == 1) then

     ixmin = int( ( x0 - length/2.0 ) * real(nx)/xl )
     ixmax = int( ( x0 + length/2.0 ) * real(nx)/xl )
     iymin = int( ( y0 - length/2.0 ) * real(ny)/yl )
     iymax = int( ( y0 + length/2.0 ) * real(ny)/yl )
     izmin = int( ( z0 - length/2.0 ) * real(nz)/zl )
     izmax = int( ( z0 + length/2.0 ) * real(nz)/zl )

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
     iradmax = int( length/2.0 * real(nx)/xl )

     do iz = ra(3), rb(3)
       do iy = ra(2), rb(2)
         do ix = ra(1), rb(1)
              irad = int( sqrt( real( (ix-ix0)**2 + (iy-iy0)**2 + (iz-iz0)**2 ) ) )
              if ( irad <= iradmax ) then
                 mask (ix, iy, iz) = 1.d0
              endif
         enddo
       enddo
     enddo

  endif
  
  if (imask == 22) then
  us = 0.d0
  
     ix0 = int( x0 * real(nx)/xl )
     iy0 = int( y0 * real(ny)/yl )
     iz0 = int( z0 * real(nz)/zl )
     iradmax = int( length/2.0 * real(nx)/xl )

     do iz = ra(3), rb(3)
       do iy = ra(2), rb(2)
         do ix = ra(1), rb(1)
              irad = int( sqrt( real( (ix-ix0)**2 + (iy-iy0)**2 + (iz-iz0)**2 ) ) )
              if ( irad <= iradmax ) then
                 mask (ix, iy, iz) = 1.d0
              endif              
              if ( ix <= 8  ) then
		mask (ix, iy, iz) = 1.d0
		us(ix,iy,iz,1)=1.d0
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
     iradmax = int( length/2.0 * real(nx)/xl )
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
  
  
  
  time_mask = time_mask + MPI_wtime() - t1  
end subroutine obst_mask




subroutine SmoothStep (f,x,t,h)
  use share_vars
  implicit none
  !-----------------------------------------------------------------
  !-- This subroutine returns the value f of a smooth step function
  !-- The sharp step function would be 1 if x<=t and 0 if x>t
  !-- h is the semi-size of the smoothing area, so
  !-- f is 1 if x<=t-h
  !-- f is 0 if x>t+h
  !-- f is variable (smooth) in between
  !-----------------------------------------------------------------
  real (kind=pr), intent (out) :: f
  real (kind=pr), intent (in)  :: x,t,h
  real (kind=pr) :: a,b,c,d, delta, GradientERF!, SmoothStep
  !--polynomial coefficients:
!   a =  1.0 / (4.0*(h**3))
!   b = -3.0*t / (4.0*(h**3))
!   c =  3.0*(t+h)*(t-h)/(4.0*(h**3))
!   d =  ((t+h)**2)*(2.0*h-t)/(4.0*(h**3))

!  if (x<=t-h) then
!    f = 1.0
!  elseif (((t-h)<x).and.(x<(t+h))) then
! !     f = a*(x**3) + b*(x**2) + c*x + d
! !     f = 1.0 - (x-t+h)/(2.0*h)
!    f = 0.5*(1.+cos((x-t+h)*pi/(2.0*h)) )
!  else
!    f = 0.0
!  endif

  !-----------------------------------
  ! version 14 - error function as non-oscilatory shape
  !-----------------------------------
  ! h - delta (gradient thickness)
  ! t - thickness (radius)
  GradientERF = abs( ( exp(-(2.0*1.0)**2)  - 1.0 )/sqrt(pi) )
  delta = h*GradientERF
  f = 0.5*( erf( (t-x)/delta ) + erf( (x+t)/delta )  )


end subroutine SmoothStep