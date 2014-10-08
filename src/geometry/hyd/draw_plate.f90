!
! Mask and velocity field of a plate, ellipse etc.
! Comparison with flapping foil, Wang, PRL 2000
!
subroutine Draw_Plate (time, mask, mask_color, us)
  use vars
  implicit none

  
  real(kind=pr),intent(in) :: time
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  integer :: ix, iy, iz, mpicode
  real (kind=pr) :: x2, y2, vx2, vy2, vx2t, vy2t, anglez2, omz2, omz2t
  real (kind=pr) :: x, y, xref, yref, xlev, ylev, tmp, N, rref, rmax, hsmth, Am, alpham
  real (kind=pr) :: Af, Sxf, Syf, Jf, forcex, forcey, torquez

  N = 3.0d0 ! smoothing coefficient
  hsmth = N*dx ! smoothing layer thickness
  rmax = 0.5d0

  ! Flapping parameters
  Am = 1.25d0
  alpham = 0.25d0*pi 

  ! Update kinematics
  x2 = x0 + Am * dcos(time/Am)
  y2 = y0
  anglez2 = 0.5d0*pi + alpham * dsin(time/Am)
  vx2 = - dsin(time/Am) 
  vy2 = 0.0d0
  omz2 = alpham/Am * dcos(time/Am)
  vx2t = - 1.0/Am * cos(time/Am)     
  vy2t = 0.0d0
  omz2t = - alpham/Am**2 * sin(time/Am)

  ! For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)

     ! Create mask: smooth indicator function of an ellipse
     x = dble(ix)*dx - x2
     y = dble(iy)*dy - y2
     xref = x*dcos(anglez2) + y*dsin(anglez2)
     yref = y*dcos(anglez2) - x*dsin(anglez2)
     rref = dsqrt( xref**2 + 4.0d0**2 * yref**2 ) ! Radius in cylindrical coordinates
     call SmoothStep (tmp, rref, rmax-0.0d0*hsmth, hsmth)
     mask(ix,iy,iz) = tmp

     ! assign color "1" where >0 indicates something "useful"
     if (tmp > 1.0e-12) mask_color(ix,iy,iz) = 1
     
     ! Velocity
     if ( rref < rmax+6.0d0*hsmth ) then
      us(ix,iy,iz,1) = -omz2*y + vx2
      us(ix,iy,iz,2) = +omz2*x + vy2
     endif
    enddo
   enddo
  enddo

  ! Calculate unsteady correction of forces
  Af = 0.0d0
  Sxf = 0.0d0
  Syf = 0.0d0
  Jf = 0.0d0 
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    ylev = real(iy)*dy - y2
     do ix =  ra(1), rb(1)
       xlev = real(ix)*dx - x2

       ! Volume, static and inertial moments
       Af = Af + mask(ix,iy,iz)
       Sxf = Sxf + xlev * mask(ix,iy,iz)
       Syf = Syf + ylev * mask(ix,iy,iz)
       Jf = Jf + ( xlev**2 + ylev**2 ) * mask(ix,iy,iz)
     enddo
   enddo
  enddo
  Af = Af * dx*dy*dz      
  Sxf = Sxf * dx*dy*dz 
  Syf = Syf * dx*dy*dz
  Jf = Jf * dx*dy*dz

  ! Unsteady correction
  forcex = vx2t * Af - omz2t * Syf - omz2**2 * Sxf
  forcey = vy2t * Af + omz2t * Sxf - omz2**2 * Syf
  torquez = omz2t * Jf

  ! Reduce 
  GlobalIntegrals%Force(:) = 0.0d0
  GlobalIntegrals%Torque(:) = 0.0d0
  call MPI_REDUCE (forcex,GlobalIntegrals%Force(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)  
  call MPI_REDUCE (forcey,GlobalIntegrals%Force(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode) 
  call MPI_REDUCE (torquez,GlobalIntegrals%Torque(3),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)                    

end subroutine Draw_Plate
