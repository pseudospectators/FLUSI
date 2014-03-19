!-------------------------------------------------------------------------------
!-- Draws a flexible plate, where the deflection line is solved using the beam 
!-- solver. 
!-------------------------------------------------------------------------------
subroutine Draw_flexible_plate (time)!, beam)
  use mpi
  use fsi_vars
  !-- use global variables for solid solver, as well as routines
  use solid_model 
  implicit none

  real(kind=pr), intent(in) :: time
  type(Solid) :: beam
  real(kind=pr) :: psi, gamma, R, tmp
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate,M1,M2,M3
  !-- for the triangles:
  real(kind=pr) :: a,b,c,alpha,beta,h
  integer :: ix,iy,iz,is
  
  type(solid), dimension(1:nBeams) :: beams
  real (kind=pr) :: time2=0.d0
  

  
  !-- angles that rotate the plate system
  psi = deg2rad(0.d0)
  beta = deg2rad(0.d0)
  gamma = deg2rad(0.d0)
  
  !-- displacement vetor
!--   x0_plate = (/ x0,y0,z0 /)
  R = 0.0
  x0_plate = (/ x0,y0+cos(psi)*R,z0+sin(psi)*R /)
  
  call Rx(M1,psi)
  call Ry(M2,beta)
  call Rz(M3,gamma)
  M_plate = matmul(M1,matmul(M2,M3))
  
  mask = 5000.d0
  
  
    !-- generate some deflection line
  x0=0.d0
  y0=0.d0
  call init_beams( beams )
  do while ((time2<=1.25))
    call SolidSolverWrapper( time2, dt_fixed , beams )
    time2 = time2+dt_fixed
  enddo
  
  beam = beams(1)
  
  
  !-- For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)
      !-- global coordinates
      x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
      !-- in the plate system
      x_plate = matmul(M_plate,x-x0_plate)
      
      !-- loop over points on the beam
      do is = 0,ns-2
        !-- a,b,c: sides of triangle
        a = dsqrt( (x_plate(1)-beam%x(is))**2 + (x_plate(2)-beam%y(is))**2 )
        !-- if the first lagrangian marker is too far away, the second
        !-- one is as well, and we can skip the whole point, greatly
        !-- reducing the computational complexity
        if ( a < 3.d0 ) then
          b = dsqrt( (x_plate(1)-beam%x(is+1))**2 + (x_plate(2)-beam%y(is+1))**2 )
          !-- c is the distance between two markers, thus ds
          c = dsqrt( (beam%x(is)-beam%x(is+1))**2 + (beam%y(is)-beam%y(is+1))**2 )
          !-- angles in the triangle
          alpha = acos ( (b**2+c**2-a**2)/(2.d0*b*c) )
          beta  = acos ( (a**2+c**2-b**2)/(2.d0*a*c) )
          
          !-- height of the triangle: this is what we were looking for!
          h = dsin(alpha)*b
          
          if ((abs(alpha)<= pi/2.0).and.(abs(beta)<=pi/2.d0)) then
            !-- we're in the area of the line segment, where the
            !-- height of the triangle defines the distance function
            mask(ix,iy,iz) = min(h, mask(ix,iy,iz))
          else
            !-- we're in a hinge zone, where it is the closest
            !-- lagrangian marker that defines the distance (circular hinges)
            mask(ix,iy,iz) = min(mask(ix,iy,iz),min(a,b))
          endif
          
        endif
      enddo
      
      
    enddo
   enddo
  enddo
  
  where (mask>1000.d0) mask=0.d0
  
  !-- at the end of thsi process, mask is the distance function from the centerline
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)  
      call smoothstep( tmp, mask(ix,iy,iz)-t_beam, t_beam, 3.d0*dx )
      mask(ix,iy,iz) = tmp
    enddo
   enddo
  enddo  
  
  !-- Calculate unsteady 
end subroutine Draw_flexible_plate
