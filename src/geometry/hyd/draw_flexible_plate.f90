!-------------------------------------------------------------------------------
! Draws a flexible plate, where the deflection line is solved using the beam 
! solver. 
!-------------------------------------------------------------------------------
subroutine Draw_flexible_plate (time)!, beam)
  use mpi
  use fsi_vars
  ! use global variables for solid solver, as well as routines
  use solid_model 
  implicit none

  real(kind=pr), intent(in) :: time
!   type(Solid), intent(in) :: beam
  real(kind=pr) :: psi, beta, gamma, R
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate,M1,M2,M3
  integer :: ix, iy, iz, mpicode, is
  integer :: ixmin,ixmax,izmin,izmax,iymin,iymax

  !-- angles that rotate the plate system
  psi = deg2rad(45.d0)
  beta = deg2rad(0.d0)
  gamma = deg2rad(0.d0)
  
  !-- displacement vetor
!   x0_plate = (/ x0,y0,z0 /)
  R = 0.5
  x0_plate = (/ x0,y0+cos(psi)*R,z0+sin(psi)*R /)
  
  call Rx(M1,psi)
  call Ry(M2,beta)
  call Rz(M3,gamma)
  M_plate = matmul(M1,matmul(M2,M3))
  
  
  ! For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)
      !-- global coordinates
      x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
      !-- in the plate system
      x_plate = matmul(M_plate,x-x0_plate)
      
      
      !-- let's create a simple rectangle
      if ((x_plate(1) >=0.d0).and.(x_plate(1)<=1.d0)) then
      if ((x_plate(2) >=-0.25d0).and.(x_plate(2)<=0.25d0)) then
      if ((x_plate(3) >=-2.d0*dx).and.(x_plate(3)<=2.d0*dx)) then
      mask(ix,iy,iz)=1.d0
      
      endif
      endif
      endif
    enddo
   enddo
  enddo

  ! Calculate unsteady 
end subroutine Draw_flexible_plate
