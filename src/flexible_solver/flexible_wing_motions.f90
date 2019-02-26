!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine flexible_wing_motions ( time, wing )
  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: i

  select case(wing%Motion)
  case ("simple_harmonic")
    call simple_harmonic_motion (time, wing)
  case ("stationary")
    continue
  case ("revolving_wing")
    continue
  case ("harmonic_ocsillation")
    wing%vt0(1) = 0.d0
    wing%vt0(2) = 0.25*(1*pi)*cos(1*pi*time)
    wing%vt0(3) = 0.d0

    wing%at0(1) = 0.d0
    wing%at0(2) = - 0.25*(1*pi)**2*sin(1*pi*time)
    wing%at0(3) = 0.d0
  end select

end subroutine

subroutine revolving_Zimmerman (time, wing)

  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: j
  real(kind=pr) :: phi_y, tau
  real(kind=pr),dimension(1:3,1:3) :: mat_Ry
  real(kind=pr),dimension(1:3) :: u

  tau = 0.4
  phi_y = 2*pi*(tau*dexp(-time/tau) + time) - 2*pi*tau

  call Ry(mat_Ry,phi_y)

    do j=1,nVeins_BC
      u = matmul(mat_Ry,(/wing%x0_BC(0,j) - wing%x0, &
                          wing%y0_BC(0,j) - wing%y0, &
                          wing%z0_BC(0,j) - wing%z0/))
      wing%x_BC(0,j) = u(1) + wing%x0
      wing%y_BC(0,j) = u(2) + wing%y0
      wing%z_BC(0,j) = u(3) + wing%z0

      u = matmul(mat_Ry,(/wing%x0_BC(-1,j) - wing%x0, &
                      wing%y0_BC(-1,j) - wing%y0, &
                      wing%z0_BC(-1,j) - wing%z0/))
      wing%x_BC(-1,j) = u(1) + wing%x0
      wing%y_BC(-1,j) = u(2) + wing%y0
      wing%z_BC(-1,j) = u(3) + wing%z0
    enddo

end subroutine

subroutine simple_harmonic_motion (time, wing)

  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: j

    do j=1,nVeins_BC
      wing%z_BC(-1,j) = wing%z0_BC(-1,j) + 0.25*sin(1*pi*time)
      wing%z_BC(0,j) = wing%z0_BC(0,j) + 0.25*sin(1*pi*time)
    enddo

end subroutine

subroutine translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), dimension (1:nWings), intent(inout) :: wings
integer :: i

do i=1,nWings
    wings(i)%at0(1) = 0.d0
    wings(i)%at0(2) = - 0.25*(1*pi)**2*sin(1*pi*time)
    wings(i)%at0(3) = 0.d0
enddo

end subroutine

subroutine moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), intent(inout) :: wings
real(kind=pr) :: c1, c2, c3, r
integer :: np

np = wings%np

!Calculate scheme coefficients for variable-step BDF2 scheme
r = dt1/dt0
c1=(1+r)**2/(1+2*r)
c2=r**2/(1+2*r)
c3=(1+r)/(1+2*r)

if (it==0) then

  wings%x(1:np) = wings%x(1:np) + dt1*wings%vt0(1) !dt1**2*wings%at0(1)
  wings%y(1:np) = wings%y(1:np) + dt1*wings%vt0(2)
  wings%z(1:np) = wings%z(1:np) + dt1*wings%vt0(3)

else

  wings%x(1:np) = wings%x(1:np) + dt1*wings%vt0(1) !dt1**2*wings%at0(1)
  wings%y(1:np) = wings%y(1:np) + dt1*wings%vt0(2)
  wings%z(1:np) = wings%z(1:np) + dt1*wings%vt0(3)

endif

end subroutine
