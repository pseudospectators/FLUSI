!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine flexible_wing_motions ( time, wings )
  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), dimension (1:nWings), intent (inout) :: wings
  integer :: i

  do i=1,nWings
  select case(wings(i)%Motion)
  case ("simple_harmonic")
    call simple_harmonic_motion (time, wings(i))
  case ("stationary")
    wings(i) = wings(i)
  end select
  enddo

end subroutine

subroutine simple_harmonic_motion (time, wings)

  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wings
  integer :: j

    do j=1,nVeins_BC
      wings%z_BC(-1,j) = wings%z0 - 0.075/5*sin(5*pi*time)
      wings%z_BC(0,j) = wings%z_BC(-1,j)
    enddo

end subroutine

subroutine translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), dimension (1:nWings), intent(inout) :: wings
integer :: i

do i=1,nWings
    wings(i)%at_inertia(1) = 0.d0
    wings(i)%at_inertia(2) = 0.d0
    wings(i)%at_inertia(3) = - 0*0.075/10*(10*pi)**2*sin(10*pi*time)
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

  wings%x(1:np) = wings%x(1:np) + dt1**2*wings%at_inertia(1)
  wings%y(1:np) = wings%y(1:np) + dt1**2*wings%at_inertia(2)
  wings%z(1:np) = wings%z(1:np) + dt1**2*wings%at_inertia(3)

else

  wings%x(1:np) = wings%x(1:np) + dt1**2*wings%at_inertia(1)
  wings%y(1:np) = wings%y(1:np) + dt1**2*wings%at_inertia(2)
  wings%z(1:np) = wings%z(1:np) + dt1**2*wings%at_inertia(3)

endif

end subroutine
