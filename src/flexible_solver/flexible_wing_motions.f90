!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine flexible_wing_motions ( time, wings )
  implicit none

  real(kind=pr),intent(in) :: time
  type(wing), dimension (1:nWings), intent (inout) :: wings
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
  type(wing), intent (inout) :: wings
  integer :: j

    do j=1,nVeins_BC
      wings%z_BC(-1,j) = wings%z0 + 0.075/10*sin(10*pi*time)
      wings%z_BC(0,j) = wings%z_BC(-1,j)
    enddo


end subroutine
