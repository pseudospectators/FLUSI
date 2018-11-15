!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct external force vector consists of gravity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_forces_construction(time,dt0,dt1, it,Wings)
! This is actually just for 1 wing

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing),dimension(1:nWings), intent(inout)  :: Wings
integer :: i, j, np

do i=1,nWings

  np = wings(i)%np
  ! Initialize
  wings(i)%Fext = 0.d0

  ! Gravitational forces
  call gravitational_forces (wings(i))

  ! Fictitious forces appear when the reference frames of wings are non-inertial frame
  call fictitious_forces_of_moving_reference_frame (time,dt0,dt1,it,wings(i))


enddo

end subroutine

subroutine gravitational_forces (wings)

implicit none

type(flexible_wing), intent (inout) :: wings
integer :: j, np

  np = wings%np

  do j=1,np
    wings%Fext(1:np)        = grav(1)*wings%m(j) !forces on the x-direction
    wings%Fext(np+1:2*np)   = grav(2)*wings%m(j) !forces on the y-direction
    wings%Fext(2*np+1:3*np) = grav(3)*wings%m(j) !forces on the z-direction
  enddo

end subroutine

subroutine fictitious_forces_of_moving_reference_frame (time,dt0,dt1,it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), intent (inout) :: wings
integer :: i,j, np

  np = wings%np

  do j=1,np
    wings%Fext(j)        = wings%Fext(j) - wings%at_inertia(1)*wings%m(j)  !forces on the x-direction
    wings%Fext(j + np)   = wings%Fext(j + np) - wings%at_inertia(2)*wings%m(j) !forces on the y-direction
    wings%Fext(j + 2*np) = wings%Fext(j + 2*np) - wings%at_inertia(3)*wings%m(j) !forces on the z-direction
  enddo

  !do i=1,nVeins_BC
  !  j = nint(wings%veins_bending_BC(1,2,i))
  !  if (it==0) write(*,*) j
  !    wings%Fext(j)        = wings%Fext(j) - wings%at_inertia(1)*wings%m(j)  !forces on the x-direction
  !    wings%Fext(j + np)   = wings%Fext(j + np) - wings%at_inertia(2)*wings%m(j) !forces on the y-direction
  !    wings%Fext(j + 2*np) = wings%Fext(j + 2*np) - wings%at_inertia(3)*wings%m(j) !forces on the z-direction

  !enddo


end subroutine
