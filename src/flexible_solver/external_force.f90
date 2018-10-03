!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct external force vector consists of gravity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_forces_construction(Wings)
! This is actually just for 1 wing

implicit none

type(Wing),dimension(1:nWings), intent(inout)  :: Wings
integer :: i, j, np

do i=1,nWings

  np = wings(i)%np
  ! Initialize
  wings(i)%Fext = 0.d0

  do j=1,np
    wings(i)%Fext(1:np)        = grav(1)*wings(i)%m(j) !forces on the x-direction
    wings(i)%Fext(np+1:2*np)   = grav(2)*wings(i)%m(j) !forces on the y-direction
    wings(i)%Fext(2*np+1:3*np) = grav(3)*wings(i)%m(j) !forces on the z-direction
  enddo

enddo

end subroutine
