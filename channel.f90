subroutine channel ( chi, ix, iy, iz )
!---------------------------------------------------------------
!     Computes penalization term due to channel walls
!---------------------------------------------------------------
  use fsi_vars
  implicit none

  integer, intent (in) :: ix, iy, iz
  real (kind=pr), intent (out) :: chi
  real (kind=pr) :: thick_wall, y, z, pos_wall, epsinv
  character(len=80), save :: iWall

  chi = 0.0d0    ! Initialize chanel mask
  epsinv = 1.0d0/eps    ! Inverse of penalization parameter

  ! Wall parameters
  iWall = "xy"
  thick_wall = 0.2d0
  pos_wall = 0.3d0

  select case (iWall)
  case ("xz")
  ! Floor - xz solid wall between y_wall-thick_wall and y_wall
  y = dble(iy)*dy
  if ( (y>=pos_wall-thick_wall) .and. (y<=pos_wall) ) then
     chi = epsinv
  endif

  case ("xy")
  ! Floor - xy solid wall between z_wall-thick_wall and z_wall
  z = dble(iz)*dz
  if ( (z>=pos_wall-thick_wall) .and. (z<=pos_wall) ) then
     chi = epsinv
  endif
  end select

end subroutine channel
