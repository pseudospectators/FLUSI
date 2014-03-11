!-------------------------------------------------------------------------------
! add a channel mask to an existing mask function
! e.g. there is an insect already in mask and we now put
! channel walls around it.
! the channel walls have the color "0" which helps excluding them 
! for example when computing the integral forces
!-------------------------------------------------------------------------------
subroutine add_channel()
  use mpi
  use fsi_vars
  implicit none
  integer :: ix,iy,iz
  real (kind=pr) :: thick_wall, y, z, pos_wall
  
  ! Wall parameters
  thick_wall = 0.2d0
  pos_wall = 0.3d0
  
  
  ! loop over the physical space
  do ix = ra(1), rb(1)
    do iy = ra(2), rb(2)
       do iz = ra(3), rb(3)
          !----------------
          select case (iChannel)
          case ("xz")
            ! Floor - xz solid wall between y_wall-thick_wall and y_wall
            y = dble(iy)*dy
            if ( (y>=pos_wall-thick_wall) .and. (y<=pos_wall) ) then
              mask(ix,iy,iz) = 1.d0
              us(ix,iy,iz,:) = 0.d0  
              ! external boxes have color 0 (important for forces)
              mask_color(ix,iy,iz) = 0
            endif

          case ("xy")
            ! Floor - xy solid wall between z_wall-thick_wall and z_wall
            z = dble(iz)*dz
            if ( (z>=pos_wall-thick_wall) .and. (z<=pos_wall) ) then
              mask(ix,iy,iz) = 1.d0
              us(ix,iy,iz,:) = 0.d0  
              ! external boxes have color 0 (important for forces)
              mask_color(ix,iy,iz) = 0
            endif
            
          case default
            write (*,*) "add_channel()::iChannel is not a known value"
            stop
          end select
          !----------------
       enddo
    enddo
  enddo
 
end subroutine add_channel
