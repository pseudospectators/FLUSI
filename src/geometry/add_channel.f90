!-------------------------------------------------------------------------------
! add a channel mask to an existing mask function
! e.g. there is an insect already in mask and we now put
! channel walls around it.
! the channel walls have the color "0" which helps excluding them 
! for example when computing the integral forces
!-------------------------------------------------------------------------------
subroutine add_channel(mask,mask_color,us)
  use vars
  implicit none
  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
    
  integer :: ix,iy,iz
  real (kind=pr) :: x,y, z, usponge, H_eff, z_chan
  
  
  ! loop over the physical space
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
       do ix = ra(1), rb(1)
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
            
          case ("turek")                  !      z
            ! Turek walls:                !      ^
            thick_wall = 0.2577143d0      !      |
            usponge = 0.5060014d0         !      |
            H_eff = zl-2.d0*thick_wall    !      -------> y
                      
            z = dble(iz)*dz
            y = dble(iy)*dy
            z_chan = z-thick_Wall
            
            !-- channel walls
            if ((z<=thick_Wall).or.(z>=zl-thick_Wall)) then
              mask(ix,iy,iz) = 1.d0
              mask_color(ix,iy,iz) = 0
              us(ix,iy,iz,:) = 0.d0  
            endif
              
            !-- velocity sponge (sharp, maybe smooth it to one side)  
            if ((y<=usponge).and.(z>=thick_Wall).and.(z<=zl-thick_Wall)) then
              mask(ix,iy,iz) = 1.d0
              mask_color(ix,iy,iz) = 0
              ! note in 2D flows, we set nx=1 and run in the y-z plane where
              ! y is the axial direction
              us(ix,iy,iz,1) = 0.d0
              us(ix,iy,iz,2) = 1.5*z_chan*(H_eff-z_chan)/((0.5*H_eff)**2)              
              us(ix,iy,iz,3) = 0.d0
            endif
            
          case default
            write (*,*) "add_channel()::iChannel is not a known value"
            call abort()
          end select
          !----------------
       enddo
    enddo
  enddo
 
end subroutine add_channel
