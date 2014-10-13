module diff

use vars
use ghosts


 contains

!------------------------------------------------------------------------------- 
 

subroutine diff_1x2cp( u, u_dx )
  ! First derivative in z-direction with 2nd order central finite differences
  ! diff_1z2cp
  !          ^-----periodic
  !         ^------centered
  !        ^-------convergence order
  !       ^--------direction
  !      ^---------derivative (1st, 2nd...)
  
  use vars
  implicit none
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dx
  ! spacing:
  real(Kind=pr) :: factor
  integer::ix,iy,iz

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (2.d0*dx)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dx(ix,iy,iz) = factor*(u(ix+1,iy,iz)-u(ix-1,iy,iz))
      enddo
    enddo
  enddo       
       
end subroutine diff_1x2cp

!-------------------------------------------------------------------------------

subroutine diff_1y2cp( u, u_dy )
  ! First derivative in z-direction with 2nd order central finite differences
  
  use vars
  implicit none
  integer::ix,iy,iz
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dy
  ! spacing:
  real(Kind=pr) :: factor

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (2.d0*dy)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dy(ix,iy,iz) = factor*(u(ix,iy+1,iz)-u(ix,iy-1,iz))
      enddo
    enddo
  enddo   
end subroutine diff_1y2cp

!-------------------------------------------------------------------------------

subroutine diff_1z2cp( u, u_dz )
  ! First derivative in z-direction with 2nd order central finite differences
  
  use vars
  implicit none
  integer::ix,iy,iz
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dz
  ! spacing:
  real(Kind=pr) :: factor

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (2.d0*dz)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dz(ix,iy,iz) = factor*(u(ix,iy,iz+1)-u(ix,iy,iz-1))
      enddo
    enddo
  enddo   
end subroutine diff_1z2cp

!-------------------------------------------------------------------------------

subroutine diff_2x2cp( u, u_dx )
  ! First derivative in z-direction with 2nd order central finite differences
  
  use vars
  implicit none
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dx
  ! spacing:
  real(Kind=pr) :: factor
  integer::ix,iy,iz

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (dx**2)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dx(ix,iy,iz) = factor*(u(ix+1,iy,iz)-2.d0*u(ix,iy,iz)+u(ix-1,iy,iz))
      enddo
    enddo
  enddo       
       
end subroutine diff_2x2cp

!-------------------------------------------------------------------------------

subroutine diff_2y2cp( u, u_dy )
  ! First derivative in z-direction with 2nd order central finite differences
  
  use vars
  implicit none
  integer::ix,iy,iz
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dy
  ! spacing:
  real(Kind=pr) :: factor

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (dy**2)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dy(ix,iy,iz) = factor*(u(ix,iy+1,iz)-2.d0*u(ix,iy,iz)+u(ix,iy-1,iz))
      enddo
    enddo
  enddo   
end subroutine diff_2y2cp

!-------------------------------------------------------------------------------

subroutine diff_2z2cp( u, u_dz )
  ! First derivative in z-direction with 2nd order central finite differences
  
  use vars
  implicit none
  integer::ix,iy,iz
  ! input: the field to be derived
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(inout)::u
  ! output: the du/dz
  real(kind=pr),dimension(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)),intent(out)::u_dz
  ! spacing:
  real(Kind=pr) :: factor

  ! use multiplication since it is faster than division. Note in later versions,
  ! this factor depends on the position because of the non-constant spacing.
  factor = 1.d0 / (dz**2)
  
  ! update ghost points. TODO: this is not always necessary, it should be done
  ! outside of this routine
  call synchronize_ghosts(u)

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        u_dz(ix,iy,iz) = factor*(u(ix,iy,iz+1)-2.d0*u(ix,iy,iz)+u(ix,iy,iz-1))
      enddo
    enddo
  enddo   
end subroutine diff_2z2cp

!-------------------------------------------------------------------------------

end module diff

