!-------------------------------------------------------------------------------
! Draws a flexible plate, where the deflection line is solved using the beam 
! solver. 
!-------------------------------------------------------------------------------
subroutine Draw_Plate (time, beam)
  use mpi
  use fsi_vars
  ! use global variables for solid solver, as well as routines
  use solid_model 
  implicit none

  real(kind=pr), intent(in) :: time
  type(Solid), intent(in) :: beam
  real(kind=pr) :: x, y, z
  integer :: ix, iy, iz, mpicode, is
  integer :: ixmin,ixmax,izmin,izmax,iymin,iymax
  
  ! loop over lagrangian points. thus the algorithm is of order
  ! ns complexity. Load balancing is a problem, CPUs not affected
  ! by the solid are finnishing earlier (their cost is not zero, 
  ! though)
  do is = 0, ns-1
    ! for each grid piece, get bounding box (include smoothing)
    ixmin = ...
    ixmax =
    iymin = 
    iymax =
    izmin =
    izmax = ...
    
    ! is this box included in our domain?    
    if (is_included) then
    
    endif
    
  enddo

  ! For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)
     x = dble(ix)*dx
     y = dble(iy)*dy
     z = dble(iz)*dz
     
    enddo
   enddo
  enddo

  ! Calculate unsteady 
end subroutine Draw_Plate
