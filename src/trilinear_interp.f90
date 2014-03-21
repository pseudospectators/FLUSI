!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! If the point does not lie in the local memory, zero is returned.
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
!-------------------------------------------------------------------------------
subroutine trilinear_interp(x,y,z,field,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),intent (in) :: x,y,z
  real(kind=pr),intent(out) :: value
  real(kind=pr),intent (in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr)::xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  integer :: ix,iy,iz
  
  xmin = ra(1)
  ymin = ra(2)
  zmin = ra(3)
  
  xmax = rb(1)
  ymax = rb(2)
  zmax = rb(3)
  
  ix = floor(x/dx)
  iy = floor(y/dy)
  iz = floor(z/dz)
  
  xd = (x-dble(ix)*dx)/dx
  yd = (y-dble(iy)*dy)/dy
  zd = (z-dble(iz)*dz)/dz
  
  if ((ix>=xmin).and.(ix<=xmax-1).and.&
      (iy>=ymin).and.(iy<=ymax-1).and.&
      (iz>=zmin).and.(iz<=zmax-1)) then
      
    c00 = field(ix  ,iy  ,iz  )*(1.d0-xd)+field(ix+1,iy  ,iz  )*xd
    c10 = field(ix  ,iy+1,iz  )*(1.d0-xd)+field(ix+1,iy+1,iz  )*xd
    c01 = field(ix  ,iy  ,iz+1)*(1.d0-xd)+field(ix+1,iy  ,iz+1)*xd
    c11 = field(ix  ,iy+1,iz+1)*(1.d0-xd)+field(ix+1,iy+1,iz+1)*xd
    
    c0 = c00*(1.d0-yd) + c10*yd
    c1 = c01*(1.d0-yd) + c11*yd
    
    value = c0*(1.d0-zd)+c1*zd
!   elseif
  
  else
    !-- point is not on the local grid, return zero
    value = 0.d0
  endif
     
  
end subroutine trilinear_interp