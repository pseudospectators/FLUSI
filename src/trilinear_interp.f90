subroutine extend_array( field, field_ext )
  use vars
  use mpi 
  implicit none
  real(kind=pr),intent(out) :: field_ext(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)+1)  
  real(kind=pr),intent (in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer ::  mpicode, destination, origin,status

  ! assume 1D decompostion for now
  field_ext(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) = field    
  
  destination = mpirank-1
  if (destination==-1) destination=mpisize-1
  origin = mpirank+1
  if (origin==mpisize) origin=0
  
  call MPI_sendrecv( field_ext(:,:,ra(3)),& ! send buffer
                     nx*ny,& ! send buffer size
                     MPI_DOUBLE_PRECISION,&
                     destination,& ! whom to send it to
                     mpirank,& ! send tag
                     field_ext(:,:,rb(3)+1),&! receive buffer
                     nx*ny,& ! recvcount
                     MPI_DOUBLE_PRECISION,&
                     origin,& !source
                     origin,& ! recv tag
                     MPI_COMM_WORLD,status,mpicode)
                     
end subroutine extend_array

!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! If the point does not lie in the local memory, zero is returned.
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
!-------------------------------------------------------------------------------
subroutine trilinear_interp(x,y,z,field_ext,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),intent (in) :: x,y,z
  real(kind=pr),intent(out) :: value
  
  real(kind=pr),intent(out) :: field_ext(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)+1) 
  
  integer::xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  integer :: ix,iy,iz,mpicode
  
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
  
  if ((ix>=xmin).and.(ix<=xmax).and.&
      (iy>=ymin).and.(iy<=ymax).and.&
      (iz>=zmin).and.(iz<=zmax)) then
      
    c00 = field_ext(ix  ,iy  ,iz  )*(1.d0-xd)+field_ext(ix+1,iy  ,iz  )*xd
    c10 = field_ext(ix  ,iy+1,iz  )*(1.d0-xd)+field_ext(ix+1,iy+1,iz  )*xd
    c01 = field_ext(ix  ,iy  ,iz+1)*(1.d0-xd)+field_ext(ix+1,iy  ,iz+1)*xd
    c11 = field_ext(ix  ,iy+1,iz+1)*(1.d0-xd)+field_ext(ix+1,iy+1,iz+1)*xd
    
    c0 = c00*(1.d0-yd) + c10*yd
    c1 = c01*(1.d0-yd) + c11*yd
    
    value = c0*(1.d0-zd)+c1*zd
  
  else
    !-- point is not on the local grid, return zero
    value = 0.d0
  endif
     
end subroutine trilinear_interp