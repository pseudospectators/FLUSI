!-------------------------------------------------------------------------------
! get plate geometry.
! the plate has an upper and lower surface, which are both 2D arrays of 3D points
! they are described in the plate coordinate system, in which their x,y coordinates
! are given by the beam, and the z coordinate is the height.
!
!-------------------------------------------------------------------------------
subroutine plate_geometry( beam, surfaces, active_points, nh, M_plate, x0_plate )
  use mpi
  use fsi_vars
  implicit none
  
  real(kind=pr),dimension(:,:,:,:),allocatable,intent(out) :: surfaces
  real(kind=pr),dimension(:,:),allocatable,intent(out) :: active_points
  real(kind=pr),dimension(1:3),intent(in) :: x0_plate
  real(kind=pr),dimension(1:3,1:3),intent(in) :: M_plate
  integer, intent(out) :: nh
  type(solid),intent (in) :: beam
  
  integer :: ih,is
  real(kind=pr) :: dh,xf,yf,s,ztop,zbot
  real(kind=pr),dimension(:),allocatable :: heights
  real(kind=pr),dimension(1:3) :: x, x_plate
  
  !-----------------------------------------------------------------------------
  !-- decide how many points in the rigid direction
  !-----------------------------------------------------------------------------
  if (nx>1) then
    !-- "true" 3D case
    !-- number of interpolation points in rigid direction (span)
    nh = nint( L_span/min(dx,dy,dz)  )
    !-- spacing in span direction (which is actually min(dx,dy,dz), you checked)
    dh = L_span/dble(nh)
    allocate(heights(0:nh))
    !-- fill heights array
    do ih=0,nh
      heights(ih) = dble(ih)*dh -0.5d0*L_span
    enddo
    
  elseif (nx==1) then
    !-- 2D case
    nh = 0
    allocate(heights(0:nh))
    !-- in 2D case, just interpolate always the same position (subsequent avg
    !-- leaves value untouched)
    heights = 0.d0
  endif
  
  
  ! this array holds the interpolation points (2 2D arrays of 3D vectors = 4 indices)
  allocate(surfaces(0:ns-1,0:nh,1:2,1:3))
  allocate(active_points(0:ns-1,0:nh))
  surfaces = 17.d0  
  active_points = 1.d0
  
  
  
  do is=0,ns-1
    s = dble(is)*ds
    ztop = +z_top(s) 
    zbot = -z_bottom(s)
    
    do ih=0,nh
      if ((heights(ih)>=zbot).and.(heights(ih)<=ztop)) then
        active_points(is,ih) = 1.d0
      else
        active_points(is,ih) = 0.d0
      endif
    enddo
  enddo
  
  !-----------------------------------------------------------------------------
  ! Compute the two 2D surfaces in 3D space (top and bottom)
  !-----------------------------------------------------------------------------
  do is=0,ns-1
    do ih=0,nh
      !-- top surface points
      xf = beam%x(is)-t_beam*dsin(beam%theta(is))
      yf = beam%y(is)+t_beam*dcos(beam%theta(is))    
      !-- point in plate coordinate system
      x_plate = (/ xf,yf,heights(ih)/)
      x = matmul( transpose(M_plate) , x_plate )
      !-- point in global coordinate system
      x = x + x0_plate
      surfaces(is,ih,1,1:3) = x

      
      !-- bottom surface points
      xf = beam%x(is)+t_beam*dsin(beam%theta(is))
      yf = beam%y(is)-t_beam*dcos(beam%theta(is))   
      !-- point in plate coordinate system
      x_plate = (/ xf,yf,heights(ih)/)
      x = matmul( transpose(M_plate) , x_plate )
      !-- point in global coordinate system
      x = x + x0_plate
      surfaces(is,ih,2,1:3) = x    
      
      if (root) then
      if((x(1)<-dx).or.(x(1)>xl+dx).or.(x(2)<-dy).or.(x(2)>yl+dy).or.(x(3)<-dz).or.(x(3)>zl+dz)) then
        write(*,*) "ERROR: Surface constuction yielded values outside valid bound."
        call abort()
      endif
      endif
    enddo
  enddo
  
  
  
  deallocate(heights)
end subroutine plate_geometry


!-------------------------------------------------------------------------------
! shape functions of the plate
! Note the maximum spanwise length is L_span for both parts,
! so each function has to provide z values smaller 0.5*L_span
! (otherwise, there is a conflict with the bounding boxes)
!-------------------------------------------------------------------------------
real(kind=pr) function z_top(s)
  use fsi_vars
  implicit none
  real(kind=pr),intent(in)::s
  
  select case(plate_shape)
  case("rectangular")
    z_top = 0.5d0*L_span
  case("fish")
    z_top = 0.25*L_span + 0.25*L_span*s**2
  case("half_circle")
!     z_top = 
!     z_top = 0.25d0*L_span * (s*(1.d0-s)/0.25) + 00.25*L_span
    z_top = (1.d0-s)* L_span/4.0 !+ L_span/4.0
  case default
    if(root) write(*,*) "ERROR. plate_shape not defined.."//plate_shape
    call abort()
  end select
end function 



real(kind=pr) function z_bottom(s)
  use fsi_vars
  implicit none
  real(kind=pr),intent(inout)::s
  
  select case(plate_shape)
  case("rectangular")
    z_bottom = 0.5d0*L_span
  case("fish")
    z_bottom = 0.25*L_span + 0.25*L_span*s**2
  case("half_circle")
    z_bottom = max(min((1.d0-s),1.d0),0.d0) * L_span/2.0
    z_bottom = (1.d0-s)* L_span/4.0 !+ L_span/4.0
  case default
    if(root) write(*,*) "ERROR. plate_shape not defined.."//plate_shape
    call abort()
  end select
end function 