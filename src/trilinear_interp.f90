module interpolation
use vars
use mpi 

integer,save :: ixmin,ixmax,iymin,iymax,izmin,izmax

 contains



subroutine extend_array_1D( field, ghosts )
  use vars
  use mpi 
  implicit none 
  real(kind=pr),intent (in) :: field(ixmin:ixmax,iymin:iymax,izmin:izmax)
  real(kind=pr),intent(out) :: ghosts(ixmin:ixmax,iymin:iymax)
  integer ::  mpicode, destination, origin,status

  if (mpidims(2)/=1) then
    write(*,*) "Fail: trying to use extend_array_1D with 2D data decompostion"
    stop
  endif
  
  destination = mpirank-1
  if (destination==-1) destination=mpisize-1
  origin = mpirank+1
  if (origin==mpisize) origin=0
  
  call MPI_sendrecv( field(:,:,izmin),&         ! send buffer
                     nx*ny,&                    ! send buffer size
                     MPI_DOUBLE_PRECISION,&
                     destination,&              ! whom to send it to
                     mpirank,&                  ! send tag
                     ghosts,&                   ! receive buffer
                     nx*ny,&                    ! recvcount
                     MPI_DOUBLE_PRECISION,&
                     origin,&                   !source
                     origin,&                   ! recv tag
                     MPI_COMM_WORLD,status,mpicode)
                     
end subroutine extend_array_1D




subroutine extend_array_2D( field, ghostsz, ghostsy )
  use vars
  use mpi 
  implicit none 
  real(kind=pr),intent (in) :: field(ixmin:ixmax,iymin:iymax,izmin:izmax)
  real(kind=pr),intent(out) :: ghostsz(ixmin:ixmax,iymin:iymax+1) ! both include the annoying part
  real(kind=pr),intent(out) :: ghostsy(ixmin:ixmax,izmin:izmax+1)
  integer ::  mpicode, destination, origin, status
  !-----------------------------------------------------------------------------
  !
  !      I receive:                     I send:
  !
  !           ? ? ? ? ? ? ?             
  ! iymax---- o o o o o o ?             S o o o o o
  !           o o o o o o ?             S o o o o o
  !           o o o o o o ?             S o o o o o
  ! iymin---- o o o o o o ?             S S S S S S
  !           |         |               
  !           izmin     izmax           
  !
  !                                     
  !           B B B B B B *
  ! iymax---- o o o o o o A
  !           o o o o o o A
  !           o o o o o o A
  ! iymin---- o o o o o o A
  !           |         |
  !           izmin     izmax
  !
  ! Step a) (look at the CPU in the middle)
  !
  !           o o o o o o     o o o o o o    o o o o o o
  !           o o o o o o     o o o o o o    o o o o o o      
  !           o o o o o o     o o o o o o    o o o o o o       
  !           o o o o o o     s s s s s s    o o o o o o
  !
  !                           r r r r r r
  !           o o o o o o     o o o o o o    o o o o o o
  !           o o o o o o     o o o o o o    o o o o o o      
  !           o o o o o o     o o o o o o    o o o o o o       
  !           o o o o o o     s s s s s s    o o o o o o
  !
  !                           r r r r r r
  !           o o o o o o     o o o o o o    o o o o o o
  !           o o o o o o     o o o o o o    o o o o o o      
  !           o o o o o o     o o o o o o    o o o o o o       
  !           o o o o o o     o o o o o o    o o o o o o  
  !
  !  
  ! Step b) (look at the CPU in the middle)
  !
  !           o o o o o o       o o o o o o     o o o o o o
  !           o o o o o o       o o o o o o     o o o o o o      
  !           o o o o o o       o o o o o o     o o o o o o       
  !           o o o o o o       o o o o o o     o o o o o o
  !                           
  !           o o o o o o r     s o o o o o r   s o o o o o
  !           o o o o o o r     s o o o o o r   s o o o o o      
  !           o o o o o o r     s o o o o o r   s o o o o o       
  !           o o o o o o r     s o o o o o r   s o o o o o
  !                           
  !           o o o o o o       o o o o o o     o o o o o o
  !           o o o o o o       o o o o o o     o o o o o o      
  !           o o o o o o       o o o o o o     o o o o o o       
  !           o o o o o o       o o o o o o     o o o o o o  
  !
  !  
  ! Step c) (look at the CPU in the middle)
  !
  !           o o o o o o       o o o o o o     o o o o o o
  !           o o o o o o       o o o o o o     o o o o o o      
  !           o o o o o o       o o o o o o     o o o o o o       
  !           o o o o o o       o o o o o o     S o o o o o
  
  !                                        R
  !           o o o o o o       o o o o o o     o o o o o o
  !           o o o o o o       o o o o o o     o o o o o o      
  !           o o o o o o       o o o o o o     o o o o o o       
  !           o o o o o o       S o o o o o     o o o o o o
  
  !                       R     
  !           o o o o o o       o o o o o o     o o o o o o
  !           o o o o o o       o o o o o o     o o o o o o      
  !           o o o o o o       o o o o o o     o o o o o o       
  !           o o o o o o       o o o o o o     o o o o o o  
  !-----------------------------------------------------------------------------
  

  ! first step is extension in z direction (points "A")
  destination = yz_plane_ranks(iymin,GetIndex(izmin-1,nz))
  origin = yz_plane_ranks(iymin,GetIndex(izmax+1,nz))
  
  call MPI_sendrecv( field(:,:,izmin),&         ! send buffer
                     nx*(iymax-iymin+1),&       ! send buffer size
                     MPI_DOUBLE_PRECISION,&
                     destination,&              ! whom to send it to
                     mpirank,&                  ! send tag
                     ghostsz(:,iymin:iymax),&   ! receive buffer, NOT full size
                     nx*(iymax-iymin+1),&       ! recvcount
                     MPI_DOUBLE_PRECISION,&
                     origin,&                   ! source
                     origin,&                   ! recv tag
                     MPI_COMM_WORLD,status,mpicode)
                     
  ! second step is extension in y direction (points "B")
  destination = yz_plane_ranks(GetIndex(iymin-1,ny),izmin)
  origin = yz_plane_ranks(GetIndex(iymax+1,ny),izmin)
  
  call MPI_sendrecv( field(:,iymin,:),&         ! send buffer
                     nx*(izmax-izmin+1),&       ! send buffer size
                     MPI_DOUBLE_PRECISION,&
                     destination,&              ! whom to send it to
                     mpirank,&                  ! send tag
                     ghostsy(:,izmin:izmax),&   ! receive buffer
                     nx*(izmax-izmin+1),&       ! recvcount
                     MPI_DOUBLE_PRECISION,&
                     origin,&                   ! source
                     origin,&                   ! recv tag
                     MPI_COMM_WORLD,status,mpicode)
                     
  ! third step is the annoying line that lies on the diagonal neigbor (point "*")
  destination = yz_plane_ranks(GetIndex(iymin-1,ny),GetIndex(izmin-1,nz))
  origin = yz_plane_ranks(GetIndex(iymax+1,ny),GetIndex(izmax+1,nz))
  
  call MPI_sendrecv( field(:,iymin,izmin),&     ! send buffer
                     nx,&                       ! send buffer size (this is a line here)
                     MPI_DOUBLE_PRECISION,&
                     destination,&              ! whom to send it to
                     mpirank,&                  ! send tag
                     ghostsy(:,izmax+1),&       ! receive buffer
                     nx,&                       ! recvcount
                     MPI_DOUBLE_PRECISION,&
                     origin,&                   ! source
                     origin,&                   ! recv tag
                     MPI_COMM_WORLD,status,mpicode)      
  ! copy anoying part to both arrays (this is redundant)                   
  ghostsz(:,iymax+1) =  ghostsy(:,izmax+1)                    
end subroutine extend_array_2D


!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! If the point does not lie in the local memory, zero is returned.
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
! Version for 1D domain decompostion
!-------------------------------------------------------------------------------
subroutine trilinear_interp_1Ddecomp(x,field,ghosts,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),dimension(1:3),intent (in) :: x
  real(kind=pr),intent(out) :: value  
  real(kind=pr),intent(in) :: field(ixmin:ixmax,iymin:iymax,izmin:izmax) 
  real(kind=pr),intent(in) :: ghosts(ixmin:ixmax,iymin:iymax)  
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  integer :: ix,iy,iz
  
  ix = floor(x(1)/dx)
  iy = floor(x(2)/dy)
  iz = floor(x(3)/dz)
  
  xd = (x(1)-dble(ix)*dx)/dx
  yd = (x(2)-dble(iy)*dy)/dy
  zd = (x(3)-dble(iz)*dz)/dz
  
  !-- point lies in the interior of the local array and we do not need the 
  !-- ghostpoints
  if ((ix>=ixmin).and.(ix<=ixmax).and.&
      (iy>=iymin).and.(iy<=iymax).and.&
      (iz>=izmin).and.(iz< izmax)) then
      
    c00 = field(ix  ,iy  ,iz  )*(1.d0-xd)+field(ix+1,iy  ,iz  )*xd
    c10 = field(ix  ,iy+1,iz  )*(1.d0-xd)+field(ix+1,iy+1,iz  )*xd
    c01 = field(ix  ,iy  ,iz+1)*(1.d0-xd)+field(ix+1,iy  ,iz+1)*xd
    c11 = field(ix  ,iy+1,iz+1)*(1.d0-xd)+field(ix+1,iy+1,iz+1)*xd
    
    c0 = c00*(1.d0-yd) + c10*yd
    c1 = c01*(1.d0-yd) + c11*yd
    
    value = c0*(1.d0-zd)+c1*zd
  
  !-- the lower z-index is the largest on the local storage, thus the next point
  !-- iz+1 is on another processor, which we synchronized previsously
  elseif ((ix>=ixmin).and.(ix<=ixmax).and.&
          (iy>=iymin).and.(iy<=iymax).and.&
          (iz>=izmin).and.(iz==izmax)) then
          
    c00 = field(ix  ,iy  ,iz  )*(1.d0-xd)+field(ix+1,iy  ,iz  )*xd
    c10 = field(ix  ,iy+1,iz  )*(1.d0-xd)+field(ix+1,iy+1,iz  )*xd
    c01 = ghosts(ix  ,iy  )*(1.d0-xd)+ghosts(ix+1,iy  )*xd
    c11 = ghosts(ix  ,iy+1)*(1.d0-xd)+ghosts(ix+1,iy+1)*xd
    
    c0 = c00*(1.d0-yd) + c10*yd
    c1 = c01*(1.d0-yd) + c11*yd
    
    value = c0*(1.d0-zd)+c1*zd
      
  !-- point is not on the local grid, return zero
  else    
    value = 0.d0
  endif
     
end subroutine trilinear_interp_1Ddecomp


!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! If the point does not lie in the local memory, zero is returned.
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
! Version for 2D domain decompostion
!-------------------------------------------------------------------------------
subroutine trilinear_interp_2Ddecomp(x,field,ghostsz,ghostsy,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),dimension(1:3),intent (in) :: x
  real(kind=pr),intent(out) :: value  
  real(kind=pr),intent(in) :: field(ixmin:ixmax,iymin:iymax,izmin:izmax) 
  real(kind=pr),intent(in) :: ghostsz(ixmin:ixmax,iymin:iymax+1)  
  real(kind=pr),intent(in) :: ghostsy(ixmin:ixmax,izmin:izmax+1)  
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  real(kind=pr)::c000,c100,c001,c101,c110,c111,c011,c010
  integer :: ix,iy,iz
  
  ix = floor(x(1)/dx)
  iy = floor(x(2)/dy)
  iz = floor(x(3)/dz)
  
  xd = (x(1)-dble(ix)*dx)/dx
  yd = (x(2)-dble(iy)*dy)/dy
  zd = (x(3)-dble(iz)*dz)/dz
  
  !-- point lies in the interior of the local array and we do not need the 
  !-- ghostpoints
  if ((ix>=ixmin).and.(ix< ixmax).and.&
      (iy>=iymin).and.(iy< iymax).and.&
      (iz>=izmin).and.(iz< izmax)) then
      
    c000 = field(ix  ,iy  ,iz  )
    c100 = field(ix+1,iy  ,iz  )
    c010 = field(ix  ,iy+1,iz  )
    c110 = field(ix+1,iy+1,iz  )
    c001 = field(ix  ,iy  ,iz+1)
    c101 = field(ix+1,iy  ,iz+1)
    c011 = field(ix  ,iy+1,iz+1)
    c111 = field(ix+1,iy+1,iz+1)
  
  elseif ((ix>=ixmin).and.(ix< ixmax).and.&
          (iy>=iymin).and.(iy< iymax).and.&
          (iz>=izmin).and.(iz==izmax)) then

    c000 = field(ix  ,iy  ,iz  )
    c100 = field(ix+1,iy  ,iz  )
    c010 = field(ix  ,iy+1,iz  )
    c110 = field(ix+1,iy+1,iz  )
    c001 = ghostsz(ix  ,iy  )
    c101 = ghostsz(ix+1,iy  )
    c011 = ghostsz(ix  ,iy+1)
    c111 = ghostsz(ix+1,iy+1)
          
  elseif ((ix>=ixmin).and.(ix< ixmax).and.&
          (iy>=iymin).and.(iy==iymax).and.&
          (iz>=izmin).and.(iz< izmax)) then
 
    c000 = field(ix  ,iy  ,iz  )
    c100 = field(ix+1,iy  ,iz  )
    c010 = ghostsy(ix,iz)
    c110 = ghostsy(ix+1,iz)
    c001 = field(ix  ,iy  ,iz+1)
    c101 = field(ix+1,iy  ,iz+1)
    c011 = ghostsy(ix,iz+1)
    c111 = ghostsy(ix+1,iz+1)
 
  elseif ((ix>=ixmin).and.(ix< ixmax).and.&
          (iy>=iymin).and.(iy==iymax).and.&
          (iz>=izmin).and.(iz==izmax)) then
          
    c000 = field(ix  ,iy  ,iz  )
    c100 = field(ix+1,iy  ,iz  )
    c010 = ghostsy(ix,iz)
    c110 = ghostsy(ix+1,iz)
    c001 = ghostsz(ix,iy)
    c101 = ghostsz(ix+1,iy)
    c011 = ghostsy(ix,iz+1)
    c111 = ghostsz(ix+1,iy+1) ! or ghostsy(ix+1.iz+1)
          
  !-- point is not on the local grid, return zero
  else    
    c000 = 0.d0
    c100 = 0.d0
    c010 = 0.d0
    c110 = 0.d0
    c001 = 0.d0
    c101 = 0.d0
    c011 = 0.d0
    c111 = 0.d0
  endif
  
  c00 = c000*(1.d0-xd)+c100*xd
  c10 = c010*(1.d0-xd)+c110*xd
  c01 = c001*(1.d0-xd)+c101*xd
  c11 = c011*(1.d0-xd)+c111*xd
  
  c0 = c00*(1.d0-yd) + c10*yd
  c1 = c01*(1.d0-yd) + c11*yd
  value = c0*(1.d0-zd)+c1*zd
end subroutine trilinear_interp_2Ddecomp



!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! If the point does not lie in the local memory, zero is returned.
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
! Serial version (can also be used in MPI when being sure that the point is
! inside the memory)
!-------------------------------------------------------------------------------
subroutine trilinear_interp(x,field,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),dimension(1:3),intent (in) :: x
  real(kind=pr),intent(out) :: value  
  real(kind=pr),intent(in) :: field(ixmin:ixmax,iymin:iymax,izmin:izmax) 
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  integer :: ix,iy,iz
  value = 0.d0
  
  ix = floor(x(1)/dx)
  iy = floor(x(2)/dy)
  iz = floor(x(3)/dz)
  
  xd = (x(1)-dble(ix)*dx)/dx
  yd = (x(2)-dble(iy)*dy)/dy
  zd = (x(3)-dble(iz)*dz)/dz
  
  c00 = field(ix  ,iy  ,iz  )*(1.d0-xd)+field(ix+1,iy  ,iz  )*xd
  c10 = field(ix  ,iy+1,iz  )*(1.d0-xd)+field(ix+1,iy+1,iz  )*xd
  c01 = field(ix  ,iy  ,iz+1)*(1.d0-xd)+field(ix+1,iy  ,iz+1)*xd
  c11 = field(ix  ,iy+1,iz+1)*(1.d0-xd)+field(ix+1,iy+1,iz+1)*xd
  
  c0 = c00*(1.d0-yd) + c10*yd
  c1 = c01*(1.d0-yd) + c11*yd
  
  value = c0*(1.d0-zd)+c1*zd
  
  if (dabs(value)>2.d0) then
    write(*,*) ix,iy,iz
    write(*,*) xd,yd,zd
    write(*,*) ixmin,iymin,izmin
    write(*,*) ixmax,iymax,izmax
    write(*,*) value
    write(*,*) c00, c10, c01, c11
    write(*,*) c0,c1
  endif
  
end subroutine trilinear_interp



subroutine init_interpolation()
!-- dealing with this kind of stuff is very difficult, so we employ a simpler
! notation
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)
end subroutine

end module interpolation
