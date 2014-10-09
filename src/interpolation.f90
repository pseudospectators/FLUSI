module interpolation
  use vars
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  contains

!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
! The array is an extended one which contains ghostpoints, call SYNCHRONIZE_GHOSTS
! prior to using this routine
! If the point x lies not on the local proc, a very small number is returned
! because afterwards we will call MPI_ALLREDUCE and take th emax
!-------------------------------------------------------------------------------
subroutine trilinear_interp_ghosts(x,field,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),dimension(1:3),intent (in) :: x
  real(kind=pr),intent(out) :: value  
  real(kind=pr),intent(in) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr)::xd,yd,zd
  real(kind=pr)::c00,c10,c01,c11,c0,c1
  integer :: ix,iy,iz
  
  ix = floor(x(1)/dx)
  iy = floor(x(2)/dy)
  iz = floor(x(3)/dz)

  xd = (x(1)-dble(ix)*dx)/dx
  yd = (x(2)-dble(iy)*dy)/dy
  zd = (x(3)-dble(iz)*dz)/dz
  
  !-- note border are ra/rb and not ga/gb
  if ((((ix>=ra(1)).and.(ix<=rb(1))).or.(nx==1)).and.&
        (iy>=ra(2)).and.(iy<=rb(2)).and.&
        (iz>=ra(3)).and.(iz<=rb(3))) then
  
      c00 = field(per(ix,nx),iy  ,iz  )*(1.d0-xd)+field(per(ix+1,nx),iy  ,iz  )*xd
      c10 = field(per(ix,nx),iy+1,iz  )*(1.d0-xd)+field(per(ix+1,nx),iy+1,iz  )*xd
      c01 = field(per(ix,nx),iy  ,iz+1)*(1.d0-xd)+field(per(ix+1,nx),iy  ,iz+1)*xd
      c11 = field(per(ix,nx),iy+1,iz+1)*(1.d0-xd)+field(per(ix+1,nx),iy+1,iz+1)*xd
      
      c0 = c00*(1.d0-yd) + c10*yd
      c1 = c01*(1.d0-yd) + c11*yd
      
      value = c0*(1.d0-zd)+c1*zd
  else  
      !-- point is not on the local grid, return a very small value 
      !-- (afterwards we take the max)
      value = -9.9d10
  endif     
end subroutine trilinear_interp_ghosts


!-------------------------------------------------------------------------------
! delta_interpolation
!-------------------------------------------------------------------------------
subroutine delta_interpolation(x,field,value)
  use vars
  use mpi
  implicit none
  real(kind=pr),dimension(1:3),intent (in) :: x
  real(kind=pr),intent(out) :: value  
  real(kind=pr),intent(in) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer :: ix,iy,iz, N_support,ix0,iy0,iz0
  real(kind=pr) :: xx,yy,zz, delx,dely,delz
  
  ix0 = nint(x(1)/dx)
  iy0 = nint(x(2)/dy)
  iz0 = nint(x(3)/dz)

  N_support = 3
  
  if ( ng/=N_support ) then
    write(*,*) "Error: number of ghostpoints not suitable for delta interp"
    call abort()
  endif
  
  
  !-- note border are ra/rb and not ga/gb
  if ((((ix0>=ra(1)).and.(ix0<=rb(1))).or.(nx==1)).and.&
        (iy0>=ra(2)).and.(iy0<=rb(2)).and.&
        (iz0>=ra(3)).and.(iz0<=rb(3))) then
      
      value = 0.d0
      
      do ix=ix0-N_support,ix0+N_support ! the box size around the point
        do iy=iy0-N_support,iy0+N_support
          do iz=iz0-N_support,iz0+N_support
          
            xx = dble(ix)*dx
            yy = dble(iy)*dy
            zz = dble(iz)*dz
            
            delx = delta(abs(xx - x(1)),dx)
            dely = delta(abs(yy - x(2)),dy)
            delz = delta(abs(zz - x(3)),dz) 
            
            value = value + delx*dely*delz*field( per(ix,nx),per(iy,ny),per(iz,nz) )
            
          enddo
        enddo
      enddo
  else  
      !-- point is not on the local grid, return a very small value 
      !-- (afterwards we take the max)
      value = -9.9d10
  endif     
end subroutine delta_interpolation




real (kind=pr) function delta(x,dx1)
  use vars
  real(kind=pr), intent(in) :: x,dx1
  real(kind=pr) :: r
  !----------------------------------
  ! This function returns a delta kernel
  ! see Yang, Zhang, Li: A smoothing technique for discrete delta functions [...] JCP 228 (2009)
  !----------------------------------  
  r=abs(x/dx1)  
  
  if (r<0.5) then
    delta = (3./8.)+(pi/32.)-0.25*r**2
  elseif ( (r>=0.5) .and. (r<=1.5)   ) then
    delta = 0.25 + (1-r)/8.  *sqrt(-2.+8.*r-4.*r**2) -asin(sqrt(2.)*(r-1.))/8.
  elseif ( (r>=1.5) .and. (r<=2.5)   ) then
    delta = (17./16.) - (pi/64.) - (3.*r/4.) + ((r**2)/8.) + (r-2.)*sqrt(-14.+16.*r-4.*r**2)/16. +asin(sqrt(2.)*(r-2.))/16.
  elseif ( (r>=2.5)    ) then
    delta = 0.0
  endif  
  
end function

end module interpolation
