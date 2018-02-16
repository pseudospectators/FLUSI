module interpolation
  use vars
  use mpi
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

  ! if the point is not on the local grid, return a very large, negative value
  ! (afterwards we take the mpimax)
  value = -9.9d10

  ! if (ng /= 2) call abort(9997110,"linear interp requires 2 ghost nodes for proper periodization")

  if (nx > 1) then
    !********
    ! 3D case
    !********
    ! note border are ra/rb and not ga/gb. However, we allow the lower limit to be one grid point
    ! below our local memory (which means it IS the ghost node). Note on the upper limit, this is not
    ! the case
    ! ascii-scheme:
    ! g o o o o o g
    !  - - - - - -
    ! the dashes are possible interpolation points
    if ((((ix>=ra(1)-1).and.(ix<=rb(1)))).and.&
          (iy>=ra(2)-1).and.(iy<=rb(2)).and.&
          (iz>=ra(3)-1).and.(iz<=rb(3))) then

        c00 = field(ix,iy  ,iz  )*(1.d0-xd)+field(ix+1,iy  ,iz  )*xd
        c10 = field(ix,iy+1,iz  )*(1.d0-xd)+field(ix+1,iy+1,iz  )*xd
        c01 = field(ix,iy  ,iz+1)*(1.d0-xd)+field(ix+1,iy  ,iz+1)*xd
        c11 = field(ix,iy+1,iz+1)*(1.d0-xd)+field(ix+1,iy+1,iz+1)*xd

        c0 = c00*(1.d0-yd) + c10*yd
        c1 = c01*(1.d0-yd) + c11*yd

        value = c0*(1.d0-zd)+c1*zd
    endif
  else
    !********
    ! 2D case
    !********
    if ((iy>=ra(2)-1).and.(iy<=rb(2)+1).and.(iz>=ra(3)-1).and.(iz<=rb(3)+1)) then

        ! ! NOTE: 28 Nov 2016, I discovered the PER function is used here but I cannot remember
        ! ! why we should so so?? the ghost nodes should take care of that.
        ! ! ANSWER: for 2D runs..
        ! c00 = field(per(ix,nx),iy  ,iz  )*(1.d0-xd)+field(per(ix+1,nx),iy  ,iz  )*xd
        ! c10 = field(per(ix,nx),iy+1,iz  )*(1.d0-xd)+field(per(ix+1,nx),iy+1,iz  )*xd
        ! c01 = field(per(ix,nx),iy  ,iz+1)*(1.d0-xd)+field(per(ix+1,nx),iy  ,iz+1)*xd
        ! c11 = field(per(ix,nx),iy+1,iz+1)*(1.d0-xd)+field(per(ix+1,nx),iy+1,iz+1)*xd
        !
        ! c0 = c00*(1.d0-yd) + c10*yd
        ! c1 = c01*(1.d0-yd) + c11*yd
        !
        ! value = c0*(1.d0-zd)+c1*zd
        ! i=int((x_target-x1_box)/dx)  ! attention on index shift because of automatic array
        ! j=int((y_target-y1_box)/dy)
        !
        !
        ! x_1= real(i)*dx + x1_box
        ! y_1= real(j)*dy + y1_box
        ! x_2= dx*real(i+1) + x1_box
        ! y_2= dy*real(j+1) + y1_box
        ! R1 = (x_2-x_target)*field2(i,j)/dx   + (x_target-x_1)*field2(i+1,j)/dx
        ! R2 = (x_2-x_target)*field2(i,j+1)/dx + (x_target-x_1)*field2(i+1,j+1)/dx
        !
        ! LinearInterpolation = (y_2-y_target)*R1/dy + (y_target-y_1)*R2/dy

        call abort(9997111,"linear interpolation for 2d simulations currently not implemented...have fun")
    endif
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
    call abort(67,"Error: number of ghostpoints not suitable for delta interp")
  endif


  !-- note border are ra/rb and not ga/gb
  if ((((ix0>=ra(1)).and.(ix0<=rb(1))).or.(nx==1)).and.&
        (iy0>=ra(2)).and.(iy0<=rb(2)).and.&
        (iz0>=ra(3)).and.(iz0<=rb(3))) then

      value = 0.d0

      do iz=iz0-N_support,iz0+N_support ! the box size around the point
        zz = dble(iz)*dz
        delz = delta(abs(zz - x(3)),dz)
        do iy=iy0-N_support,iy0+N_support
          yy = dble(iy)*dy
          dely = delta(abs(yy - x(2)),dy)
          do ix=ix0-N_support,ix0+N_support
            xx = dble(ix)*dx
            delx = delta(abs(xx - x(1)),dx)
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

  if (r<0.5d0) then
    delta = (3.0d0/8.0d0)+(pi/32.0d0)-0.250d0*r**2
  elseif ( (r>=0.50d0) .and. (r<=1.50d0)   ) then
    delta = 0.25d0 + (1.0d0-r)/8.0d0  *sqrt(-2.0d0 + 8.0d0*r - 4.0d0*r**2) -asin(sqrt(2.0d0)*(r-1.0d0))/8.0d0
  elseif ( (r>=1.5d0) .and. (r<=2.5d0)   ) then
    delta = (17.0d0/16.0d0) - (pi/64.0d0) - (3.0d0*r/4.0d0) + ((r**2)/8.0d0) + &
            (r-2.0d0)*sqrt(-14.0d0 + 16.0d0*r - 4.0d0*r**2)/16.0d0 +asin(sqrt(2.0d0)*(r-2.0d0))/16.0d0
  elseif ( (r>=2.5d0)    ) then
    delta = 0.0d0
  endif

end function


real (kind=pr) function interp2_nonper (x_target, y_target, field2, axis)
!  LINEAR Interpolation in a field. The field is of automatic size, indices starting with 0 both. The domain is
!  defined by x1_box,y1_box and x2_box,y2_box. The target coordinates should lie within that box.
!  NOTE: attention on the upper point of the box. In the rest of the code, which is periodic, the grid is 0:nx-1
!        but the lattice spacing is yl/nx. This means that the point (nx-1) has NOT the coordinate yl but yl-dx
!        (otherwise this point would exist two times!)
  use vars, only : pr
  implicit none
  integer :: i,j
  real (kind=pr) :: x,y,x_1,y_1,x_2,y_2,dx, dy, R1,R2
  real (kind=pr), intent (in) :: field2(0:,0:), x_target, y_target, axis(1:4)
  real(kind=pr) :: x1_box, y1_box, x2_box, y2_box

  x1_box = axis(1)
  x2_box = axis(2)
  y1_box = axis(3)
  y2_box = axis(4)


  dx = (x2_box-x1_box) / dble(size(field2,1)-1 )
  dy = (y2_box-y1_box) / dble(size(field2,2)-1 )


  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    ! return zero if point lies outside valid bounds
    interp2_nonper = 0.0d0
    return
  endif

  i = int((x_target-x1_box)/dx)
  j = int((y_target-y1_box)/dy)

  x_1 = dble(i)*dx + x1_box
  y_1 = dble(j)*dy + y1_box
  x_2 = dx*dble(i+1) + x1_box
  y_2 = dy*dble(j+1) + y1_box
  R1 = (x_2-x_target)*field2(i,j)/dx   + (x_target-x_1)*field2(i+1,j)/dx
  R2 = (x_2-x_target)*field2(i,j+1)/dx + (x_target-x_1)*field2(i+1,j+1)/dx

  interp2_nonper = (y_2-y_target)*R1/dy + (y_target-y_1)*R2/dy

end function interp2_nonper


end module interpolation
