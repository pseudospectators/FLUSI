module interpolation
  contains

!-------------------------------------------------------------------------------
! delta_interpolation
!-------------------------------------------------------------------------------
subroutine delta_interpolation(x,field,value)
  use vars
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
  if (  ( ((ix0>=ra(1)).and.(ix0<=rb(1))) .or. (nx==1) ) .and. &
           (iy0>=ra(2)).and.(iy0<=rb(2)) .and.&
           (iz0>=ra(3)).and.(iz0<=rb(3)) ) then

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


!-------------------------------------------------------------------------------
! Trilinear interpolation of the field at the location (x,y,z)
! follows closely http://en.wikipedia.org/wiki/Trilinear_interpolation
!-------------------------------------------------------------------------------
! Case 1: one rank has the complete array and wishes to interpolate on it
! periodic=.true. or .false.
! Case 2: the global array is mpi-distributed, so here we get a NON-PERIODIC subset
! of a greater array: periodic=.false.
! Do not forget to specify the origin of the chunk (that is the coordinate of the lower GHOST NODE)
!-------------------------------------------------------------------------------
! The grid looks like this:
!
! x x x x x
! x x x x x
! x x x x x
! x x x x x
! Y x x x x
!
! The point Y=x0 has the indices (0,0,0) in the field
!-------------------------------------------------------------------------------
! The routin returns -9.9d10 if the specified x_target is not on the grid (which
! also happens if this mpirank is not the owner of the chunk of data required for
! interpolation)
!-------------------------------------------------------------------------------
real(kind=pr) function trilinear_interp(x0, dx, field, x_target, periodic)
  use vars, only: pr, abort
  implicit none
  ! origin of array grid. Note this is a compatibility issue with WABBIT. In FLUSI
  ! we usually have only one grid start starts at 0,0,0 and then the variables ra(:) and rb(:)
  ! inidicate which part is on the mpirank. In WABBIT we have bocks, always 1:n,1:n,1:n but
  ! each with its own origin and spacing. For transferability, this routine is written in the
  ! general style. In FLUSI:
  !     x0=(/dble(ra(1))*dx, dble(ra(2))*dy, dble(ra(3))*dz/) and dx=(/dx,dy,dz/)
  ! or, with ghost nodes:
  !     x0=(/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/) and dx=(/dx,dy,dz/)
  real(kind=pr),dimension(1:3),intent(in) :: x0
  ! spacing of array grid
  real(kind=pr),dimension(1:3),intent(in) :: dx
  ! the point at which we wish to interpolate
  real(kind=pr),dimension(1:3),intent(in) :: x_target
  ! actual field. zero-based indexing.
  real(kind=pr),intent(inout) :: field(0:,0:,0:)
  ! assume periodicity of field or not?
  ! ATTENTION: this means we suppose the array FIELDS to be PERIODIC ON ITS OWN
  ! the global field may well be PERIODIC, but if you pass mpi-chunks to this routine
  ! you MUST set periodic=.false. even if the global field is periodic.
  logical, intent(in) :: periodic

  ! array bounds of the field array
  integer,dimension(1:3) :: lbounds, ubounds
  real(kind=pr) :: xd, yd, zd
  real(kind=pr) :: c00, c10, c01, c11, c0, c1
  integer :: ix, iy, iz, ix1, iy1, iz1

  lbounds = 0
  ubounds = (/size(field,1), size(field,2), size(field,3)/) - 1

  ! indices of cube containing the target point, lower end
  ix = floor( (x_target(1)-x0(1))/dx(1) )
  iy = floor( (x_target(2)-x0(2))/dx(2) )
  iz = floor( (x_target(3)-x0(3))/dx(3) )

  ! distance to lower point, normalized (0..1)
  xd = ( x_target(1)-(dble(ix)*dx(1)+x0(1)) ) / dx(1)
  yd = ( x_target(2)-(dble(iy)*dx(2)+x0(2)) ) / dx(2)
  zd = ( x_target(3)-(dble(iz)*dx(3)+x0(3)) ) / dx(3)

  ! if the point is not on the grid, return a very large, negative value
  trilinear_interp = -9.9d10

  if (periodic) then
    ! *** periodic case ***
    if ( (ix>=lbounds(1)).and.(ix<=ubounds(1)) ) then
      if ( (iy>=lbounds(2)).and.(iy<=ubounds(2)) ) then
        if ( (iz>=lbounds(3)).and.(iz<=ubounds(3)) ) then
          ix1 = ix+1
          iy1 = iy+1
          iz1 = iz+1

          ! periodization. note ix,iy,iz can be ubounds at most, so the next point
          ! ix+1 would be the first point again.
          if (ix1>ubounds(1)) ix1 = lbounds(1)
          if (iy1>ubounds(2)) iy1 = lbounds(2)
          if (iz1>ubounds(3)) iz1 = lbounds(3)

          c00 = field(ix,iy  ,iz )*(1.d0-xd)+field(ix1 ,iy  ,iz )*xd
          c10 = field(ix,iy1 ,iz )*(1.d0-xd)+field(ix1 ,iy1 ,iz )*xd
          c01 = field(ix,iy  ,iz1)*(1.d0-xd)+field(ix1 ,iy  ,iz1)*xd
          c11 = field(ix,iy1 ,iz1)*(1.d0-xd)+field(ix1 ,iy1 ,iz1)*xd

          c0 = c00*(1.d0-yd) + c10*yd
          c1 = c01*(1.d0-yd) + c11*yd

          trilinear_interp = c0*(1.d0-zd)+c1*zd
        endif
      endif
    endif
  else
    ! *** non-periodic case ***
    if ( (ix>=lbounds(1)).and.(ix<ubounds(1)) ) then
      if ( (iy>=lbounds(2)).and.(iy<ubounds(2)) ) then
        if ( (iz>=lbounds(3)).and.(iz<ubounds(3)) ) then
          ix1 = ix+1
          iy1 = iy+1
          iz1 = iz+1

          c00 = field(ix,iy  ,iz )*(1.d0-xd)+field(ix1 ,iy  ,iz )*xd
          c10 = field(ix,iy1 ,iz )*(1.d0-xd)+field(ix1 ,iy1 ,iz )*xd
          c01 = field(ix,iy  ,iz1)*(1.d0-xd)+field(ix1 ,iy  ,iz1)*xd
          c11 = field(ix,iy1 ,iz1)*(1.d0-xd)+field(ix1 ,iy1 ,iz1)*xd

          c0 = c00*(1.d0-yd) + c10*yd
          c1 = c01*(1.d0-yd) + c11*yd

          trilinear_interp = c0*(1.d0-zd)+c1*zd
        endif
      endif
    endif
  endif

  ! with periodic it can work, but non-periodic not
  if (size(field,1)<=1) then
    call abort(9997111,"linear interpolation for 2d simulations currently not implemented...have fun")
  end if
end function trilinear_interp




subroutine interp_testing()
  use vars
  real(kind=pr), dimension(1:10,1:10,1:10) :: field
  integer :: ix, iy, iz
  real(kind=pr) :: x,y,z,v,vex

  x0 = 10.d0
  y0 = 12.d0
  z0 = 33.33d0

  dx = 0.70d0
  dy = 0.80d0
  dz = 1.11d0

  do ix =1,10
    do iy =1,10
      do iz =1,10
        x = x0 + dble(ix)*dx
        y = y0 + dble(iy)*dy
        z = z0 + dble(iz)*dz

        field(ix,iy,iz) = 2.0d0*x + 5.0d0*y + 9.0d0*z

      enddo
    enddo
  enddo

  ! do ix =2,3
  !   do iy =2,3
  !     do iz =2,3
  ix = 6
  iy = 9
  iz = 4
        x = x0 + dble(ix)*dx + 0.23d0*dx
        y = y0 + dble(iy)*dy + 0.22d0*dy
        z = z0 + dble(iz)*dz + 0.21d0*dz

        v = trilinear_interp( (/x0,y0,z0/), (/dx,dy,dz/), field, (/x,y,z/), .false. )

        vex = 2.0d0*x + 5.0d0*y + 9.0d0*z

        write(*,*) "error:", vex-v

  !     enddo
  !   enddo
  ! enddo



end subroutine

end module interpolation
