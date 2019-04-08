subroutine length_calculation_wrapper(x,y,z,springs)
!This subroutine calculates actually for the whole vein or membrane

implicit none

real(kind=pr), intent(in) :: x(:),y(:),z(:)
real(kind=pr), intent(inout) :: springs(:,:)
integer :: ind

do ind=1,nint(maxval(springs(:,1)))

    call length_calculation(x(nint(springs(ind,2))), &
                            x(nint(springs(ind,3))), &
                            y(nint(springs(ind,2))), &
                            y(nint(springs(ind,3))), &
                            z(nint(springs(ind,2))), &
                            z(nint(springs(ind,3))), &
                            springs(ind,5))

enddo

end subroutine

subroutine length_calculation(x1,x2,y1,y2,z1,z2,length)

  real(kind=pr), intent(in) :: x1,x2,y1,y2,z1,z2
  real(kind=pr), intent(inout) :: length

  length = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

end subroutine

subroutine angle_calculation_wrapper(x,y,z,vein)
!Calculate relative angles between two consecutive segments belonging to a vein
!   x(vein(ind,2))          x(vein(ind,3))                x(vein(ind,4))
!        O-----------------------O----------------------------O
!                  X1                          X2
!This subroutine calculates actually for the whole vein

implicit none

real(kind=pr), intent(in) :: x(:),y(:),z(:)
real(kind=pr), intent(inout) :: vein(:,:)
!real(kind=pr) :: x1,x2,y1,y2,z1,z2
integer :: ind

do ind=1,nint(maxval(vein(:,1)))

    call angle_calculation(x(nint(vein(ind,2))),x(nint(vein(ind,3))), &
                           x(nint(vein(ind,4))),y(nint(vein(ind,2))), &
                           y(nint(vein(ind,3))),y(nint(vein(ind,4))), &
                           z(nint(vein(ind,2))),z(nint(vein(ind,3))), &
                           z(nint(vein(ind,4))),vein(ind,7),vein(ind,8))

enddo

end subroutine

subroutine angle_calculation(x1,x2,x3,y1,y2,y3,z1,z2,z3,theta,phi)

implicit none

  real(kind=pr), intent(in) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
  real(kind=pr), intent(inout) :: theta, phi
  real(kind=pr) :: dx1, dy1, dz1, dx2, dy2, dz2

  dx1 = x2 - x1
  dy1 = y2 - y1
  dz1 = z2 - z1
  dx2 = x3 - x2
  dy2 = y3 - y2
  dz2 = z3 - z2

  theta = atan2(dx1*dy2-dx2*dy1,dx1*dx2+dy1*dy2)

    phi = atan2(dx1*dz2-dx2*dz1,dx1*dx2+dz1*dz2)

end subroutine

subroutine convert_flexural_rigidity_into_spring_stiffness(EIy, EIz, kby0, kbz0, extension_springs)

  implicit none
  real(kind=pr), intent(in) :: extension_springs(1:,:)
  real(kind=pr), intent(in) :: EIy, EIz
  real(kind=pr), intent(inout) :: kby0, kbz0
  real(kind=pr) :: L, ns

  ns = maxval(extension_springs(:,1)) + 1
  L = sum(extension_springs(1:(nint(ns)-1),4))

  kby0 = EIz/L*ns*(2*ns+1)/(2*ns+2)
  kbz0 = EIy/L*ns*(2*ns+1)/(2*ns+2)

end subroutine

subroutine solve_linear_system_wing ( A, b, x )
  implicit none
  real(kind=pr),dimension(1:,1:), intent(in) :: A
  real(kind=pr),dimension(1:), intent(inout) :: x
  real(kind=pr),dimension(1:), intent(in) :: b
  real(kind=pr) :: t0
  integer :: error, nn
  integer, allocatable, dimension(:) :: ipiv

  nn = size(b)

  if (size(b) /= size(x) ) call abort(91,'Solve_linear_system: size(b)/=size(x)')

  allocate(ipiv(1:nn))

  t0 = MPI_wtime()
  
  call dgetrf( nn, nn, A , nn, ipiv, error )
  if (error .ne. 0) then
    write(*,*) "!!! mmCrutial: dgetrf error.", error ,nn
    call abort(92, "Error in solve liner system (dgetrf)")
  endif

  x = b
  call dgetrs( 'N', nn, 1, A, nn, ipiv, x, nn, error )
  if (error .ne. 0) then
    write(*,*) "!!! mmCrutial: dgetrs error.", error ,nn
    call abort(93, "Error in solve liner system (dgetrs)")
  endif

  call toc("Flexible wing (solve_linear_system_wing)", MPI_wtime() - t0)

end subroutine

subroutine Moving_boundary_point(wings)


    implicit none
    type(flexible_wing), dimension(1:nWings), intent(inout) :: wings
    integer :: i,j

    do i=1,nWings
        do j=1,nVeins_BC

          wings(i)%x_BC(-1,j) = wings(i)%x_BC(-1,j) + 1.0d-2
          wings(i)%y_BC(-1,j) = wings(i)%y_BC(-1,j) + 1.0d-2
          wings(i)%z_BC(-1,j) = wings(i)%z_BC(-1,j) + 1.0d-2

          wings(i)%x_BC(0,j) = wings(i)%x_BC(-1,j)
          wings(i)%y_BC(0,j) = wings(i)%y_BC(-1,j)
          wings(i)%z_BC(0,j) = wings(i)%z_BC(-1,j)

        enddo
      enddo
end subroutine

subroutine rotate_wing(wing)

  implicit none
  type(flexible_wing), intent(inout) :: wing
  real(kind=pr), dimension(1:3,1:3) :: mat_Rx, mat_Ry, mat_Rz
  integer :: i
  real(kind=pr), dimension(1:3) :: u


  !The wing is rotated based on conventional Euler angles. It is then rotated in
  !order: around z axis (yaw) first, then y axis (pitch) and finally x axis (roll)

    call Rx(mat_Rx,wing%WingAngle_x)
    call Ry(mat_Ry,wing%WingAngle_y)
    call Rz(mat_Rz,wing%WingAngle_z)

  ! Rotate wing around x axis
  do i=1,wing%np
    u = matmul(mat_Rx,(/wing%u_old(i), wing%u_old(i+wing%np), wing%u_old(i+2*wing%np)/))
    wing%x(i) = u(1)
    wing%y(i) = u(2)
    wing%z(i) = u(3)
  enddo

  ! Rotate wing around y axis
  do i=1,wing%np
    u = matmul(mat_Ry,(/wing%x(i), wing%y(i), wing%z(i)/))
    wing%x(i) = u(1)
    wing%y(i) = u(2)
    wing%z(i) = u(3)
  enddo

  ! Rotate wing around z axis
  do i=1,wing%np
    u = matmul(mat_Rz,(/wing%x(i), wing%y(i), wing%z(i)/))
    wing%x(i) = u(1)
    wing%y(i) = u(2)
    wing%z(i) = u(3)
  enddo

end subroutine

subroutine translate_wing(wing)

  implicit none
  type(flexible_wing), intent(inout) :: wing
  integer :: j, np
  real(kind=pr), dimension(1:3) :: u,v

  np = wing%np

  call rotate_vector_into_global_system(wing,wing%vr0)

  do j=1,np

    u = cross((/wing%x(j),wing%y(j),wing%z(j)/)&
              ,wing%vr0(1:3))


    v(1) = wing%u_old(3*wing%np+j)
    v(2) = wing%u_old(4*wing%np+j)
    v(3) = wing%u_old(5*wing%np+j)

    call rotate_vector_into_global_system(wing,v)

    wing%vx(j) = wing%vt0(1) + v(1) + u(1)
    wing%vy(j) = wing%vt0(2) + v(2) + u(2)
    wing%vz(j) = wing%vt0(3) + v(3) + u(3)

    wing%x(j) = wing%x(j) + wing%x0
    wing%y(j) = wing%y(j) + wing%y0
    wing%z(j) = wing%z(j) + wing%z0

  enddo

end subroutine

subroutine rotate_vector_into_wing_system(wing,Vector)

implicit none

!real(kind=pr),intent(in) :: time
type(flexible_wing), intent (inout) :: wing
real(kind=pr), intent(inout) :: Vector(1:3)
real(kind=pr), dimension(1:3,1:3) :: mat_Rx, mat_Ry, mat_Rz
real(kind=pr), dimension(1:3) :: u


!The wing is rotated based on conventional Euler angles. It is then rotated in
!order: around z axis (yaw) first, then y axis (pitch) and finally x axis (roll)

  call Rx(mat_Rx,-wing%WingAngle_x)
  call Ry(mat_Ry,-wing%WingAngle_y)
  call Rz(mat_Rz,-wing%WingAngle_z)

! Rotate wing around x axis
  u = matmul(mat_Rz,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)

! Rotate wing around x axis
  u = matmul(mat_Ry,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)

! Rotate wing around x axis
  u = matmul(mat_Rx,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)


end subroutine

subroutine rotate_vector_into_global_system(wing,Vector)

implicit none

!real(kind=pr),intent(in) :: time
type(flexible_wing), intent (inout) :: wing
real(kind=pr), intent(inout) :: Vector(1:3)
real(kind=pr), dimension(1:3,1:3) :: mat_Rx, mat_Ry, mat_Rz
real(kind=pr), dimension(1:3) :: u


!The wing is rotated based on conventional Euler angles. It is then rotated in
!order: around z axis (yaw) first, then y axis (pitch) and finally x axis (roll)

  call Rx(mat_Rx,wing%WingAngle_x)
  call Ry(mat_Ry,wing%WingAngle_y)
  call Rz(mat_Rz,wing%WingAngle_z)

! Rotate wing around x axis
  u = matmul(mat_Rx,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)

! Rotate wing around x axis
  u = matmul(mat_Ry,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)

! Rotate wing around x axis
  u = matmul(mat_Rz,(/Vector(1), Vector(2), Vector(3)/))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)


end subroutine


logical function Vector_isNAN(f)
  real(kind=pr), intent(in) :: f(:)
  integer :: a, i
  a = size(f)
  Vector_isNAN = .false.
  do i=1,a
    if (.not.(f(i).eq.f(i))) then
      Vector_isNAN = .true.
      return
    endif
  enddo
end function

logical function Matrix_isNAN(f,i,j)
  real(kind=pr), intent(in) :: f(:,:)
  integer, intent(out) :: i, j

  Matrix_isNAN = .false.
  do i=1,size(f,DIM=1)
    do j=1,size(f,DIM=2)
      if (.not.(f(i,j).eq.f(i,j))) then
        Matrix_isNAN = .true.
        if (root) then
          write(*,*) i, j
        endif
        !return
      endif
    enddo
  enddo
end function

subroutine truncated_triangular_prism_centroid(centroid,tri1,tri2,tri3,press_tri1,press_tri2,press_tri3)
!This function calculates the centroid of a truncated triangular prism by simple
!barycentric interpolation
! NOTICE: this can also be done by calculating directly the centroid of the truncated
! formed by the triangle and the pressure at the three vertices by deviding it
! into an ordinary triangular prism "minus" a pyramid

  implicit none
  real(kind=pr), dimension(1:3), intent(inout) :: centroid
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,tri3
  !array contains the pressure at the three vertices of the triangle tri1, tri2, tri3
  real(kind=pr), intent(in) :: press_tri1,press_tri2,press_tri3
  !the normal vector of the triangle
!  real(kind=pr), dimension(1:3), intent(in) :: normal

  centroid(1:3) =(press_tri1*tri1(1:3) + &
                                             press_tri2*tri2(1:3) + &
                                             press_tri3*tri3(1:3))/ &
                                            (press_tri1 + press_tri2 + press_tri3)

end subroutine

real(kind=pr) function pressure_trilinear_interp(x0, dx, field, x_target, periodic)
  use vars, only: pr, abort, mpirank, ga, gb
  use mpi
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
  real(kind=pr),intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
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

  if (ga(1)==0) then
    lbounds(1) = ga(1)
  else
    lbounds(1) = ga(1) - 1
  endif

  if (ga(2)==0) then
    lbounds(2) = ga(2)
  else
    lbounds(2) = ga(2) - 1
  endif

  if (ga(3)==0) then
    lbounds(3) = ga(3)
  else
    lbounds(3) = ga(3) - 1
  endif

  ubounds = gb(1:3)!(/size(field,1), size(field,2), size(field,3)/) - 1

  ! indices of cube containing the target point, lower end
  ix = floor( (x_target(1)-x0(1))/dx(1) ) + ga(1)
  iy = floor( (x_target(2)-x0(2))/dx(2) ) + ga(2)
  iz = floor( (x_target(3)-x0(3))/dx(3) ) + ga(3)

  ! distance to lower point, normalized (0..1)
  xd = ( x_target(1)-(dble(ix)*dx(1)+x0(1)) ) / dx(1)
  yd = ( x_target(2)-(dble(iy)*dx(2)+x0(2)) ) / dx(2)
  zd = ( x_target(3)-(dble(iz)*dx(3)+x0(3)) ) / dx(3)

  ! if the point is not on the grid, return a very large, negative value
  pressure_trilinear_interp = -9.9d10

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

          pressure_trilinear_interp = c0*(1.d0-zd)+c1*zd
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

          pressure_trilinear_interp = c0*(1.d0-zd)+c1*zd
        endif
      endif
    endif

    !if (pressure_trilinear_interp==-9.9d10) then
    !write(*,'("I am ",i5," CPU in interp")') mpirank
    !write(*,*) "bounds"
    !write(*,*) lbounds(1), ix, ubounds(1), lbounds(2), iy, ubounds(2), lbounds(3), iz, ubounds(3)
    !write(*,*) dx(3), x0(3), x_target(3)
    !endif
    !write(*,*) "Interpolated pressure before"
    !write(*,*) pressure_trilinear_interp
  endif

  ! with periodic it can work, but non-periodic not
  if (size(field,1)<=1) then
    call abort(9997111,"linear interpolation for 2d simulations currently not implemented...have fun")
  end if
end function pressure_trilinear_interp
