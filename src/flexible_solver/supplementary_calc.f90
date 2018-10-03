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

  time_LAPACK = time_LAPACK + MPI_wtime() - t0

  write(*,*) maxval(x)

end subroutine

subroutine Moving_boundary_point(wings)


    implicit none
    type(wing), dimension(1:nWings), intent(inout) :: wings
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

subroutine rotate_wing(wings)

  implicit none
  type(wing), intent(inout) :: wings
  real(kind=pr), dimension(1:3,1:3) :: mat_Ry, mat_Rz
  integer :: i
  real(kind=pr), dimension(1:3) :: u

    call Ry(mat_Ry,wings%Anglewing_y)
    call Rz(mat_Rz,wings%Anglewing_z)

  ! Rotate wing around y axis
  do i=1,wings%np
    u = matmul(mat_Ry,(/wings%x(i), wings%y(i), wings%z(i)/))
    wings%x(i) = u(1)
    wings%y(i) = u(2)
    wings%z(i) = u(3)
  enddo
  write(*,*) 'rotate'
  ! Rotate wing around z axis
  do i=1,wings%np
    u = matmul(mat_Rz,(/wings%x(i), wings%y(i), wings%z(i)/))
    wings%x(i) = u(1)
    wings%y(i) = u(2)
    wings%z(i) = u(3)
  enddo

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
