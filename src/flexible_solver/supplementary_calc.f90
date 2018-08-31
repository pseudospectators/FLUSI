subroutine length_calculation(x,y,z,springs,springs_length)

real(kind=pr), intent(in) :: x,y,z
real(kind=pr), intent(in) :: springs
real(kind=pr), allocatable, intent(out) :: springs_length(:)
integer :: ind

do ind=1:nint(maxval(springs(:,1)))
    springs_length(ind,1) = sqrt((x(nint(springs(ind,2)))-x(nint(springs(ind,3))))^2 + &
                                 (y(nint(springs(ind,2)))-y(nint(springs(ind,3))))^2 + &
                                 (z(nint(springs(ind,2)))-z(nint(springs(ind,3))))^2)
enddo

end subroutine

subroutine angle_calculation(x,y,z,vein,theta,phi)
!Calculate relative angles between two consecutive segments belonging to a vein
!   x(vein(ind,2))          x(vein(ind,3))                x(vein(ind,4))
!        O-----------------------O----------------------------O
!                  X1                          X2
!

real(kind=pr), intent(in) :: x,y,z
real(kind=pr), intent(in) :: vein
real(kind=pr), allocatable, intent(inout) :: theta(:), phi(:)
integer :: ind

do ind=1:nint(maxval(vein(:,1)))
    x1 = x(nint(vein(ind,3))) - x(nint(vein(ind,2)))
    y1 = y(nint(vein(ind,3))) - y(nint(vein(ind,2)))
    z1 = z(nint(vein(ind,3))) - z(nint(vein(ind,2)))
    x2 = x(nint(vein(ind,4))) - x(nint(vein(ind,3)))
    y2 = y(nint(vein(ind,4))) - y(nint(vein(ind,3)))
    z2 = z(nint(vein(ind,4))) - z(nint(vein(ind,3)))

    theta(ind) = atan2(x1*y2-x2*y1,x1*x2+y1*y2)

      phi(ind) = atan2(x1*z2-x2*z1,x1*x2+z1*z2)
enddo

end subroutine

subroutine convert_flexural_rigidity_into_spring_stiffness(EIy, EIz, kby0, kbz0, L, ns)

  implicit none
  real(kind=pr), intent(in) :: EIy, EIz, L, ns
  real(kind=pr), intent(out) :: kby0, kbz0

  kby0 = EIy*(ns-1)/L
  kbz0 = EIz*(ns-1)/L

end subroutine
