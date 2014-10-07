!-----------------------
! Rotation Matrices
!-----------------------
subroutine Rx (R,angle)
  use vars
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R
  R(1,:) = (/ 1.d0, 0.d0, 0.d0/)
  R(2,:) = (/ 0.d0, cos(angle), sin(angle) /)
  R(3,:) = (/ 0.d0, -sin(angle), cos(angle) /)
end subroutine 


subroutine Ry (R,angle)
  use vars
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R
  R(1,:) = (/ cos(angle), 0.d0, -sin(angle)/)
  R(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
  R(3,:) = (/ +sin(angle), 0.d0, cos(angle) /)
end subroutine 


subroutine Rz (R,angle)
  use vars
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R
  R(1,:) = (/ cos(angle), +sin(angle), 0.d0/)
  R(2,:) = (/ -sin(angle), cos(angle), 0.d0/)
  R(3,:) = (/ 0.d0, 0.d0, 1.d0  /)  
end subroutine 
