!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine prescribed_wing ( time, wing )
  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: i


  select case (wing%Motion)
  case ("prescribed_revolving_wing")
    call prescribed_revolving_wing (time, wing)
  case ("stationary")
    continue
  end select
end subroutine

subroutine prescribed_revolving_wing (time, wing)

  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: j
  real(kind=pr),dimension(1:3,1:3) :: mat_Rx, mat_Ry
  real(kind=pr),dimension(1:3) :: u,v
  real(kind=pr) :: tau,phi_z

  !call Rx(mat_Rx,wing%WingAngle_x)

  !do j=1,wing%np
  !  u = matmul(mat_Rx,(/wing%u_old(j) - wing%x0, &
  !                      wing%u_old(wing%np + j) -wing%y0, &
  !                      wing%u_old(2*wing%np + j) - wing%z0/))
  !  wing%x(j) = u(1) + wing%x0
  !  wing%y(j) = u(2) + wing%y0
  !  wing%z(j) = u(3) + wing%z0
  !enddo

  tau = 4.d-1
  wing%vr0 = (/0.d0,0.d0,0.d0/)

  phi_z= -1.0d0*(tau*dexp(-time/tau) + time) + 1.0d0*tau
  wing%vr0(3) = 1.0d0*(dexp(-time/tau) - 1)

  call Rz(mat_Ry,phi_z)

  do j=1,wing%np
    u = matmul(mat_Ry,(/wing%u_old(j) - wing%x0, &
                        wing%u_old(wing%np + j) -wing%y0, &
                        wing%u_old(2*wing%np + j) - wing%z0/))
    !u = matmul(mat_Ry,(/wing%x(j) - wing%x0, &
    !                    wing%y(j) -wing%y0, &
    !                    wing%z(j) - wing%z0/))
    wing%x(j) = u(1) + wing%x0
    wing%y(j) = u(2) + wing%y0
    wing%z(j) = u(3) + wing%z0

    v = -cross(wing%vr0(1:3),u(1:3))
    wing%vx(j) = v(1)
    wing%vy(j) = v(2)
    wing%vz(j) = v(3)
enddo

end subroutine
