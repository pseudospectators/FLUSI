subroutine integrate_position (time, beam) !attention! parameter time added (for mouvement)
! this function integrates the incoming beam(:,5) and beam(:,6) (angle and angle_dot) and returns the new positions
! in the same array
  implicit none
  integer :: i,s
  real(kind=pr) :: xrel, yrel, vxrel, vyrel
  type(solid), intent(inout) :: beam
  real(kind=pr) :: alpha, alpha_t, alpha_tt
  real(kind=pr), intent(in) :: time
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)

  call mouvement ( time, alpha, alpha_t, alpha_tt, LeadingEdge, beam )

  ! delete old positions + velocities
  beam%x = 0.d0
  beam%y = 0.d0
  beam%vx = 0.d0
  beam%vy = 0.d0

  do s=0,ns-1
    xrel = 0.d0
    yrel = 0.d0
    vxrel = 0.d0
    vyrel = 0.d0

    do i=0,s-1
      xrel = xrel + dcos( beam%theta(i) )*ds
      yrel = yrel + dsin( beam%theta(i) )*ds
      vxrel = vxrel - beam%theta_dot(i) * dsin( beam%theta(i) )*ds
      vyrel = vyrel + beam%theta_dot(i) * dcos( beam%theta(i) )*ds
    enddo
    ! these sums now contain relative pos / vel
    ! convert them to absolute values:
    beam%x(s)  = LeadingEdge(1) + xrel*dcos(alpha) - yrel*dsin(alpha)
    beam%y(s)  = LeadingEdge(2) + yrel*dcos(alpha) + xrel*dsin(alpha)
    beam%vx(s) = LeadingEdge(3) + dcos(alpha)*(vxrel - alpha_t*yrel) - dsin(alpha)*(vyrel + alpha_t*xrel)
    beam%vy(s) = LeadingEdge(4) + dcos(alpha)*(vyrel + alpha_t*xrel) + dsin(alpha)*(vxrel - alpha_t*yrel)
  enddo


end subroutine integrate_position
