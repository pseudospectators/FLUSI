subroutine integrate_position (time, beam) !attention! parameter time added (for mouvement)
! this function integrates the incoming beam(:,5) and beam(:,6) (angle and angle_dot) and returns the new positions
! in the same array
  implicit none
  integer :: i,s
  real (kind=pr) :: xrel, yrel, vxrel, vyrel
  type(solid), intent(inout) :: beam
  real (kind=pr) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  
  call mouvement ( time, alpha, alpha_t, alpha_tt, LeadingEdge, beam )
  
  xrel = 0.d0
  yrel = 0.d0
  vxrel = 0.d0
  vyrel = 0.d0
  
  ! delete old positions + velocities
  beam%x = 0.d0
  beam%y = 0.d0
  beam%vx = 0.d0
  beam%vy = 0.d0
  
  do s=0,ns-1
    do i=0,s-1
      xrel = xrel + cos( beam%theta(i) )*ds
      yrel = yrel + sin( beam%theta(i) )*ds
      vxrel = vxrel - beam%theta_dot(i) * sin( beam%theta(i) )*ds 
      vyrel = vyrel + beam%theta_dot(i) * cos( beam%theta(i) )*ds 
    enddo 
    ! these sums now contain relative pos / vel 
    ! convert them to absolute values:
    beam%x(s)  = LeadingEdge(1) + xrel*cos(alpha) - yrel*sin(alpha)
    beam%y(s)  = LeadingEdge(2) + yrel*cos(alpha) + xrel*sin(alpha)
    beam%vx(s) = LeadingEdge(3) + cos(alpha)*(vxrel - alpha_t*yrel) - sin(alpha)*(vyrel + alpha_t*xrel)
    beam%vy(s) = LeadingEdge(4) + cos(alpha)*(vyrel + alpha_t*xrel) + sin(alpha)*(vxrel - alpha_t*yrel)
    xrel  = 0.d0
    yrel  = 0.d0
    vxrel = 0.d0
    vyrel = 0.d0    
  enddo
  

end subroutine integrate_position
