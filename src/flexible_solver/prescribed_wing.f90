!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine prescribed_wing ( time, wing , Insect)
  implicit none

  real(kind=pr),intent(in) :: time
  type(diptera), intent(inout) :: Insect
  type(flexible_wing), intent (inout) :: wing
  integer :: i


  select case (wing%Motion)
  case ("revolving_wing","from_file")
    call prescribed_wing_motion (time, wing, Insect)
  case ("stationary")
    continue
  end select
end subroutine

subroutine prescribed_wing_motion (time, wing, Insect)

  implicit none

  real(kind=pr),intent(in) :: time
  type(diptera), intent(inout) :: Insect
  type(flexible_wing), intent (inout) :: wing
  integer :: j,i
  real(kind=pr),dimension(1:3,1:3) :: mat_Rx, mat_Ry
  real(kind=pr),dimension(1:3) :: u,v
  !real(kind=pr) :: ttau, phi, phi_dt, alpha, alpha_dt, theta, theta_dt

  !-----------------------------------------------------------------------------
  ! fetch current motion state
  !-----------------------------------------------------------------------------
  call BodyMotion (time, Insect)
  call StrokePlane (time, Insect)
  if (wing%ID == "left") then
    call Flexible_wing_motions ( time, wing, Insect%kine_wing_l )
  elseif (wing%ID == "right") then
    call Flexible_wing_motions ( time, wing, Insect%kine_wing_r )
  endif

  !-----------------------------------------------------------------------------
  ! define the rotation matrices to change between coordinate systems
  !-----------------------------------------------------------------------------
  call body_rotation_matrix( Insect, Insect%M_body )
  Insect%M_body_inv = transpose(Insect%M_body)
  call MSM_solver_rotation_matrix( Wing, wing%M_solver )
  Wing%M_solver_inv = transpose(Wing%M_solver)
  call flexible_wing_rotation_matrix( Wing, Insect, Wing%M_wing )
  Wing%M_wing_inv = transpose(Wing%M_wing)

  ! rel+abs wing angular velocities in the w/b/g coordinate system
  call flexible_wing_angular_velocities (time, Wing, Insect, Insect%M_body )

  call rotate_and_translate_wing_into_global_system(wing, Insect)

  call construct_total_velocity(wing,Insect%M_body,Insect%M_body_inv)

end subroutine
