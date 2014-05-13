module kine
  implicit none

  integer, parameter :: prk = 8
  integer, save :: nk
  real (kind=prk), dimension (:), allocatable, save :: vec_t, &
    vec_phi,vec_alpha,vec_theta,vec_pitch,vec_vert,vec_horz,  &
    vec_phi_dt,vec_alpha_dt,vec_theta_dt,vec_pitch_dt,vec_vert_dt,vec_horz_dt
  ! Kinematic variables for the dynamic solver
  ! In this version (4 Dec 2013), only x and z coordinates are solved
  ! RHS at 2 time layers are required for AB
  ! kine_now contains RHS at current time step
  ! kine_old contains RHS at previous time step
  ! The system of ODEs:
  ! dx/dt = v_x
  ! dz/dt = v_z
  ! dv_x/dt = (f_x_external+f_x_aero_steadypart)/(m_solid-m_fluid)
  ! dv_z/dt = (f_z_external+f_z_aero_steadypart)/(m_solid-m_fluid) 
  ! The unsteady correction is thus treated implicitly
  real (kind=prk), dimension (4), save :: kine_now,kine_old
  ! All peremeters are normalized as everywhere else in the code
end module kine

