module share_kine
  implicit none

  integer, parameter :: prk = 8
  integer, save :: nk
  real (kind=prk), dimension (:), allocatable, save :: vec_t, &
    vec_phi,vec_alpha,vec_theta,vec_pitch,vec_vert,vec_horz,  &
    vec_phi_dt,vec_alpha_dt,vec_theta_dt,vec_pitch_dt,vec_vert_dt,vec_horz_dt
end module share_kine

