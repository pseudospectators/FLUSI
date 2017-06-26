! datatype for wing kinematics, if described by a Fourier series or kineloader
! For both wings such a datatype is contained in the insect.
type wingkinematics
  ! Fourier coefficients
  real(kind=pr) :: a0_alpha, a0_phi, a0_theta
  real(kind=pr), dimension(1:nfft_max) :: ai_phi, bi_phi, ai_theta, bi_theta, ai_alpha, bi_alpha
  integer :: nfft_phi, nfft_alpha, nfft_theta
  ! coefficients are read only once from file (or set differently)
  logical :: initialized = .false.
  ! some details about the file, if reading from ini file
  character(len=strlen) :: infile_convention="", infile_type="", infile_units="", infile=""
  ! variables for kineloader (which uses non-periodic hermite interpolation)
  integer :: nk
  real(kind=pr), dimension (1:nhrmt_max) :: vec_t, &
    vec_phi,vec_alpha,vec_theta,vec_pitch,vec_vert,vec_horz,  &
    vec_phi_dt,vec_alpha_dt,vec_theta_dt,vec_pitch_dt,vec_vert_dt,vec_horz_dt
end type


!-----------------------------------------------------------------------------
! derived datatype for insect parameters (for readability)
type diptera
  !-------------------------------------------------------------
  ! Body motion state, wing motion state and characteristic points on insect
  !-------------------------------------------------------------
  ! position of logical center, and translational velocity
  real(kind=pr), dimension(1:3) :: xc_body_g, vc_body_g
  ! initial or tethered position, velocity and yawpitchroll angles:
  real(kind=pr), dimension(1:3) :: x0, v0, yawpitchroll_0
  ! roll pitch yaw angles and their time derivatives
  real(kind=pr) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt, eta0
  ! body pitch angle, if it is constant (used in forward flight and hovering)
  real(kind=pr) :: body_pitch_const
  ! angles of the wings (left and right)
  real(kind=pr) :: phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r
  real(kind=pr) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l
  ! stroke plane angle
  real(kind=pr) :: eta_stroke
  ! angular velocity vectors (wings L+R, body)
  real(kind=pr), dimension(1:3) :: rot_body_b, rot_body_g
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_w, rot_rel_wing_r_w
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_b, rot_rel_wing_r_b
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_g, rot_rel_wing_r_g
  real(kind=pr), dimension(1:3) :: rot_abs_wing_l_g, rot_abs_wing_r_g
  ! angular acceleration vectors (wings L+R)
  real(kind=pr), dimension(1:3) :: rot_dt_wing_l_w, rot_dt_wing_r_w
  real(kind=pr), dimension(1:3) :: rot_dt_wing_l_g, rot_dt_wing_r_g
  ! Vector from body centre to pivot points in global reference frame
  real(kind=pr), dimension(1:3) :: x_pivot_l_g, x_pivot_r_g
  ! vectors desribing the positoions of insect's key elements
  ! in the body coordinate system
  real(kind=pr), dimension(1:3) :: x_head,x_eye_r,x_eye_l,x_pivot_l,x_pivot_r
  ! moments of inertia in the body reference frame
  real(kind=pr) :: Jroll_body, Jyaw_body, Jpitch_body
  ! total mass of insect:
  real(kind=pr) :: mass, gravity
  !-------------------------------------------------------------
  ! for free flight solver
  !-------------------------------------------------------------
  real(kind=pr) :: time
  real(kind=pr), dimension(1:20) :: RHS_old, RHS_this
  real(kind=pr), dimension(1:20) :: STATE
  ! STATE(1) : x-position of body
  ! STATE(2) : y-position of body
  ! STATE(3) : z-position of body
  ! STATE(4) : x-velocity of body
  ! STATE(5) : y-velocity of body
  ! STATE(6) : z-velocity of body
  ! STATE(7) : 1st component of body quaternion
  ! STATE(8) : 2nd component of body quaternion
  ! STATE(9) : 3rd component of body quaternion
  ! STATE(10) : 4th component of body quaternion
  ! STATE(11) : x-angular velocity of body (in body system)
  ! STATE(12) : y-angular velocity of body (in body system)
  ! STATE(13) : z-angular velocity of body (in body system)
  ! STATE(14) : 1st component of left wing quaternion
  ! STATE(15) : 2nd component of left wing quaternion
  ! STATE(16) : 3rd component of left wing quaternion
  ! STATE(17) : 4th component of left wing quaternion
  ! STATE(18) : x-angular velocity of left wing
  ! STATE(19) : y-angular velocity of left wing
  ! STATE(20) : z-angular velocity of left wing
  real(kind=pr), dimension(1:6) :: DoF_on_off
  character(len=strlen) :: startup_conditioner
  !-------------------------------------------------------------
  ! for wing fsi solver
  !-------------------------------------------------------------
  character(len=strlen) :: wing_fsi
  real(kind=pr), dimension(1:3) :: torque_muscle_l_w, torque_muscle_r_w
  real(kind=pr), dimension(1:3) :: torque_muscle_l_b, torque_muscle_r_b
  real(kind=pr), dimension(1:3) :: init_alpha_phi_theta
  !-------------------------------------------------------------
  ! wing shape parameters
  !-------------------------------------------------------------
  ! wing shape fourier coefficients. Note notation:
  ! R = a0/2 + SUM ( ai cos(2pi*i) + bi sin(2pi*i)  )
  ! to avoid compatibility issues, the array is of fixed size, although only
  ! the first nftt_wings entries will be used
  real(kind=pr), dimension(1:nfft_max) :: ai_wings, bi_wings
  real(kind=pr) :: a0_wings
  ! fill the R0(theta) array once, then only table-lookup instead of Fseries
  real(kind=pr), dimension(1:25000) :: R0_table
  ! describes the origin of the wings system
  real(kind=pr) :: xc,yc
  ! number of fft coefficients for wing geometry
  integer :: nfft_wings
  logical :: wingsetup_done = .false.
  logical :: wings_radius_table_ready = .false.
  ! wing bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
  real(kind=pr) :: wing_bounding_box(1:6) = 0.d0
  ! wing inertia
  real(kind=pr) :: Jxx,Jyy,Jzz,Jxy
  character(len=strlen) :: wing_thickness_distribution = "constant"
  ! a non-constant thickness is stored here:
  real(kind=pr), allocatable, dimension(:,:) :: wing_thickness_profile
  ! wing corrugation (i.e. the deviation from a flat wing)
  logical :: corrugated
  real(kind=pr), allocatable, dimension(:,:) :: corrugation_profile

  !--------------------------------------------------------------
  ! Wing kinematics
  !--------------------------------------------------------------
  ! wing kinematics Fourier coefficients
  type(wingkinematics) :: kine_wing_l, kine_wing_r

  !-------------------------------------------------------------
  ! parameters that control shape of wings,body, and motion
  !-------------------------------------------------------------
  character(len=strlen) :: WingShape, BodyType, BodyMotion, HasDetails
  character(len=strlen) :: FlappingMotion_right, FlappingMotion_left
  character(len=strlen) :: infile, LeftWing, RightWing
  ! parameters for body:
  real(kind=pr) :: L_body, b_body, R_head, R_eye
  ! parameters for wing shape:
  real(kind=pr) :: b_top, b_bot, L_chord, L_span, WingThickness
  ! this is a safety distance for smoothing:
  real(kind=pr) :: safety, smooth
  ! parameter for hovering:
  real(kind=pr) :: distance_from_sponge
  ! Wings and body forces (1:body,2:left wing,3:right wing)
  type(Integrals), dimension(1:3) :: PartIntegrals

end type diptera
!-----------------------------------------------------------------------------
