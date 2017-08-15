! datatype for wing kinematics, if described by a Fourier series or kineloader
! For both wings such a datatype is contained in the insect.
type wingkinematics
  ! Fourier coefficients
  real(kind=pr) :: a0_alpha, a0_phi, a0_theta=0.d0
  real(kind=pr), dimension(1:nfft_max) :: ai_phi=0.d0, bi_phi=0.d0, ai_theta=0.d0, &
    bi_theta=0.d0, ai_alpha=0.d0, bi_alpha=0.d0
  integer :: nfft_phi=0, nfft_alpha=0, nfft_theta=0
  ! coefficients are read only once from file (or set differently)
  logical :: initialized = .false.
  ! some details about the file, if reading from ini file
  character(len=strlen) :: infile_convention="", infile_type="", infile_units="", infile=""
  ! variables for kineloader (which uses non-periodic hermite interpolation)
  integer :: nk=0
  real(kind=pr), dimension (1:nhrmt_max) :: vec_t=0.d0, &
    vec_phi=0.d0,vec_alpha=0.d0,vec_theta=0.d0,vec_pitch=0.d0,vec_vert=0.d0,vec_horz=0.d0,  &
    vec_phi_dt=0.d0,vec_alpha_dt=0.d0,vec_theta_dt=0.d0,vec_pitch_dt=0.d0,vec_vert_dt=0.d0, &
    vec_horz_dt=0.d0
end type


!-----------------------------------------------------------------------------
! derived datatype for insect parameters (for readability)
type diptera
  !-------------------------------------------------------------
  ! Body motion state, wing motion state and characteristic points on insect
  !-------------------------------------------------------------
  ! position of logical center, and translational velocity
  real(kind=pr), dimension(1:3) :: xc_body_g=0.d0, vc_body_g=0.d0
  ! initial or tethered position, velocity and yawpitchroll angles:
  real(kind=pr), dimension(1:3) :: x0=0.d0, v0=0.d0, yawpitchroll_0=0.d0
  ! roll pitch yaw angles and their time derivatives
  real(kind=pr) :: psi=0.d0, beta=0.d0, gamma=0.d0, psi_dt=0.d0, beta_dt=0.d0, gamma_dt=0.d0, eta0=0.d0
  ! body pitch angle, if it is constant (used in forward flight and hovering)
  real(kind=pr) :: body_pitch_const=0.d0
  ! angles of the wings (left and right)
  real(kind=pr) :: phi_r=0.d0, alpha_r=0.d0, theta_r=0.d0, phi_dt_r=0.d0, alpha_dt_r=0.d0, theta_dt_r=0.d0
  real(kind=pr) :: phi_l=0.d0, alpha_l=0.d0, theta_l=0.d0, phi_dt_l=0.d0, alpha_dt_l=0.d0, theta_dt_l=0.d0
  ! stroke plane angle
  real(kind=pr) :: eta_stroke=0.d0
  ! angular velocity vectors (wings L+R, body)
  real(kind=pr), dimension(1:3) :: rot_body_b=0.d0, rot_body_g=0.d0
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_w=0.d0, rot_rel_wing_r_w=0.d0
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_b=0.d0, rot_rel_wing_r_b=0.d0
  real(kind=pr), dimension(1:3) :: rot_rel_wing_l_g=0.d0, rot_rel_wing_r_g=0.d0
  real(kind=pr), dimension(1:3) :: rot_abs_wing_l_g=0.d0, rot_abs_wing_r_g=0.d0
  ! angular acceleration vectors (wings L+R)
  real(kind=pr), dimension(1:3) :: rot_dt_wing_l_w=0.d0, rot_dt_wing_r_w=0.d0
  real(kind=pr), dimension(1:3) :: rot_dt_wing_l_g=0.d0, rot_dt_wing_r_g=0.d0
  ! Vector from body centre to pivot points in global reference frame
  real(kind=pr), dimension(1:3) :: x_pivot_l_g=0.d0, x_pivot_r_g=0.d0
  ! vectors desribing the positoions of insect's key elements
  ! in the body coordinate system
  real(kind=pr), dimension(1:3) :: x_head=0.d0,x_eye_r=0.d0,x_eye_l=0.d0,x_pivot_l=0.d0,x_pivot_r=0.d0
  ! moments of inertia in the body reference frame
  real(kind=pr) :: Jroll_body=0.d0, Jyaw_body=0.d0, Jpitch_body=0.d0
  ! total mass of insect:
  real(kind=pr) :: mass, gravity=0.d0
  !-------------------------------------------------------------
  ! for free flight solver
  !-------------------------------------------------------------
  real(kind=pr) :: time
  real(kind=pr), dimension(1:20) :: RHS_old=0.d0, RHS_this=0.d0
  real(kind=pr), dimension(1:20) :: STATE=0.d0
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
  real(kind=pr), dimension(1:6) :: DoF_on_off=0.d0
  character(len=strlen) :: startup_conditioner=""
  !-------------------------------------------------------------
  ! for wing fsi solver
  !-------------------------------------------------------------
  character(len=strlen) :: wing_fsi="no"
  real(kind=pr), dimension(1:3) :: torque_muscle_l_w=0.d0, torque_muscle_r_w=0.d0
  real(kind=pr), dimension(1:3) :: torque_muscle_l_b=0.d0, torque_muscle_r_b=0.d0
  real(kind=pr), dimension(1:3) :: init_alpha_phi_theta=0.d0
  !-------------------------------------------------------------
  ! wing shape parameters
  !-------------------------------------------------------------
  ! wing shape fourier coefficients. Note notation:
  ! R = a0/2 + SUM ( ai cos(2pi*i) + bi sin(2pi*i)  )
  ! to avoid compatibility issues, the array is of fixed size, although only
  ! the first nftt_wings entries will be used
  real(kind=pr), dimension(1:nfft_max) :: ai_wings=0.d0, bi_wings=0.d0
  real(kind=pr) :: a0_wings=0.d0
  ! fill the R0(theta) array once, then only table-lookup instead of Fseries
  real(kind=pr), dimension(1:25000) :: R0_table=0.d0
  ! describes the origin of the wings system
  real(kind=pr) :: xc=0.d0,yc=0.d0
  ! number of fft coefficients for wing geometry
  integer :: nfft_wings=0
  logical :: wingsetup_done = .false.
  logical :: wings_radius_table_ready = .false.
  ! wing bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
  real(kind=pr) :: wing_bounding_box(1:6) = 0.d0
  ! wing inertia
  real(kind=pr) :: Jxx=0.d0,Jyy=0.d0,Jzz=0.d0,Jxy=0.d0
  character(len=strlen) :: wing_thickness_distribution = "constant"
  ! a non-constant thickness is stored here:
  real(kind=pr), allocatable, dimension(:,:) :: wing_thickness_profile
  ! wing corrugation (i.e. the deviation from a flat wing)
  logical :: corrugated = .false.
  real(kind=pr), allocatable, dimension(:,:) :: corrugation_profile

  !--------------------------------------------------------------
  ! Wing kinematics
  !--------------------------------------------------------------
  ! wing kinematics Fourier coefficients
  type(wingkinematics) :: kine_wing_l, kine_wing_r

  !-------------------------------------------------------------
  ! parameters that control shape of wings,body, and motion
  !-------------------------------------------------------------
  character(len=strlen) :: WingShape="", BodyType="", BodyMotion="", HasDetails=""
  character(len=strlen) :: FlappingMotion_right="", FlappingMotion_left=""
  character(len=strlen) :: infile="", LeftWing="", RightWing=""
  ! parameters for body:
  real(kind=pr) :: L_body=0.d0, b_body=0.d0, R_head=0.d0, R_eye=0.d0
  ! parameters for wing shape:
  real(kind=pr) :: b_top=0.d0, b_bot=0.d0, L_chord=0.d0, L_span=0.d0, WingThickness=0.d0
  ! this is a safety distance for smoothing:
  real(kind=pr) :: safety=0.d0, smooth=0.d0
  ! parameter for hovering:
  real(kind=pr) :: distance_from_sponge=0.d0
  ! Wings and body forces (1:body,2:left wing,3:right wing)
  type(Integrals), dimension(1:3) :: PartIntegrals

end type diptera
!-----------------------------------------------------------------------------
