!-------------------------------------------------------------------------------
! Body motion protocoll, different choices.
! Input:
!      time (self explanatory)
! Output:
!      Insect% psi:     roll angle
!      Insect%beta:     pitch angle
!      Insect%gamma:    yaw angle
!      Insect%psi_dt:   roll angular velocity
!      Insect%beta_dt:  pitch angular velocity
!      Insect%gamma_dt: yaw angular velocity
!      Insect%xc:       center of gravity coordinate
!      Insect%vc:       translational velocity of the body
! The actual motion depends on the choices in the parameter file, namely
! Insect%BodyMotion, and sub-parameters that may further precise a given motion
! protocoll
! Note that in new versions, all the angles and positions are stored in one
! datastructure, which is then the only output variable of this routine.
!-------------------------------------------------------------------------------
subroutine BodyMotion(time, Insect)
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt
  real(kind=pr) :: xc(1:3), vc(1:3), ep(0:3)
  real(kind=pr) :: T,R
  character(len=strlen) :: dummy

  ! the tag body_moves is used to draw the insect's body only once, if the body
  ! does not move (body_moves=="no"). For safety, we initialize the body as moving
  ! so if you forget to specify (body_moves=="no"), the body is drawn every time
  body_moves = "yes"


  select case (Insect%BodyMotion)
  case ("command-line")
    ! use values specified in global variables. this is used to draw a single
    ! mask with a dry run, where all 12 parameters are specified in the command line
    if(root) write(*,*) "Reading left wing kinematics (params 4,5,6:x,y,z and 7,8,9:psi,beta,gamma)"
    if(root) write(*,*) "note body does NOT move (no velocity field)"
    call get_command_argument(4,dummy)
    read (dummy,*) xc(1)
    call get_command_argument(5,dummy)
    read (dummy,*) xc(2)
    call get_command_argument(6,dummy)
    read (dummy,*) xc(3)
    call get_command_argument(7,dummy)
    read (dummy,*) psi
    call get_command_argument(8,dummy)
    read (dummy,*) beta
    call get_command_argument(9,dummy)
    read (dummy,*) gamma

    body_moves = "no"
    psi_dt   = 0.d0
    beta_dt  = 0.d0
    gamma_dt = 0.d0
    vc = (/0.d0, 0.d0, 0.d0/) ! tethered: no velocity

    if(root) write(*,'("x=",g12.4,"y=",g12.4,"z=",g12.4)') xc
    if(root) write(*,'("psi=",g12.4,"beta=",g12.4,"gamma=",g12.4)') psi,beta,gamma

  case ("roll")
    psi      = 30.d0*sin(2.d0*pi*time)
    beta     = 0.d0 ! pitch
    gamma    = 0.d0 ! yaw
    psi_dt   = 30.d0*cos(2.d0*pi*time)*2.d0*pi
    beta_dt  = 0.d0
    gamma_dt = 0.d0
    xc = Insect%x0
    vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
    body_moves = "yes"

  case ("pitch")
      psi      = 0.d0
      beta     = 30.d0*sin(2.d0*pi*time)
      gamma    = 0.d0 ! yaw
      psi_dt   = 0.d0
      beta_dt  = 30.d0*cos(2.d0*pi*time)*2.d0*pi
      gamma_dt = 0.d0
      xc = Insect%x0
      vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
      body_moves = "yes"

  case ("yaw")
    psi      = 0.d0
    beta     = 0.d0
    gamma    = 30.d0*sin(2.d0*pi*time)
    psi_dt   = 0.d0
    beta_dt  = 0.d0
    gamma_dt = 30.d0*cos(2.d0*pi*time)*2.d0*pi
    xc = Insect%x0
    vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
    body_moves = "yes"

  case ("tethered")
    psi      = Insect%yawpitchroll_0(3) ! roll
    beta     = Insect%yawpitchroll_0(2) ! pitch
    gamma    = Insect%yawpitchroll_0(1) ! yaw
    psi_dt   = 0.d0  ! tethered: angles const
    beta_dt  = 0.d0
    gamma_dt = 0.d0
    xc = Insect%x0
    vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
    body_moves = "no" ! tethered: body does not move

    if (Insect%BodyType=="suzuki_thin_rod") then
      body_moves = "yes"
    endif

  case ("free_flight")
    ! The case "free_flight" is different from the others. The free flight solver
    ! computes the current state of the insect in INSECT%STATE, which is a 13
    ! component vector (6 translation, 4 quaternions, 3 angular velocity)
    ! in this case, the position is dynamically computed, and quaternions are used
    body_moves = "yes"

    ! copy data from insect state vector
    xc              = Insect%STATE(1:3)
    vc              = Insect%STATE(4:6)
    ep              = Insect%STATE(7:10)
    Insect%rot_body = Insect%STATE(11:13)

    ! compute yaw pitch roll from the quaternion. attention: this is just for
    ! information to dump in the log file (may be useful), but NOT for the rotation
    ! matrix. otherwise, we could just omit the quaternion
    psi      = atan2(2*(ep(2)*ep(3) + ep(0)*ep(1)), ep(0)*ep(0) - ep(1)*ep(1) - ep(2)*ep(2) + ep(3)*ep(3))
    beta     = asin(-2*(ep(1)*ep(3) - ep(0)*ep(2)))
    gamma    = atan2(2*(ep(1)*ep(2) + ep(0)*ep(3)), ep(0)*ep(0) + ep(1)*ep(1) - ep(2)*ep(2) - ep(3)*ep(3))
    ! these values cannot easily be computed, but they are not really necessary
    psi_dt   = 0.0d0
    beta_dt  = 0.0d0
    gamma_dt = 0.0d0

    ! note yaw,pitch and roll do NOT enter the body rotation matrix, it is computed
    ! from the orientation quaternion
    call rotation_matrix_from_quaternion( ep , Insect%M_body_quaternion )

    if(root) write(*,'(f12.4,2x,9(es12.4,1x))') time, &
    xc,vc,Insect%rot_body

  case default
    if (mpirank==0) then
      write(*,*) Insect%BodyMotion
      write(*,*) "body_motion.f90::BodyMotion: motion case (Insect%BodyMotion) undefined"
      call abort()
    endif
  end select




  if ((mpirank==0).and.(maxval(vc)>0.0d0).and.(body_moves=="no")) then
    write(*,*) "error in body_motion.f90: I found maxval(vc)>0 but the body_moves"
    write(*,*) "flag is set to no, which means we will draw the body only once"
    write(*,*) "This is probably not intented - you should look into it."
    call abort()
  endif


  ! save above values in the insect
  Insect%psi      = psi
  Insect%beta     = beta
  Insect%gamma    = gamma
  Insect%psi_dt   = psi_dt
  Insect%beta_dt  = beta_dt
  Insect%gamma_dt = gamma_dt
  Insect%xc_body  = xc
  Insect%vc_body  = vc


  ! for compability, we update the x0,y0,z0 also
  ! this is used e.g. for torque computation
  x0 = xc(1)
  y0 = xc(2)
  z0 = xc(3)

end subroutine BodyMotion
