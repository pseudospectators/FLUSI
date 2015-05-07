!-------------------------------------------------------------------------------
! Rigid solid time stepping routines
! AB2 method with Euler startup
!-------------------------------------------------------------------------------
subroutine rigid_solid_time_step(time,dt0,dt1,it,Insect)
  use vars
  implicit none

  real (kind=pr),intent (in) :: time,dt1,dt0
  type(diptera),intent(inout)::Insect
  integer,intent (in) :: it
  real (kind=pr) :: b10,b11

  ! select scheme
  if(it == 0) then
    ! EULER startup scheme
    ! compute rhs at this time step (updates Insect%RHS_this)
    call rigid_solid_rhs(time,it,Insect)

    ! Euler step
    Insect%STATE = Insect%STATE + dt1*Insect%RHS_this
  else
    ! ADAMS-BASHFORTH 2 scheme
    ! update vectors
    Insect%RHS_OLD = Insect%RHS_this

    ! compute rhs at this time step (updates Insect%RHS_this)
    call rigid_solid_rhs(time,it,Insect)

    ! adaptive AB2 coefficients
    b10 = dt1/dt0*(0.5*dt1 + dt0)
    b11 = -0.5*dt1*dt1/dt0

    ! adaptive AB2 step
    Insect%STATE = Insect%STATE + b10*Insect%RHS_this + b11*Insect%RHS_old
  endif
  Insect%time = Insect%time + dt1
end subroutine rigid_solid_time_step



!-------------------------------------------------------------------------------
! Insect free flight dynamics.
! RHS of the ODE system.
! TASK: from INSECT%STATE_THIS compute INSECT%RHS_THIS
!-------------------------------------------------------------------------------
subroutine rigid_solid_rhs(time,it,Insect)
  use mpi
  use vars
  implicit none

  integer, intent(in) :: it
  real (kind=pr), intent (in) :: time
  type(diptera),intent(inout)::Insect
  real (kind=pr) :: m,g, Jx, Jy,Jz
  real (kind=pr), dimension(0:3) :: ep
  real (kind=pr), dimension(1:3) :: ROT, torque_body
  real (kind=pr), dimension(1:3,1:3) :: C_mat

  ! initialization
  Insect%RHS_this=0.d0

  if ((root).and.(Insect%BodyMotion/="free_flight")) then
    write(*,*) "We have a problem in Insect%BodyMotion"
    call abort()
  endif


  g  = Insect%gravity
  m  = Insect%mass
  Jx = Insect%Jroll_body
  Jy = Insect%Jpitch_body
  Jz = Insect%Jyaw_body

  ! extract rotation (unit) quaternion from state vector, and create the rotation
  ! matrix from it.
  ep = Insect%STATE(7:10)
  call rotation_matrix_from_quaternion( ep,Insect%M_body_quaternion )
  ! extract angular velocity vector from state vector
  ROT = Insect%STATE(11:13)

  ! The equations of motion for the rotation are written in the body reference
  ! frame, thus the fluid torque has to be transformed to the body system
  torque_body = matmul( Insect%M_body_quaternion, GlobalIntegrals%torque)

  !-----------------------------------------------------------------------------
  ! Integrate the 13 equations of motion (written in quaternion formulation)
  !-----------------------------------------------------------------------------
  ! integrates coordinates (dx/dt = vx)
  Insect%RHS_this(1) = Insect%STATE(4)
  Insect%RHS_this(2) = Insect%STATE(5)
  Insect%RHS_this(3) = Insect%STATE(6)
  ! integrates velocities (dvx/dt = F)
  Insect%RHS_this(4) = GlobalIntegrals%force(1)/m
  Insect%RHS_this(5) = GlobalIntegrals%force(2)/m
  Insect%RHS_this(6) = GlobalIntegrals%force(3)/m + g
  ! integrate quaternion attitudes
  Insect%RHS_this(7) = -ep(1)*ROT(1)-ep(2)*ROT(2)-ep(3)*ROT(3)
  Insect%RHS_this(8) = +ep(0)*ROT(1)-ep(3)*ROT(2)+ep(2)*ROT(3)
  Insect%RHS_this(9) = +ep(3)*ROT(1)+ep(0)*ROT(2)-ep(1)*ROT(3)
  Insect%RHS_this(10) =-ep(2)*ROT(1)+ep(1)*ROT(2)+ep(0)*ROT(3)

  Insect%RHS_this(11) = ( (Jy-Jz)*ROT(2)*ROT(3)+torque_body(1) )/Jx
  Insect%RHS_this(12) = ( (Jz-Jx)*ROT(3)*ROT(1)+torque_body(2) )/Jy
  Insect%RHS_this(13) = ( (Jx-Jy)*ROT(1)*ROT(2)+torque_body(3) )/Jz

  ! there is a factor 1/2 missing!
  Insect%RHS_this(7:13) = Insect%RHS_this(7:13) / 2.d0


end subroutine rigid_solid_rhs




!-------------------------------------------------------------------------------
! Initialize insect free flight dynamics solver.
! Task: set up the current state of the insect in INSECT%STATE
!-------------------------------------------------------------------------------
subroutine rigid_solid_init(time, Insect)
  use mpi
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real (kind=pr) :: yaw,pitch,roll
  integer :: mpicode
  real (kind=pr), dimension(0:3) :: ep
  real (kind=pr), dimension(1:3,1:3) :: C_mat

  Insect%time = time

  ! free flight solver based on quaternions. the task here is to initialize
  ! the "attitude" quaternion (Insect%quaternion) from yaw, pitch and roll
  ! angles. Note that for the free flight solver, this is the last time
  ! body yaw,pitch,roll are meaningful. from now on, quaternions are used

  ! initialization, these values are read from parameter file
  Insect%gamma = Insect%yawpitchroll_0(1)
  Insect%beta  = Insect%yawpitchroll_0(2)
  Insect%psi   = Insect%yawpitchroll_0(3)
  Insect%xc_body = Insect%x0
  Insect%vc_body = Insect%v0
  Insect%rot_body = 0.d0

  ! create initial value for attitude quaternion
  yaw   =  Insect%gamma / 2.d0
  pitch =  Insect%beta  / 2.d0
  roll  =  Insect%psi   / 2.d0
  ep(0) = cos(roll)*cos(pitch)*cos(yaw) + sin(roll)*sin(pitch)*sin(yaw)
  ep(1) = sin(roll)*cos(pitch)*cos(yaw) - cos(roll)*sin(pitch)*sin(yaw)
  ep(2) = cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*cos(pitch)*sin(yaw)
  ep(3) = cos(roll)*cos(pitch)*sin(yaw) - sin(roll)*sin(pitch)*cos(yaw)

  ! to draw the very first insect, we also need the rotation matix
  call rotation_matrix_from_quaternion( ep,Insect%M_body_quaternion )


  ! Assemble the state vector
  Insect%STATE = (/ Insect%xc_body(1),Insect%xc_body(2),Insect%xc_body(3), &
                    Insect%vc_body(1),Insect%vc_body(2),Insect%vc_body(3), &
                   ep(0),ep(1),ep(2),ep(3),&
                   Insect%rot_body(1),Insect%rot_body(2),Insect%rot_body(3) /)

  Insect%RHS_this = 0.0
  Insect%RHS_old = 0.0

  write (*,*) "Insect%STATE", Insect%STATE


end subroutine rigid_solid_init
