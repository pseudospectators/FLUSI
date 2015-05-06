! Rigid solid time stepping routines
! AB2 method with Euler startup
subroutine rigid_solid_time_step(time,dt0,dt1,it,Insect)
  use vars
  implicit none

  real (kind=pr),intent (in) :: time,dt1,dt0
  type(diptera),intent(inout)::Insect
  integer,intent (in) :: it
  real (kind=pr) :: b10,b11

  ! select scheme
  if(it == 0) then
    ! compute rhs at this time step
    call rigid_solid_rhs(time,it,Insect)

    ! Euler step
    Insect%STATE = Insect%STATE + dt1*Insect%RHS_this
  else
    ! update vectors
    Insect%RHS_OLD = Insect%RHS_this

    ! compute rhs at this time step (updates Insect%RHS_this)
    call rigid_solid_rhs(time,it,Insect)

    ! adaptive AB2 coefficients
    b10=dt1/dt0*(0.5*dt1 + dt0)
    b11=-0.5*dt1*dt1/dt0

    ! adaptive AB2 step
    Insect%STATE = Insect%STATE + b10*Insect%RHS_this + b11*Insect%RHS_old
  endif

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
  real (kind=pr), dimension(1:3) :: ROT
  real (kind=pr), dimension(1:3,1:3) :: C_mat

  ! initialization
  Insect%RHS_this=0.d0

  if ((root).and.(Insect%BodyMotion/="free_flight")) then
    write(*,*) "We have a problem in Insect%BodyMotion"
    call abort()
  endif

  Insect%mass_solid = 2.0
  Insect%gravity = 1.0

  ! gravity acceleration
  g = Insect%gravity
  m = Insect%mass_solid
  Jx = 10.d0
  Jy = 10.d0
  Jz = 0.1

  ep = Insect%STATE(7:10)
  ROT = Insect%STATE(11:13)


  call rotation_matrix_from_quaternion( ep,Insect%M_body_quaternion )

  ! integrates coordinates (dx/dt = vx)
  Insect%RHS_this(1) = Insect%STATE(4)
  Insect%RHS_this(2) = Insect%STATE(5)
  Insect%RHS_this(3) = Insect%STATE(6)
  ! integrates velocities (dvx/dt = F)
  Insect%RHS_this(4) = 0.0!GlobalIntegrals%force(1)/m
  Insect%RHS_this(5) = 0.0!GlobalIntegrals%force(2)/m
  Insect%RHS_this(6) = GlobalIntegrals%force(3)/m - g ! gravity always in neg. z
  ! integrate quaternion attitudes
  Insect%RHS_this(7) = -ep(1)*ROT(1)-ep(2)*ROT(2)-ep(3)*ROT(3)
  Insect%RHS_this(8) = +ep(0)*ROT(1)-ep(3)*ROT(2)+ep(2)*ROT(3)
  Insect%RHS_this(9) = +ep(3)*ROT(1)+ep(0)*ROT(2)-ep(1)*ROT(3)
  Insect%RHS_this(10) =-ep(2)*ROT(1)+ep(1)*ROT(2)+ep(0)*ROT(3)

  Insect%RHS_this(11) = 0.0!( (Jy-Jz)*ROT(2)*ROT(3)+GlobalIntegrals%torque(1) )/Jx
  Insect%RHS_this(12) = 0.0!( (Jz-Jx)*ROT(3)*ROT(1)+GlobalIntegrals%torque(2) )/Jy
  Insect%RHS_this(13) = 0.0!( (Jx-Jy)*ROT(1)*ROT(2)+GlobalIntegrals%torque(3) )/Jz

  ! there is a factor 1/2 missing!
  Insect%RHS_this(7:13) = Insect%RHS_this(7:13) / 2.d0


end subroutine rigid_solid_rhs


! Initialize insect free flight dynamics solver.
! Task: set up the current state of the insect in INSECT%STATE
subroutine rigid_solid_init(Insect)
  use mpi
  use vars
  implicit none

  type(diptera),intent(inout)::Insect
  real (kind=pr) :: yaw,pitch,roll
  integer :: mpicode
  real (kind=pr), dimension(0:3) :: ep
  real (kind=pr), dimension(1:3,1:3) :: C_mat


  ! free flight solver based on quaternions. the task here is to initialize
  ! the "attitude" quaternion (Insect%quaternion) from yaw, pitch and roll
  ! angles. Note that for the free flight solver, this is the last time
  ! body yaw,pitch,roll are meaningful. from now on, quaternions are used

  Insect%gamma = deg2rad(0.0d0)
  Insect%beta  = deg2rad(0.0d0)
  Insect%psi   = deg2rad(0.0d0)

  Insect%rot_body = 0.d0

  yaw   =  Insect%gamma / 2.d0
  pitch =  Insect%beta  / 2.d0
  roll  =  Insect%psi   / 2.d0

  ep(0) = cos(roll)*cos(pitch)*cos(yaw) + sin(roll)*sin(pitch)*sin(yaw)
  ep(1) = sin(roll)*cos(pitch)*cos(yaw) - cos(roll)*sin(pitch)*sin(yaw)
  ep(2) = cos(roll)*sin(pitch)*cos(yaw) + sin(roll)*cos(pitch)*sin(yaw)
  ep(3) = cos(roll)*cos(pitch)*sin(yaw) - sin(roll)*sin(pitch)*cos(yaw)

  ! to draw the very first insect, we also need the rotation matix
  call rotation_matrix_from_quaternion( ep,Insect%M_body_quaternion )

!!!!!!!!!!
Insect%xc_body=(/x0,y0,z0/)
Insect%vc_body=(/0.0,0.0,0.0/)
!!!!!!!!!!

  ! Assemble the state vector
  Insect%STATE = (/ Insect%xc_body(1),Insect%xc_body(2),Insect%xc_body(3), &
                    Insect%vc_body(1),Insect%vc_body(2),Insect%vc_body(3), &
                   ep(0),ep(1),ep(2),ep(3),&
                   Insect%rot_body(1),Insect%rot_body(2),Insect%rot_body(3) /)

  Insect%RHS_this = 0.0
  Insect%RHS_old = 0.0

  write (*,*) "Insect%STATE", Insect%STATE


end subroutine rigid_solid_init
