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

  ! Periodization: if the insect moves out of the box, eg, x<0.0 or x>xl then
  ! we
  if (periodic) then
    if (Insect%STATE(1)<0.0) Insect%STATE(1) = Insect%STATE(1) + xl
    if (Insect%STATE(2)<0.0) Insect%STATE(2) = Insect%STATE(2) + yl
    if (Insect%STATE(3)<0.0) Insect%STATE(3) = Insect%STATE(3) + zl

    if (Insect%STATE(1)>=xl) Insect%STATE(1) = Insect%STATE(1) - xl
    if (Insect%STATE(2)>=yl) Insect%STATE(2) = Insect%STATE(2) - yl
    if (Insect%STATE(3)>=zl) Insect%STATE(3) = Insect%STATE(3) - zl
  endif

  if (mpirank==0) then
    open  (17,file='rigidsolidsolver.t',status='unknown',position='append')
    write (17,'(14(es15.8,1x))') time, Insect%STATE
    close (17)
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
  real (kind=pr) :: m,g, Jx, Jy,Jz, s
  real (kind=pr), dimension(0:3) :: ep
  real (kind=pr), dimension(1:3) :: ROT, torque_body
  real (kind=pr), dimension(1:3,1:3) :: C_mat

  ! initialization
  Insect%RHS_this=0.d0

  if ((root).and.(Insect%BodyMotion/="free_flight")) then
    call abort("Insect%BodyMotion"//trim(adjustl(Insect%BodyMotion))//" but using free-flight?")
  endif

  ! copy some shortcuts (this is easier to code)
  g  = Insect%gravity
  m  = Insect%mass
  Jx = Insect%Jroll_body
  Jy = Insect%Jpitch_body
  Jz = Insect%Jyaw_body

  ! extract rotation (unit) quaternion from state vector, and create the rotation
  ! matrix from it.
  ep = Insect%STATE(7:10)
  call rotation_matrix_from_quaternion( ep,Insect%M_body_quaternion )
  ! The equations of motion for the rotation are written in the body reference
  ! frame, thus the fluid torque has to be transformed to the body system
  torque_body = matmul( Insect%M_body_quaternion, GlobalIntegrals%torque+GlobalIntegrals%torque_unst)

  ! extract angular velocity vector from state vector
  ROT = Insect%STATE(11:13)

  ! startup conditioner (to avoid problems with impulsively started motion)
  if (Insect%startup_conditioner=="yes") then
    s = startup_conditioner(time, 0.1d0, 0.5d0)
  else
    s = 1.d0
  endif

  !-----------------------------------------------------------------------------
  ! Integrate the 13 equations of motion (written in quaternion formulation)
  ! The underlying eqns can be found in
  ! M. Maeda et al. (2012) A free-flight Simulation of Insect flapping flight
  ! J. Aero Aqua Bio-Mech(1):1,71-79
  !-----------------------------------------------------------------------------
  ! To avoid confusion with the famous unsteady corrections, the actual translation eqn reads:
  !   rho_s*V*u_dot = (rho_s-rho_f)*V*g + F_fluid
  ! but the latter term is
  !   rho_s*V*u_dot = (rho_s-rho_f)*V*g + Integral(Penal) + V*rho_f*u_dot
  ! so we can put this, as uhlmann does, on the right hand side:
  !   (rho_s-rho_f)*V*u_dot = (rho_s-rho_f)*V*g + Integral(Penal)
  ! We have thus two options regarding this, either
  !   m_corrected * u_dot = m_corrected * g + Integral(penal)       [implicit unst corrections]
  ! or
  !   m * u_dot = m_corrected * g + Integral(penal) + force_unst    [explicit unst corrections]
  ! the last equation can also be rewritten using twice the same m
  !   m * u_dot = m * g_corrected + Integral(penal) + force_unst    [explicit unst corrections]
  ! where g_corrected = (rho_s-rho_f)/rho_s * g
  !-----------------------------------------------------------------------------
  ! For the torque, the unsteady corrections are explicitly taken into account
  ! above, thus cal_unst=1; is required in the parameter file.
  !-----------------------------------------------------------------------------
  ! integrate coordinates (dx/dt = vx) Note: this is in global reference frame
  Insect%RHS_this(1) = Insect%STATE(4)
  Insect%RHS_this(2) = Insect%STATE(5)
  Insect%RHS_this(3) = Insect%STATE(6)
  ! integrate velocities (dvx/dt = F) Note: this is in global reference frame
  Insect%RHS_this(4) = s*(GlobalIntegrals%force(1)+GlobalIntegrals%force_unst(1))/m
  Insect%RHS_this(5) = s*(GlobalIntegrals%force(2)+GlobalIntegrals%force_unst(2))/m
  Insect%RHS_this(6) = s*(GlobalIntegrals%force(3)+GlobalIntegrals%force_unst(3))/m + g
  ! integrate quaternion attitudes
  Insect%RHS_this(7)  = 0.5d0*(-ep(1)*ROT(1)-ep(2)*ROT(2)-ep(3)*ROT(3))
  Insect%RHS_this(8)  = 0.5d0*(+ep(0)*ROT(1)-ep(3)*ROT(2)+ep(2)*ROT(3))
  Insect%RHS_this(9)  = 0.5d0*(+ep(3)*ROT(1)+ep(0)*ROT(2)-ep(1)*ROT(3))
  Insect%RHS_this(10) = 0.5d0*(-ep(2)*ROT(1)+ep(1)*ROT(2)+ep(0)*ROT(3))
  ! integrate angular velocities
  Insect%RHS_this(11) = ( (Jy-Jz)*ROT(2)*ROT(3) + s*torque_body(1) )/Jx
  Insect%RHS_this(12) = ( (Jz-Jx)*ROT(3)*ROT(1) + s*torque_body(2) )/Jy
  Insect%RHS_this(13) = ( (Jx-Jy)*ROT(1)*ROT(2) + s*torque_body(3) )/Jz

  ! turn on or off degrees of freedom for free flight solver. The string from
  ! ini file contains 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll
  ! degrees of freedom by multiplying the respective RHS by zero, keeping the
  ! value thus constant
  Insect%RHS_this(4) = Insect%RHS_this(4) * Insect%DoF_on_off(1)   ! x translation
  Insect%RHS_this(5) = Insect%RHS_this(5) * Insect%DoF_on_off(2)   ! y translation
  Insect%RHS_this(6) = Insect%RHS_this(6) * Insect%DoF_on_off(3)   ! z translation
  Insect%RHS_this(13) = Insect%RHS_this(13) * Insect%DoF_on_off(4) ! yaw rotation
  Insect%RHS_this(12) = Insect%RHS_this(12) * Insect%DoF_on_off(5) ! pitch rotation
  Insect%RHS_this(11) = Insect%RHS_this(11) * Insect%DoF_on_off(6) ! roll rotation

  ! Insect%RHS_this(4:6) = Insect%RHS_this(4:6) * startup_conditioner(time,0.0d0,0.25d0)
  ! Insect%RHS_this(11:13) = Insect%RHS_this(11:13) * startup_conditioner(time,0.0d0,0.25d0)

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
  if (mpirank==0) write(*,*) "rigid solid init at time=", Insect%time

  if (inicond(1:8)=="backup::") then
    ! resuming the rigid solid solver from a backup
    if (mpirank==0) write(*,*) "Rigid solid solver is resuming from backup..."
    if (mpirank==0) write(*,*) inicond(9:len_trim(inicond))//".rigidsolver"

    if (mpirank==0) then
      open(10, file=inicond(9:len_trim(inicond))//".rigidsolver", &
          form='formatted', status='old')
      read(10, *) Insect%time, Insect%STATE, Insect%RHS_old, &
          Insect%RHS_this, Insect%M_body_quaternion
      close(10)
    endif

    call MPI_BCAST( Insect%time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%STATE, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%RHS_this, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%RHS_old, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,2), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,3), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

  else
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
  endif
  if(root) write (*,*) "Insect%STATE", Insect%STATE


end subroutine rigid_solid_init
