!-------------------------------------------------------------------------------
! Rigid solid time stepping routines
! AB2 method with Euler startup
!-------------------------------------------------------------------------------
subroutine rigid_solid_time_step(time,dt0,dt1,it,Insect)
  use vars
  implicit none

  real (kind=pr),intent (in) :: time,dt1,dt0
  type(diptera),intent(inout) :: Insect
  integer,intent (in) :: it
  real (kind=pr) :: b10,b11

  ! select scheme
  if(it == 0) then
    ! EULER startup scheme
    ! compute rhs at this time step (updates Insect%RHS_this)
    call rigid_solid_rhs(time,dt1,it,Insect)

    ! Euler step
    Insect%STATE = Insect%STATE + dt1*Insect%RHS_this
    if (Insect%wing_fsi == "feathering") &
      Insect%WING_STATE = Insect%WING_STATE + dt1*Insect%WING_RHS_this
  else
    ! ADAMS-BASHFORTH 2 scheme
    ! update vectors
    Insect%RHS_OLD = Insect%RHS_this
    if (Insect%wing_fsi == "feathering") &
      Insect%WING_RHS_OLD = Insect%WING_RHS_this

    ! compute rhs at this time step (updates Insect%RHS_this)
    call rigid_solid_rhs(time,dt1,it,Insect)

    ! adaptive AB2 coefficients
    b10 = dt1/dt0*(0.5*dt1 + dt0)
    b11 = -0.5*dt1*dt1/dt0

    ! adaptive AB2 step
    Insect%STATE = Insect%STATE + b10*Insect%RHS_this + b11*Insect%RHS_old
    if (Insect%wing_fsi == "feathering") &
      Insect%WING_STATE = Insect%WING_STATE &
      + b10*Insect%WING_RHS_this + b11*Insect%WING_RHS_old
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
subroutine rigid_solid_rhs(time,dt,it,Insect)
  use mpi
  use vars
  implicit none

  integer, intent(in) :: it
  real (kind=pr), intent (in) :: time,dt
  type(diptera),intent(inout) :: Insect
  real (kind=pr) :: m,g,Jx,Jy,Jz,s
  real (kind=pr) :: Maero_l,Jyy_l,Jxy_l,stif_l,damp_l,ang0_l
  real (kind=pr) :: Maero_r,Jyy_r,Jxy_r,stif_r,damp_r,ang0_r
  real (kind=pr) :: phi_dt_i1,phi_dt_i2,theta_dt_i1,theta_dt_i2
  real (kind=pr) :: phi_l,theta_l,phi_dt_l,theta_dt_l,phi_dt2_l,theta_dt2_l
  real (kind=pr) :: phi_r,theta_r,phi_dt_r,theta_dt_r,phi_dt2_r,theta_dt2_r
  real (kind=pr), dimension(0:3) :: ep
  real (kind=pr), dimension(1:3) :: ROT,torque_body
  real (kind=pr), dimension(1:3,1:3) :: C_mat
  real (kind=pr), dimension(1:3) :: mom_aero_l,mom_aero_r
  real (kind=pr), dimension(1:5) :: foo

  !-------------------------------------------------------------------------------
  ! body dynamics
  !-------------------------------------------------------------------------------
  ! initialization
  Insect%RHS_this=0.d0

  if ((root).and.(Insect%BodyMotion/="free_flight")) then
    call abort(9030,"Insect%BodyMotion"//trim(adjustl(Insect%BodyMotion))//" but using free-flight?")
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

  !-------------------------------------------------------------------------------
  ! wing dynamics
  !-------------------------------------------------------------------------------
  if (Insect%wing_fsi == "feathering") then
    ! initialize parameters of the feathering rotation model.
    ! at this point we assume that both wings have the same properties
    Jyy_l = Insect%Jyy
    Jxy_l = Insect%Jxy
    stif_l = Insect%stif
    damp_l = Insect%damp
    ang0_l = Insect%ang0
    Jyy_r = Insect%Jyy
    Jxy_r = Insect%Jxy
    stif_r = Insect%stif
    damp_r = Insect%damp
    ang0_r = Insect%ang0

    ! Left wing
    ! read phi_l, theta_l, phi_dt_l, theta_dt_l,
    ! and data for evaluating phi_dt2_l numerically using FD2
    call FlappingMotion ( time+2.0d0*dt, Insect, Insect%FlappingMotion_left, &
    foo(1),foo(2),foo(3),phi_dt_i2,foo(4),theta_dt_i2,Insect%kine_wing_l )
    call FlappingMotion ( time+dt, Insect, Insect%FlappingMotion_left, &
    foo(1),foo(2),foo(3),phi_dt_i1,foo(4),theta_dt_i1,Insect%kine_wing_l )
    call FlappingMotion ( time, Insect, Insect%FlappingMotion_left, &
    phi_l,foo(2),theta_l,phi_dt_l,foo(4),theta_dt_l,Insect%kine_wing_l )
    ! second derivative of the positional angle, which is prescribed
    phi_dt2_l = (-phi_dt_i2+4.0d0*phi_dt_i1-3.0d0*phi_dt_l) / (2.0d0*dt)
    theta_dt2_l = (-theta_dt_i2+4.0d0*theta_dt_i1-3.0d0*theta_dt_l) / (2.0d0*dt)

    ! Right wing
    ! read phi_r, theta_r, phi_dt_r, theta_dt_r,
    ! and data for evaluating phi_dt2_r numerically using FD2
    call FlappingMotion ( time+2.0d0*dt, Insect, Insect%FlappingMotion_right, &
    foo(1),foo(2),foo(3),phi_dt_i2,foo(4),theta_dt_i2,Insect%kine_wing_r )
    call FlappingMotion ( time+dt, Insect, Insect%FlappingMotion_right, &
    foo(1),foo(2),foo(3),phi_dt_i1,foo(4),theta_dt_i1,Insect%kine_wing_r )
    call FlappingMotion ( time, Insect, Insect%FlappingMotion_right, &
    phi_r,foo(2),theta_r,phi_dt_r,foo(4),theta_dt_r,Insect%kine_wing_r )
    ! second derivative of the positional angle, which is prescribed
    phi_dt2_r = (-phi_dt_i2+4.0d0*phi_dt_i1-3.0d0*phi_dt_r) / (2.0d0*dt)
    theta_dt2_r = (-theta_dt_i2+4.0d0*theta_dt_i1-3.0d0*theta_dt_r) / (2.0d0*dt)

    ! Calculate aerodynamic torques of the wings, in the wing reference frames
    call wing_aero_torques( time, Insect, mom_aero_l, mom_aero_r )
    ! RHS of the passive feathering equations
    ! Use only the 2nd (y) component of the aerodynamic torques
    ! WARNING: Note the sign of theta changed. In the passive feathering equation,
    ! positive theta is up.
    call rhs_alpha(mom_aero_l(2),Jyy_l,Jxy_l,stif_l,damp_l,ang0_l,&
                   phi_l,-theta_l,phi_dt_l,-theta_dt_l,phi_dt2_l,-theta_dt2_l,&
                   Insect%WING_STATE(1:2),Insect%WING_RHS_this(1:2))
    call rhs_alpha(-mom_aero_r(2),Jyy_r,Jxy_r,stif_r,damp_r,ang0_r,&
                   phi_r,-theta_r,phi_dt_r,-theta_dt_r,phi_dt2_r,-theta_dt2_r,&
                   Insect%WING_STATE(3:4),Insect%WING_RHS_this(3:4))
    ! Apply startup conditionning on accelerations
    Insect%WING_RHS_this(2) = Insect%WING_RHS_this(2) * s
    Insect%WING_RHS_this(4) = Insect%WING_RHS_this(4) * s
  endif

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

  ! body dynamics
  if (inicond(1:8)=="backup::") then
    ! resuming the rigid solid solver from a backup
    if (mpirank==0) write(*,*) "Rigid solid solver is resuming from backup..."
    if (mpirank==0) write(*,*) inicond(9:len_trim(inicond))//".rigidsolver"

    if (mpirank==0) then
      open(10, file=inicond(9:len_trim(inicond))//".rigidsolver", &
          form='formatted', status='old')
      ! Select rigid model
      if (Insect%wing_fsi == "feathering") then
        ! With passive feathring rotation
        read(10, *) Insect%time, Insect%STATE, Insect%RHS_old, &
            Insect%RHS_this, Insect%M_body_quaternion, &
            Insect%WING_STATE, Insect%WING_RHS_old, Insect%WING_RHS_this
      else
        ! Without passive feathering rotation
        read(10, *) Insect%time, Insect%STATE, Insect%RHS_old, &
            Insect%RHS_this, Insect%M_body_quaternion
      endif
      close(10)
    endif

    ! Broadcast to all ranks
    call MPI_BCAST( Insect%time, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%STATE, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%RHS_this, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%RHS_old, 13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,1), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,2), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_BCAST( Insect%M_body_quaternion(:,3), 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    if (Insect%wing_fsi == "feathering") then
      call MPI_BCAST( Insect%WING_STATE, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
      call MPI_BCAST( Insect%WING_RHS_this, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
      call MPI_BCAST( Insect%WING_RHS_old, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    endif

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


    ! assemble the state vector
    Insect%STATE = (/ Insect%xc_body(1),Insect%xc_body(2),Insect%xc_body(3), &
    Insect%vc_body(1),Insect%vc_body(2),Insect%vc_body(3), &
    ep(0),ep(1),ep(2),ep(3),&
    Insect%rot_body(1),Insect%rot_body(2),Insect%rot_body(3) /)

    Insect%RHS_this = 0.0
    Insect%RHS_old = 0.0

    ! wing dynamics
    if (Insect%wing_fsi == "feathering") then
      ! Fetch the prescribed wing kinematics
      call FlappingMotion_right (time, Insect)
      call FlappingMotion_left (time, Insect)
      ! Assemble the state vector
      Insect%WING_STATE = (/ Insect%alpha_l,Insect%alpha_dt_l, &
                             Insect%alpha_r,Insect%alpha_dt_r /)
      Insect%WING_RHS_this = 0.0
      Insect%WING_RHS_old = 0.0
    endif

  endif
  if(root) write (*,*) "Insect%STATE", Insect%STATE
  if(root) write (*,*) "Insect%WING_STATE", Insect%WING_STATE

end subroutine rigid_solid_init


!-------------------------------------------------------------------------------
! Evaluate RHS of passive feathering rotation equations for one wing.
!-------------------------------------------------------------------------------
subroutine rhs_alpha(Maero,Jyy,Jxy,stif,damp,ang0,phi,theta,phi_dt,theta_dt,phi_dt2,theta_dt2,u,du)
  implicit none

  real(kind=pr),intent(in) :: Maero,Jyy,Jxy,stif,damp,ang0,phi,theta,phi_dt,theta_dt,phi_dt2,theta_dt2,u(1:2)
  real(kind=pr),intent(out) :: du(1:2)
  real(kind=pr) :: Miner

!!! THIS FORMULATION ASSUMES THAT THE BODY DOES NOT MOVE !!!

  ! Inertial torque
!  Miner = Jxy*phi_dt2*dcos(u(1)) ! First version
  Miner = Jyy*( 0.5d0*(phi_dt*phi_dt*dcos(theta)*dcos(theta)-theta_dt*theta_dt)*dsin(2.0d0*u(1)) - &
                phi_dt2*dsin(theta) - phi_dt*theta_dt*dcos(theta)*(1.0d0+dcos(2.0d0*u(1))) ) + &
          Jxy*( phi_dt2*dcos(theta)*dcos(u(1)) + theta_dt2*sin(u(1)) + &
                0.5d0*phi_dt*phi_dt*dsin(2.0d0*theta)*dsin(u(1))-2.0d0*theta_dt*phi_dt*dsin(theta)*dcos(u(1)) )
!  Miner = 0

  ! rhs of the pitching equation
  du(1) = u(2)
  du(2) = (Miner+Maero-stif*(u(1)-ang0)-damp*u(2))/Jyy

! debug : only inertia
!du(2) = (Miner-stif*(u(1)-ang0)-damp*u(2))/Jyy

end subroutine rhs_alpha
