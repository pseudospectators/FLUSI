!-------------------------------------------------------------------------------
! 2nd generation insect module
!-------------------------------------------------------------------------------
! contains now motion of the body (translation + yaw,pitch,roll
! contains a stroke plane angle for direct implementation of literature results
!-------------------------------------------------------------------------------

module insect_module
  use vars
  use helpers
  implicit none

  ! this will hold the surface markers and their normals used for particles:
  real(kind=pr), allocatable, dimension(:,:) :: particle_points
  ! variables to decide whether to draw the body or not.
  logical :: body_already_drawn = .false.
  character(len=strlen) :: body_moves="yes"
  ! do we use periodicity?
  logical, parameter :: periodic = .true.
  ! arrays for fourier coefficients are fixed size (avoiding issues with allocatable
  ! elements in derived datatypes) this is their length:
  integer, parameter :: nfft_max = 1024


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
    real(kind=pr), dimension (:), allocatable :: vec_t, &
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
    real(kind=pr), dimension(1:3) :: xc_body, vc_body
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
    real(kind=pr), dimension(1:3) :: rot_l, rot_r, rot_body
    ! angular acceleration vectors (wings L+R)
    real(kind=pr), dimension(1:3) :: rot_dt_l, rot_dt_r
    ! Angular velocities of wings and body in global system
    real(kind=pr), dimension(1:3) :: rot_body_glob, rot_l_glob, rot_r_glob
    ! Vector from body centre to pivot points in global reference frame
    real(kind=pr), dimension(1:3) :: x_pivot_l_glob, x_pivot_r_glob
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
    real(kind=pr), dimension(1:3,1:3) :: M_body_quaternion
    real(kind=pr), dimension(1:13) :: RHS_old, RHS_this, STATE
    real(kind=pr), dimension(1:6) :: DoF_on_off
    character(len=strlen) :: startup_conditioner

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
    ! wing inertia
    real(kind=pr) :: Jxx,Jyy,Jzz,Jxy

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

  real(kind=pr), save :: smoothing

contains


  !---------------------------------------
  ! note these include files also have to be specified as dependencies in the
  ! Makefile for make to check if one of them changed
  include "periodization.f90"
  include "body_geometry.f90"
  include "body_motion.f90"
  include "rigid_solid_time_stepper.f90"
  include "wings_geometry.f90"
  include "wings_motion.f90"
  include "stroke_plane.f90"
  include "kineloader.f90"
  !---------------------------------------


  !-------------------------------------------------------------------------------
  ! Main routine for drawing insects. Loops over the entire domain, computes
  ! coordinates in various systems (global-, body-, stroke-, wing-) and calls
  ! subroutines doing the actual job of defining the mask. Note all surfaces are
  ! smoothed.
  !-------------------------------------------------------------------------------
  subroutine Draw_Insect ( time, Insect, mask, mask_color, us)
    use vars
    implicit none

    real(kind=pr), intent(in) :: time
    type(diptera),intent(inout) :: Insect
    real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

    real(kind=pr),dimension(1:3) :: x, x_body, x_wing_l, x_wing_r, x_eye_r,&
    x_eye_l,x_head, v_tmp, rot_b_psi,rot_b_beta,rot_b_gamma
    real(kind=pr), dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, &
    M1_b, M2_b, M3_b, M1, M2, M3, M_stroke_l, M_stroke_r, M_body_inv, &
    M_wing_l_inv, M_wing_r_inv
    real(kind=pr)::t1
    integer :: ix, iy, iz
    integer, save :: counter = 0
    integer(kind=2) :: color_body, color_l, color_r

    !-----------------------------------------------------------------------------
    ! colors for Diptera (one body, two wings)
    !-----------------------------------------------------------------------------
    color_body = 1
    color_l = 2
    color_r = 3

    if ((dabs(Insect%time-time)>1.0d-10).and.(mpirank==0).and.(Insect%BodyMotion=="free_flight")) then
      write (*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
    endif

    !-----------------------------------------------------------------------------
    ! delete old mask
    !-----------------------------------------------------------------------------
    if (body_moves=="no") then
      ! the body is at rest, so we will not draw it. Delete everything EXCEPT the
      ! body, which is marked by its specific color
      where (mask_color/=color_body)
        mask = 0.d0
        mask_color = 0
      end where
      ! the value in the mask array is divided by eps already and will be divided
      ! by eps again in the main mask wrapper. we thus multiply the existing body
      ! by eps.
      where (mask_color==color_body)
        mask = mask*eps
      end where
      ! as the body rests it has no solid body velocity, which means we can safely
      ! reset the velocity everywhere (this step is actually unnessesary, but for
      ! safety we do it as well)
      us = 0.d0
    else
      ! the body of the insect moves, so we will construct the entire insect in this
      ! (and all other) call, and therefore we can safely reset the entire mask to zeros.
      mask = 0.d0
      mask_color = 0
      us = 0.d0
    endif

    if (nx/=1) then
      Insect%safety = 2.5d0*max(dz,dy,dx)
      Insect%smooth = 1.0d0*max(dz,dy,dx)
    else
      Insect%safety = 2.5d0*max(dz,dy)
      Insect%smooth = 1.0d0*max(dz,dy)
    endif
    smoothing = Insect%smooth

    ! some checks
    if ((mpirank==0).and.((iMoving.ne.1).or.(iPenalization.ne.1))) then
      write (*,*) "insects.f90::DrawInsect: the parameters iMoving or iPenalization are wrong."
      call abort()
    endif

    !-----------------------------------------------------------------------------
    ! fetch current motion state
    !-----------------------------------------------------------------------------
    call BodyMotion (time, Insect)
    call FlappingMotion_right (time, Insect)
    call FlappingMotion_left (time, Insect)
    call StrokePlane (time, Insect)

    !-----------------------------------------------------------------------------
    ! define the rotation matrices to change between coordinate systems
    !-----------------------------------------------------------------------------
    if (Insect%BodyMotion=="free_flight") then
    !   call Rx(M1_b,Insect%psi)
    !   call Ry(M2_b,Insect%beta)
    !   call Rz(M3_b,Insect%gamma)
    !   M_body = matmul(M1_b,matmul(M2_b,M3_b))
    !   if(root) then
    !   write(*,*) "--rotmat-"
    !   write(*,'(3(es12.4,1x))') M_body(1,:)
    !   write(*,'(3(es12.4,1x))') M_body(2,:)
    !   write(*,'(3(es12.4,1x))') M_body(3,:)
    !   write(*,*) "--quat-"
    !   write(*,'(3(es12.4,1x))') Insect%M_body_quaternion(1,:)
    !   write(*,'(3(es12.4,1x))') Insect%M_body_quaternion(2,:)
    !   write(*,'(3(es12.4,1x))') Insect%M_body_quaternion(3,:)
    !   write(*,*) "---"
    ! endif
      !
      ! ! QUATERNIONS
      M_body = Insect%M_body_quaternion
      !
      ! if(maxval(M_body)==0.d0) call abort()

    else
      ! conventional yaw, pitch, roll. Note the order of matrices is important.
      ! first we yaw, then we pitch, then we roll the insect. Note that when the
      ! free-flight solver is used, this matrix is obtained from quaternions, and
      ! not as a product of simple rotaion matrices. The latter can cause "gimbal-lock"
      call Rx(M1_b,Insect%psi)
      call Ry(M2_b,Insect%beta)
      call Rz(M3_b,Insect%gamma)
      M_body = matmul(M1_b,matmul(M2_b,M3_b))
    endif

    call Ry(M1,Insect%eta_stroke)
    M_stroke_l = M1

    call Rx(M1,pi)
    call Ry(M2,Insect%eta_stroke)
    M_stroke_r = matmul(M1,M2)

    call Ry(M1,Insect%alpha_l)
    call Rz(M2,Insect%theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,Insect%phi_l)
    M_wing_l = matmul(M1,matmul(M2,matmul(M3,M_stroke_l)))

    ! note the coordinate system is rotated so we don't need to inverse the sign
    ! of theta, and the wings still rotate in opposite direction
    call Ry(M1,-Insect%alpha_r)
    call Rz(M2, Insect%theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,-Insect%phi_r)
    M_wing_r = matmul(M1,matmul(M2,matmul(M3,M_stroke_r)))

    !-----------------------------------------------------------------------------
    ! angular velocity vectors
    !-----------------------------------------------------------------------------
    ! If Insect%BodyMotion=="free_flight", then
    ! the angular velocity is a state variable of the insect free flight
    ! solver, and the value is copied in body-motion to the insect variable
    ! in line 61 of body_motion.f90.
    ! Otherwise, the following block is active:
    if (Insect%BodyMotion/="free_flight") then
      ! NOTE: when using the quaternion based free-flight solver, the angular
      ! velocity of the body is computed dynamically, and the rotation matrix that
      ! brings us from global to body system is computed with quaternions
      rot_b_psi   = (/ Insect%psi_dt, 0.0d0, 0.0d0   /)
      rot_b_beta  = (/ 0.0d0, Insect%beta_dt, 0.0d0  /)
      rot_b_gamma = (/ 0.0d0, 0.0d0, Insect%gamma_dt /)

      ! body angular velocity in the body coodrinate system
      Insect%rot_body = matmul(M_body,matmul(transpose(M3_b),rot_b_gamma+ &
      matmul(transpose(M2_b),rot_b_beta+matmul(transpose(M1_b),rot_b_psi))))
    endif

    ! wing angular velocities in the wing coordinate system
    call angular_velocities ( time, Insect )

    !-----------------------------------------------------------------------------
    ! angular acceleration for wings (required for inertial power)
    !-----------------------------------------------------------------------------
    call angular_accel( time, Insect )

    !-----------------------------------------------------------------------------

    ! inverse of the rotation matrices
    !-----------------------------------------------------------------------------
    M_body_inv = transpose(M_body)
    M_wing_l_inv = transpose(M_wing_l)
    M_wing_r_inv = transpose(M_wing_r)

    !-----------------------------------------------------------------------------
    ! angular velocity in the global reference frame
    !-----------------------------------------------------------------------------
    Insect%rot_body_glob = matmul(M_body_inv, Insect%rot_body)
    Insect%rot_l_glob = matmul(M_body_inv, matmul(M_wing_l_inv, Insect%rot_l) + Insect%rot_body)
    Insect%rot_r_glob = matmul(M_body_inv, matmul(M_wing_r_inv, Insect%rot_r) + Insect%rot_body)

    !-----------------------------------------------------------------------------
    ! vector from body centre to left/right pivot point in global reference frame,
    ! for aerodynamic power
    !-----------------------------------------------------------------------------
    Insect%x_pivot_l_glob = matmul(M_body_inv, Insect%x_pivot_l)
    Insect%x_pivot_r_glob = matmul(M_body_inv, Insect%x_pivot_r)


    !-----------------------------------------------------------------------------
    ! write kinematics to disk (Dmitry, 28 Oct 2013)
    ! do so only every itdrag time steps (Thomas, 8 Jul 2014)
    !-----------------------------------------------------------------------------
    counter = counter + 1
    if ((mpirank == 0).and.(mod(counter,itdrag)==0)) then
      open  (17,file='kinematics.t',status='unknown',position='append')
      write (17,'(26(es15.8,1x))') time, Insect%xc_body, Insect%psi, Insect%beta, &
      Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
      Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
      Insect%rot_l, Insect%rot_r, Insect%rot_dt_l, Insect%rot_dt_r
      close (17)
    endif

    !-----------------------------------------------------------------------------
    ! Draw indivudual parts of the Diptera. Separate loops are faster
    ! since the compiler can optimize them better
    !-----------------------------------------------------------------------------
    ! BODY. Now the body is special: if the insect does not move (or rotate), the
    ! body does not change in time. On the other hand, it is quite expensive to
    ! compute, since it involves a lot of points (volume), and it is a source of
    ! load balancing problems, since many cores do not draw the body at all.
    ! We thus try to draw it only once and then simply not to erase it later.
    !-----------------------------------------------------------------------------
    if (body_moves=="no") then
      if (body_already_drawn .eqv. .false.) then
        ! the body is at rest, but it is the first call to this routine, so
        ! draw it now.
        if (root) write(*,*) "Flag body_moves is no and we did not yet draw"
        if (root) write(*,*) "the body once: we do that now, and skip draw_body"
        if (root) write(*,*) "from now on."
        call draw_body( mask, mask_color, us, Insect, color_body, M_body)
        body_already_drawn = .true.
      endif
    else
      ! the body moves, draw it
      call draw_body( mask, mask_color, us, Insect, color_body, M_body)
    endif

    !-----------------------------------------------------------------------------
    ! Wings
    !-----------------------------------------------------------------------------
    if (Insect%RightWing == "yes") then
      call draw_wing(mask,mask_color,us,Insect,color_r,M_body,M_wing_r,&
      Insect%x_pivot_r,Insect%rot_r )
    endif

    if (Insect%LeftWing == "yes") then
      call draw_wing(mask,mask_color,us,Insect,color_l,M_body,M_wing_l,&
      Insect%x_pivot_l,Insect%rot_l )
    endif

    !-----------------------------------------------------------------------------
    ! Add solid body rotation (i.e. the velocity field that originates
    ! from the body rotation and translation. Until now, the wing velocities
    ! were the only ones set plus they are in the body reference frame
    !-----------------------------------------------------------------------------
    t1 = MPI_wtime()
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
          x = periodize_coordinate(x - Insect%xc_body)
          x_body = matmul(M_body,x)
          ! add solid body rotation in the body-reference frame, if color
          ! indicates that this part of the mask belongs to the insect
          if ((mask(ix,iy,iz) > 0.d0).and.(mask_color(ix,iy,iz)>0)) then

            ! translational part. we compute the rotational part in the body
            ! reference frame, therefore, we must transform the body translation
            ! velocity Insect%vc (which is in global coordinates) to the body frame
            v_tmp = matmul(M_body,Insect%vc_body)

            ! add solid body rotation to the translational velocity field. Note
            ! that rot_body and x_body are in the body reference frame
            v_tmp(1) = v_tmp(1)+Insect%rot_body(2)*x_body(3)-Insect%rot_body(3)*x_body(2)
            v_tmp(2) = v_tmp(2)+Insect%rot_body(3)*x_body(1)-Insect%rot_body(1)*x_body(3)
            v_tmp(3) = v_tmp(3)+Insect%rot_body(1)*x_body(2)-Insect%rot_body(2)*x_body(1)

            ! the body motion is added to the wing motion, which is already in us
            ! and they are also in the body refrence frame. However, us has to be
            ! in the global reference frame, so M_body_inverse is applied
            us(ix,iy,iz,1:3) = matmul( M_body_inv, us(ix,iy,iz,1:3)+v_tmp )
          endif
        enddo
      enddo
    enddo
    time_insect_vel = time_insect_vel + MPI_wtime() - t1
  end subroutine Draw_Insect


  !-------------------------------------------------------
  ! short for the smooth step function.
  ! the smooting is defined in Insect%smooth, here we need only x, and the
  ! thickness (i.e., in the limit, steps=1 if x<t and steps=0 if x>t
  !-------------------------------------------------------
  real(kind=pr) function steps(x,t)
    use vars
    implicit none
    real(kind=pr) :: f,x,t
    call smoothstep(f,x,t,smoothing)
    steps=f
  end function


  !-------------------------------------------------------
  ! Compute angle from coefficients provided by Maeda
  !-------------------------------------------------------
  subroutine get_dangle( angles, F, a, b, shift_phase, initial_phase, dangle, dangle_dt )
    use vars
    implicit none
    integer, intent(in) :: F  ! wavenumber (Dmitry, 7 Nov 2013)
    real(kind=pr), intent(in) :: angles ! 2*pi*F*time (Dmitry, 7 Nov 2013)
    real(kind=pr), intent(in) :: a
    real(kind=pr), intent(in) :: b
    real(kind=pr), intent(in) :: shift_phase
    real(kind=pr), intent(in) :: initial_phase
    real(kind=pr), intent(out) :: dangle
    real(kind=pr), intent(out) :: dangle_dt ! velocity increment (Dmitry, 7 Nov 2013)
    real(kind=pr) :: dAmp
    real(kind=pr) :: factor_amp = 1.0d0  ! Dmitry, 7 Nov 2013
    real(kind=pr) :: phase
    !!----------------------------

    !! d_amplitude
    dAmp = dsqrt(a**2 +b**2)*factor_amp

    !! phase
    if( b>0.0d0 ) then
      phase = datan(a/b)
    elseif( b<0.0d0 ) then
      phase = datan(a/b) +pi
    else !! b == 0 -> avoid division by zero
      phase = pi*0.5d0 !! sin(PI/2) = cos
    endif

    phase = phase + (shift_phase +initial_phase*2.0d0*pi)*dble(F)

    !! d_angle
    dangle = dAmp*dsin( angles +phase )

    !! velocity increment (Dmitry, 7 Nov 2013)
    dangle_dt = 2.0d0*pi*dble(F) * dAmp*dcos( angles +phase )

    return
  end subroutine get_dangle



  ! Compute aerodynamic power
  subroutine aero_power(Insect,apowtotal)
    use vars
    implicit none

    integer :: color_body, color_l, color_r
    real(kind=pr), dimension(1:3) :: omrel, momrel
    real(kind=pr), intent(out) :: apowtotal
    type(diptera),intent(inout)::Insect

    ! colors for Diptera (one body, two wings)
    color_body = 1
    color_l = 2
    color_r = 3

    ! body is not driven directly, therefore the power is set to zero
    Insect%PartIntegrals(color_body)%APow = 0.0d0

    ! left wing
    ! relative angular velocity
    omrel = Insect%rot_l_glob - Insect%rot_body_glob

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_l)%Torque + &
    Insect%PartIntegrals(color_l)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_l)%APow = - sum( momrel * omrel )

    ! right wing
    ! relative angular velocity
    omrel = Insect%rot_r_glob - Insect%rot_body_glob

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_r)%Torque + &
    Insect%PartIntegrals(color_r)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_r)%APow = - sum( momrel * omrel )

    ! Total aerodynamic power
    apowtotal = Insect%PartIntegrals(color_body)%APow + &
    Insect%PartIntegrals(color_l)%APow + Insect%PartIntegrals(color_r)%APow

  end subroutine aero_power


  !-------------------------------------------------------------------------------
  ! Compute interial power, i.e. the power the insect would have to invest
  ! when flapping its wings in vacuum.
  ! OUTPUT:
  !       ipowtotal: total inertial power
  !       Insect%PartIntegrals%IPow: (global): individual inertial power
  ! INPUT:
  !       Insect%rot_dt_l (global): left wing angular acceleration
  !       Insect%rot_dt_r (global): right wing angular acceleration
  !       Insect%Jxx,Jyy,Jxy,Jzz (global) Wing inertia
  ! MATHEMATICS:
  !       P_inertia = omega*( J*omega_dt + omega \cross (J*omega) )
  !                 = omega*( a + omega \cross b )
  !       The interia tensor is (it is specified in the PARAMS file)
  !           / Jxx Jxy 0   \
  !       J = | Jxy Jyy 0   |
  !           \ 0   0   Jzz /
  ! SEE ALSO
  !       Berman, Wang: Energy minimizing kinematics in hovering insect flight
  !       (JFM 582, 2007), eqn 2.22 (looks a bit different)
  !-------------------------------------------------------------------------------
  subroutine inert_power(Insect,ipowtotal)
    use vars
    implicit none

    real(kind=pr), intent(out) :: ipowtotal
    real(kind=pr), dimension(1:3) :: a,b
    integer(kind=2) :: color_body, color_l,color_r
    type(diptera),intent(inout)::Insect

    !-- colors for Diptera (one body, two wings)
    color_body = 1
    color_l = 2
    color_r = 3

    !-- LEFT WING
    a(1) = Insect%Jxx * Insect%rot_dt_l(1) + Insect%Jxy * Insect%rot_dt_l(2)
    a(2) = Insect%Jxy * Insect%rot_dt_l(1) + Insect%Jyy * Insect%rot_dt_l(2)
    a(3) = Insect%Jzz * Insect%rot_dt_l(3)

    b(1) = Insect%Jxx * Insect%rot_l(1) + Insect%Jxy * Insect%rot_l(2)
    b(2) = Insect%Jxy * Insect%rot_l(1) + Insect%Jyy * Insect%rot_l(2)
    b(3) = Insect%Jzz * Insect%rot_l(3)

    Insect%PartIntegrals(color_l)%IPow = &
    Insect%rot_l(1) * (a(1)+Insect%rot_l(2)*b(3)-Insect%rot_l(3)*b(2)) +&
    Insect%rot_l(2) * (a(2)+Insect%rot_l(3)*b(1)-Insect%rot_l(1)*b(3)) +&
    Insect%rot_l(3) * (a(3)+Insect%rot_l(1)*b(2)-Insect%rot_l(2)*b(1))

    !-- RIGHT WING
    a(1) = Insect%Jxx * Insect%rot_dt_r(1) + Insect%Jxy * Insect%rot_dt_r(2)
    a(2) = Insect%Jxy * Insect%rot_dt_r(1) + Insect%Jyy * Insect%rot_dt_r(2)
    a(3) = Insect%Jzz * Insect%rot_dt_r(3)

    b(1) = Insect%Jxx * Insect%rot_r(1) + Insect%Jxy * Insect%rot_r(2)
    b(2) = Insect%Jxy * Insect%rot_r(1) + Insect%Jyy * Insect%rot_r(2)
    b(3) = Insect%Jzz * Insect%rot_r(3)

    Insect%PartIntegrals(color_r)%IPow = &
    Insect%rot_r(1) * (a(1)+Insect%rot_r(2)*b(3)-Insect%rot_r(3)*b(2)) +&
    Insect%rot_r(2) * (a(2)+Insect%rot_r(3)*b(1)-Insect%rot_r(1)*b(3)) +&
    Insect%rot_r(3) * (a(3)+Insect%rot_r(1)*b(2)-Insect%rot_r(2)*b(1))

    ipowtotal = Insect%PartIntegrals(color_r)%IPow + Insect%PartIntegrals(color_l)%IPow

  end subroutine inert_power


  !-------------------------------------------------------------------------------
  ! given the angles of each wing (and their time derivatives), compute
  ! the angular velocity vectors for both wings.
  ! output Insect%rot_r, Insect%rot_l is understood in the WING coordinate system.
  !-------------------------------------------------------------------------------
  subroutine angular_velocities ( time, Insect )
    use vars
    implicit none

    real(kind=pr), intent(in) :: time
    type(diptera), intent(inout) :: Insect

    real(kind=pr) :: eta_stroke
    real(kind=pr) :: phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r
    real(kind=pr) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l
    real(kind=pr), dimension(1:3) :: rot_l_alpha, rot_l_theta, rot_l_phi, &
    rot_r_alpha, rot_r_theta, rot_r_phi
    real(kind=pr), dimension(1:3,1:3) :: M_wing_l, M_wing_r, &
    M1_tmp, M2_tmp, M1_l, M2_l, M3_l, M1_r, M2_r, M3_r, &
    M_stroke_l, M_stroke_r


    phi_r      = Insect%phi_r
    alpha_r    = Insect%alpha_r
    theta_r    = Insect%theta_r
    phi_dt_r   = Insect%phi_dt_r
    alpha_dt_r = Insect%alpha_dt_r
    theta_dt_r = Insect%theta_dt_r

    phi_l      = Insect%phi_l
    alpha_l    = Insect%alpha_l
    theta_l    = Insect%theta_l
    phi_dt_l   = Insect%phi_dt_l
    alpha_dt_l = Insect%alpha_dt_l
    theta_dt_l = Insect%theta_dt_l

    eta_stroke = Insect%eta_stroke

    !-----------------------------------------------------------------------------
    ! define the rotation matrices to change between coordinate systems
    !-----------------------------------------------------------------------------
    call Ry(M1_tmp,eta_stroke)
    M_stroke_l = M1_tmp

    call Rx(M1_tmp,pi)
    call Ry(M2_tmp,eta_stroke)
    M_stroke_r = matmul(M1_tmp,M2_tmp)

    call Ry(M1_l,alpha_l)
    call Rz(M2_l,theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3_l,phi_l)
    M_wing_l = matmul(M1_l,matmul(M2_l,matmul(M3_l,M_stroke_l)))

    ! note the coordinate system is rotated so we don't need to inverse the sign
    ! of theta, and the wings still rotate in opposite direction
    call Ry(M1_r,-alpha_r)
    call Rz(M2_r,theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3_r,-phi_r)
    M_wing_r = matmul(M1_r,matmul(M2_r,matmul(M3_r,M_stroke_r)))

    !-----------------------------------------------------------------------------
    ! angular velocity vectors
    !-----------------------------------------------------------------------------
    rot_l_alpha = (/ 0.0d0, alpha_dt_l, 0.0d0 /)
    rot_l_theta = (/ 0.0d0, 0.0d0, theta_dt_l /)
    rot_l_phi   = (/ phi_dt_l, 0.0d0, 0.0d0   /)
    rot_r_alpha = (/ 0.0d0, -alpha_dt_r, 0.0d0/)
    rot_r_theta = (/ 0.0d0, 0.0d0, theta_dt_r /)
    rot_r_phi   = (/ -phi_dt_r, 0.0d0, 0.0d0  /)

    ! in the wing coordinate system
    Insect%rot_l = matmul(M_wing_l,matmul(transpose(M_stroke_l),matmul(transpose(M3_l), &
    rot_l_phi+matmul(transpose(M2_l),rot_l_theta+matmul(transpose(M1_l), &
    rot_l_alpha)))))
    Insect%rot_r = matmul(M_wing_r,matmul(transpose(M_stroke_r),matmul(transpose(M3_r), &
    rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
    rot_r_alpha)))))

  end subroutine angular_velocities



  !-------------------------------------------------------------------------------
  ! Numerically estimate (it's a very precise estimation) the angular acceleration
  ! vectors for both wings, using centered finite differences
  !-------------------------------------------------------------------------------
  subroutine angular_accel( time, Insect )
    use vars
    implicit none
    real(kind=pr), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    type(diptera) :: Insect1, Insect2

    real(kind=pr) :: dt
    real(kind=pr), dimension(1:3) :: rot_l_1, rot_l_2, rot_r_1, rot_r_2
    real(kind=pr) :: eta_stroke
    real(kind=pr) :: phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r
    real(kind=pr) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l

    dt =1.0d-7

    ! note we cannot go to negative times
    if (time > 1.0d-7) then
      Insect1 = Insect
      Insect2 = Insect

      ! fetch motion state at time+dt
      call FlappingMotion_right(time+dt, Insect2)
      call FlappingMotion_left (time+dt, Insect2)
      call StrokePlane (time+dt, Insect2)
      call angular_velocities ( time+dt, Insect2 )

      ! fetch motion state at time-dt
      call FlappingMotion_right(time-dt, Insect1)
      call FlappingMotion_left (time-dt, Insect1)
      call StrokePlane (time-dt, Insect1)
      call angular_velocities (time-dt, Insect1)

      ! use centered finite differences
      Insect%rot_dt_l = (Insect2%rot_l - Insect1%rot_l)/(2.d0*dt)
      Insect%rot_dt_r = (Insect2%rot_r - Insect1%rot_r)/(2.d0*dt)
    else
      Insect1 = Insect
      Insect2 = Insect

      ! fetch motion state at time+dt
      call FlappingMotion_right(time+dt, Insect2)
      call FlappingMotion_left (time+dt, Insect2)
      call StrokePlane (time+dt, Insect2)
      call angular_velocities ( time+dt, Insect2 )

      ! fetch motion state at time-dt
      call FlappingMotion_right(time, Insect1)
      call FlappingMotion_left (time, Insect1)
      call StrokePlane (time, Insect1)
      call angular_velocities (time, Insect1)

      ! use centered finite differences
      Insect%rot_dt_l = (Insect2%rot_l - Insect1%rot_l)/dt
      Insect%rot_dt_r = (Insect2%rot_r - Insect1%rot_r)/dt
    endif
  end subroutine angular_accel



  subroutine insect_clean(Insect)
    use vars
    implicit none
    type(diptera),intent(inout)::Insect

    if (allocated(particle_points)) deallocate ( particle_points )
    call load_kine_clean( Insect%kine_wing_l )
    call load_kine_clean( Insect%kine_wing_r )
  end subroutine insect_clean


end module
