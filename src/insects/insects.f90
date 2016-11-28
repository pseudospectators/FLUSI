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
  logical, save :: body_already_drawn = .false.
  character(len=strlen), save :: body_moves="yes"
  ! arrays for fourier coefficients are fixed size (avoiding issues with allocatable
  ! elements in derived datatypes) this is their length:
  integer, parameter :: nfft_max = 1024
  ! Maximum number of Hermite interpolation nodes (hardcoded because of sxf90 compiler requirements)
  integer, parameter :: nhrmt_max = 10000
  ! only used in steps (to reduce number of arguments)
  real(kind=pr), save :: smoothing

  include "type_definitions.f90"

contains


  !---------------------------------------
  ! note these include files also have to be specified as dependencies in the
  ! Makefile for make to check if one of them changed
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

    real(kind=pr) :: t1
    real(kind=pr),dimension(1:3) :: x, x_body, v_tmp
    real(kind=pr),dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, M_body_inv
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

    ! some checks
    if ((mpirank==0).and.((iMoving.ne.1).or.(iPenalization.ne.1))) then
      write (*,*) "insects.f90::DrawInsect: the parameters iMoving or iPenalization are wrong."
      call abort(4453)
    endif


    if (nx/=1) then
      Insect%smooth = 1.0d0*max(dz,dy,dx)
    else
      Insect%smooth = 1.0d0*max(dz,dy)
    endif
    Insect%safety = 3.5d0*Insect%smooth
    smoothing = Insect%smooth


    ! delete old mask
    call delete_old_mask( time, mask, mask_color, us, Insect )
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
    call body_rotation_matrix( Insect, M_body )
    call wing_right_rotation_matrix( Insect, M_wing_r )
    call wing_left_rotation_matrix( Insect, M_wing_l )

    ! inverse of the body rotation matrices
    M_body_inv = transpose(M_body)

    ! body angular velocity vector in b/g coordinate system
    call body_angular_velocity( Insect, Insect%rot_body_b, Insect%rot_body_g, M_body )
    ! rel+abs wing angular velocities in the w/b/g coordinate system
    call wing_angular_velocities ( time, Insect, M_body )
    ! angular acceleration for wings (required for inertial power)
    call wing_angular_accel( time, Insect )

    !-----------------------------------------------------------------------------
    ! vector from body centre to left/right pivot point in global reference frame,
    ! for aerodynamic power
    !-----------------------------------------------------------------------------
    Insect%x_pivot_l_g = matmul(M_body_inv, Insect%x_pivot_l)
    Insect%x_pivot_r_g = matmul(M_body_inv, Insect%x_pivot_r)

    !-----------------------------------------------------------------------------
    ! write kinematics to disk (Dmitry, 28 Oct 2013)
    ! do so only every itdrag time steps (Thomas, 8 Jul 2014)
    !-----------------------------------------------------------------------------
    counter = counter + 1
    if ((mpirank == 0).and.(mod(counter,itdrag)==0)) then
      open  (17,file='kinematics.t',status='unknown',position='append')
      write (17,'(26(es15.8,1x))') time, Insect%xc_body_g, Insect%psi, Insect%beta, &
      Insect%gamma, Insect%eta_stroke, Insect%alpha_l, Insect%phi_l, &
      Insect%theta_l, Insect%alpha_r, Insect%phi_r, Insect%theta_r, &
      Insect%rot_rel_wing_l_w, Insect%rot_rel_wing_r_w, &
      Insect%rot_dt_wing_l_w, Insect%rot_dt_wing_r_w
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
      Insect%x_pivot_r,Insect%rot_rel_wing_r_w )
    endif

    if (Insect%LeftWing == "yes") then
      call draw_wing(mask,mask_color,us,Insect,color_l,M_body,M_wing_l,&
      Insect%x_pivot_l,Insect%rot_rel_wing_l_w )
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
          x = periodize_coordinate(x - Insect%xc_body_g)
          x_body = matmul(M_body,x)
          ! add solid body rotation in the body-reference frame, if color
          ! indicates that this part of the mask belongs to the insect
          if ((mask(ix,iy,iz) > 0.d0).and.(mask_color(ix,iy,iz)>0)) then

            ! translational part. we compute the rotational part in the body
            ! reference frame, therefore, we must transform the body translation
            ! velocity Insect%vc (which is in global coordinates) to the body frame
            v_tmp = matmul(M_body,Insect%vc_body_g)

            ! add solid body rotation to the translational velocity field. Note
            ! that rot_body_b and x_body are in the body reference frame
            v_tmp(1) = v_tmp(1)+Insect%rot_body_b(2)*x_body(3)-Insect%rot_body_b(3)*x_body(2)
            v_tmp(2) = v_tmp(2)+Insect%rot_body_b(3)*x_body(1)-Insect%rot_body_b(1)*x_body(3)
            v_tmp(3) = v_tmp(3)+Insect%rot_body_b(1)*x_body(2)-Insect%rot_body_b(2)*x_body(1)

            ! the body motion is added to the wing motion, which is already in us
            ! and they are also in the body refrence frame. However, us has to be
            ! in the global reference frame, so M_body_inverse is applied
            us(ix,iy,iz,1:3) = matmul( M_body_inv, us(ix,iy,iz,1:3)+v_tmp )
          endif
        enddo
      enddo
    enddo
    time_insect_vel = time_insect_vel + MPI_wtime() - t1

    ! print some important numbers, routine exectutes only once during a simulation
    call print_insect_reynolds_numbers( Insect )

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

    !-----------
    ! left wing
    !-----------
    ! relative angular velocity, in global system
    omrel = Insect%rot_rel_wing_l_g

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_l)%Torque + &
             Insect%PartIntegrals(color_l)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_l)%APow = - sum( momrel * omrel )

    !-----------
    ! right wing
    !-----------
    ! relative angular velocity, in global system
    omrel = Insect%rot_rel_wing_r_g

    ! compute moment with respect to the pivot point
    ! initialize it as the moment with respect to insect's centre point
    momrel = Insect%PartIntegrals(color_r)%Torque + &
             Insect%PartIntegrals(color_r)%Torque_unst

    ! aerodynamic power
    Insect%PartIntegrals(color_r)%APow = - sum( momrel * omrel )

    !-----------
    ! Total aerodynamic power
    !-----------
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
  !       Insect%rot_dt_wing_l_w (global): left wing angular acceleration
  !       Insect%rot_dt_wing_r_w (global): right wing angular acceleration
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
    a(1) = Insect%Jxx * Insect%rot_dt_wing_l_w(1) + Insect%Jxy * Insect%rot_dt_wing_l_w(2)
    a(2) = Insect%Jxy * Insect%rot_dt_wing_l_w(1) + Insect%Jyy * Insect%rot_dt_wing_l_w(2)
    a(3) = Insect%Jzz * Insect%rot_dt_wing_l_w(3)

    b(1) = Insect%Jxx * Insect%rot_rel_wing_l_w(1) + Insect%Jxy * Insect%rot_rel_wing_l_w(2)
    b(2) = Insect%Jxy * Insect%rot_rel_wing_l_w(1) + Insect%Jyy * Insect%rot_rel_wing_l_w(2)
    b(3) = Insect%Jzz * Insect%rot_rel_wing_l_w(3)

    Insect%PartIntegrals(color_l)%IPow = &
    Insect%rot_rel_wing_l_w(1) * (a(1)+Insect%rot_rel_wing_l_w(2)*b(3)-Insect%rot_rel_wing_l_w(3)*b(2)) +&
    Insect%rot_rel_wing_l_w(2) * (a(2)+Insect%rot_rel_wing_l_w(3)*b(1)-Insect%rot_rel_wing_l_w(1)*b(3)) +&
    Insect%rot_rel_wing_l_w(3) * (a(3)+Insect%rot_rel_wing_l_w(1)*b(2)-Insect%rot_rel_wing_l_w(2)*b(1))

    !-- RIGHT WING
    a(1) = Insect%Jxx * Insect%rot_dt_wing_r_w(1) + Insect%Jxy * Insect%rot_dt_wing_r_w(2)
    a(2) = Insect%Jxy * Insect%rot_dt_wing_r_w(1) + Insect%Jyy * Insect%rot_dt_wing_r_w(2)
    a(3) = Insect%Jzz * Insect%rot_dt_wing_r_w(3)

    b(1) = Insect%Jxx * Insect%rot_rel_wing_r_w(1) + Insect%Jxy * Insect%rot_rel_wing_r_w(2)
    b(2) = Insect%Jxy * Insect%rot_rel_wing_r_w(1) + Insect%Jyy * Insect%rot_rel_wing_r_w(2)
    b(3) = Insect%Jzz * Insect%rot_rel_wing_r_w(3)

    Insect%PartIntegrals(color_r)%IPow = &
    Insect%rot_rel_wing_r_w(1) * (a(1)+Insect%rot_rel_wing_r_w(2)*b(3)-Insect%rot_rel_wing_r_w(3)*b(2)) +&
    Insect%rot_rel_wing_r_w(2) * (a(2)+Insect%rot_rel_wing_r_w(3)*b(1)-Insect%rot_rel_wing_r_w(1)*b(3)) +&
    Insect%rot_rel_wing_r_w(3) * (a(3)+Insect%rot_rel_wing_r_w(1)*b(2)-Insect%rot_rel_wing_r_w(2)*b(1))

    ipowtotal = Insect%PartIntegrals(color_r)%IPow + Insect%PartIntegrals(color_l)%IPow

  end subroutine inert_power

  !-----------------------------------------------------------------------------
  ! Body angular velocity vector
  !-----------------------------------------------------------------------------
  ! Variant (a) : free flight with quaternion solver
  !
  !    when using the quaternion based free-flight solver, the angular
  !    velocity of the body is computed dynamically, and the rotation matrix that
  !    brings us from global to body system is computed with quaternions. In body_motion
  !    the solver sets Insect%rot_body_b, so here we compute only rot_body_g
  !
  ! Variant (b) : imposed (prescribed) body dynamics
  !
  !    If the free flight solver is not active, body yaw,pitch,roll are known, so
  !    given yaw.pitch roll angles and their time derivatives, return the bodies
  !    angular velocity vector in global and body frame
  !-----------------------------------------------------------------------------
  subroutine body_angular_velocity( Insect, rot_body_b, rot_body_g, M_body )
    implicit none

    type(diptera), intent(inout) :: Insect
    real(kind=pr), intent(in) :: M_body(1:3,1:3)
    real(kind=pr), dimension(1:3), intent(out) :: rot_body_b, rot_body_g
    real(kind=pr) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt

    psi = Insect%psi
    beta = Insect%beta
    gamma = Insect%gamma
    psi_dt = Insect%psi_dt
    beta_dt = Insect%beta_dt
    gamma_dt = Insect%gamma_dt

    if ( Insect%BodyMotion == "free_flight" ) then
      ! variant (a)
      rot_body_b = Insect%rot_body_b ! copy (useless, actually, but required for interface)
      rot_body_g = matmul( transpose(M_body), rot_body_b)
    else
      ! variant (b)
      ! in global frame
      rot_body_g = (/ psi_dt*cos(beta)*cos(gamma)-beta_dt*sin(gamma) ,&
                      beta_dt*cos(gamma)+psi_dt*cos(beta)*sin(gamma) ,&
                      gamma_dt-psi_dt*sin(beta) /)
      ! in body frame
      rot_body_b = (/ psi_dt-gamma_dt*sin(beta) ,&
                      beta_dt*cos(psi)+gamma_dt*cos(beta)*sin(psi) ,&
                      gamma_dt*cos(beta)*cos(psi)-beta_dt*sin(psi) /)
    endif
  end subroutine body_angular_velocity



  !-------------------------------------------------------------------------------
  ! given the angles of each wing (and their time derivatives), compute
  ! the angular velocity vectors for both wings.
  ! output:
  !    Insect%rot_rel_wing_r_w    relative angular velocity of right wing (wing frame)
  !    Insect%rot_rel_wing_r_b    relative angular velocity of right wing (body frame)
  !    Insect%rot_rel_wing_r_g    relative angular velocity of right wing (glob frame)
  !    Insect%rot_abs_wing_r_g    absolute angular velocity of right wing (glob frame)
  !-------------------------------------------------------------------------------
  subroutine wing_angular_velocities ( time, Insect, M_body )
    use vars
    implicit none

    real(kind=pr), intent(in) :: time
    real(kind=pr), intent(in) :: M_body(1:3,1:3)
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
    ! angular velocity vectors (in wing system)
    !-----------------------------------------------------------------------------
    rot_l_alpha = (/ 0.0d0, alpha_dt_l, 0.0d0 /)
    rot_l_theta = (/ 0.0d0, 0.0d0, theta_dt_l /)
    rot_l_phi   = (/ phi_dt_l, 0.0d0, 0.0d0   /)
    rot_r_alpha = (/ 0.0d0, -alpha_dt_r, 0.0d0/)
    rot_r_theta = (/ 0.0d0, 0.0d0, theta_dt_r /)
    rot_r_phi   = (/ -phi_dt_r, 0.0d0, 0.0d0  /)

    ! in the wing coordinate system
    Insect%rot_rel_wing_l_w = matmul(M_wing_l,matmul(transpose(M_stroke_l),matmul(transpose(M3_l), &
    rot_l_phi+matmul(transpose(M2_l),rot_l_theta+matmul(transpose(M1_l), &
    rot_l_alpha)))))
    Insect%rot_rel_wing_r_w = matmul(M_wing_r,matmul(transpose(M_stroke_r),matmul(transpose(M3_r), &
    rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
    rot_r_alpha)))))

    ! direct definition, equivalent to what is above.
    ! Insect%rot_rel_wing_l_w = (/phi_dt_l*cos(alpha_l)*cos(theta_l)-theta_dt_l*sin(alpha_l),&
    !   alpha_dt_l-phi_dt_l*sin(theta_l),&
    !   theta_dt_l*cos(alpha_l)+phi_dt_l*sin(alpha_l)*cos(theta_l)/)

    ! prior to the call of this routine, the routine body_angular_velocity has
    ! computed the body angular velocity (both g/b frames) so here now we can also
    ! compute global and absolute wing angular velocities.
    Insect%rot_rel_wing_l_b = matmul( transpose(M_wing_l), Insect%rot_rel_wing_l_w )
    Insect%rot_rel_wing_r_b = matmul( transpose(M_wing_r), Insect%rot_rel_wing_r_w )

    Insect%rot_rel_wing_l_g = matmul( transpose(M_body), Insect%rot_rel_wing_l_b )
    Insect%rot_rel_wing_r_g = matmul( transpose(M_body), Insect%rot_rel_wing_r_b )

    Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
    Insect%rot_abs_wing_r_g = Insect%rot_body_g + Insect%rot_rel_wing_r_g


    if (Insect%wing_fsi == "yes") then
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! overwrite the left wing
      Insect%rot_rel_wing_l_w = Insect%STATE(18:20)
      Insect%rot_rel_wing_l_b = matmul( transpose(M_wing_l), Insect%rot_rel_wing_l_w )
      Insect%rot_rel_wing_l_g = matmul( transpose(M_body), Insect%rot_rel_wing_l_b )
      ! the last guy is actually unused, as we have non-rotating body
      Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
    endif

  end subroutine wing_angular_velocities



  !-------------------------------------------------------------------------------
  ! Numerically estimate (it's a very precise estimation) the angular acceleration
  ! vectors for both wings, using one-sided finite differences (in the future)
  ! NOTE: this routine requires us to be able to evaluate both body and wing state
  !       at arbitrary times.
  !-------------------------------------------------------------------------------
  subroutine wing_angular_accel( time, Insect )
    use vars
    implicit none
    real(kind=pr), intent(in) :: time
    type(diptera), intent(inout) :: Insect

    real(kind=pr) :: M_body(1:3,1:3), rot_dt_wing_g(1:3), M_wing_r(1:3,1:3), M_wing_l(1:3,1:3)
    type(diptera) :: Insect2
    real(kind=pr) :: dt,t

    dt = 1.0d-8
    Insect2 = Insect

    Insect%rot_dt_wing_l_w = 0.d0
    Insect%rot_dt_wing_r_w = 0.d0

    ! fetch motion state at time+dt
    call BodyMotion (time+dt, Insect2)
    call FlappingMotion_right(time+dt, Insect2)
    call FlappingMotion_left (time+dt, Insect2)
    call StrokePlane (time+dt, Insect2)
    call body_rotation_matrix( Insect2, M_body )
    call wing_angular_velocities ( time+dt, Insect2, M_body )
    ! this is the current state:
    call body_rotation_matrix( Insect, M_body )
    call wing_right_rotation_matrix( Insect, M_wing_r )
    call wing_left_rotation_matrix( Insect, M_wing_l )

    ! use one-sided finite differences to derive the absolute angular velocity with
    ! respect to time. NOte in older code versions, this was wrong, as we derived
    ! the ang. vel. in the wing coordinate system, which is a moving reference frame.

    ! now happily, the old results are still correct, as long as the body does not rotate
    ! see, e.g., this document https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjzl-XP6_LNAhWoC5oKHUdDCHwQFggoMAI&url=http%3A%2F%2Focw.mit.edu%2Fcourses%2Faeronautics-and-astronautics%2F16-07-dynamics-fall-2009%2Flecture-notes%2FMIT16_07F09_Lec08.pdf&usg=AFQjCNHzEB-n_NMm6K3J1eRpIaGnuKpW0Q&sig2=yEPNin3bL5DnWauNJk2hcw&bvm=bv.126993452,d.bGs&cad=rjt
    ! however, if the body moves, an additional term occurs, and this was indeed missing
    ! in previous results.
    rot_dt_wing_g = (Insect2%rot_rel_wing_l_g - Insect%rot_rel_wing_l_g) / dt
    Insect%rot_dt_wing_l_w = matmul(M_wing_l,matmul(M_body, rot_dt_wing_g))

    rot_dt_wing_g = (Insect2%rot_rel_wing_r_g - Insect%rot_rel_wing_r_g) / dt
    Insect%rot_dt_wing_r_w = matmul(M_wing_r,matmul(M_body, rot_dt_wing_g))


    ! if (root) then
    !   write(*,*) "L new code", Insect%rot_dt_wing_l_w
    !   write(*,*) "L old code", (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt
    !   write(*,*) "rot_rel_wing_l_g", Insect%rot_rel_wing_l_g
    !   write(*,*) "rot_body_g", Insect%rot_body_g
    !   write(*,*) "rot_body_b", Insect%rot_body_b
    !   write(*,*) "rot_dt_body_g", (Insect2%rot_body_g-Insect%rot_body_g)/dt
    ! endif
    !
    ! if (root) then
    !   t = 0.d0
    !   open  (17,file='test.t',status='replace')
    !   do while (t<=6.d0)
    !     call FlappingMotion_left ( t, Insect)
    !     write (17,'(7(es15.8,1x))') t,  &
    !     Insect%alpha_l, Insect%phi_l, Insect%theta_l, &
    !     Insect%alpha_dt_l, Insect%phi_dt_l, Insect%theta_dt_l
    !
    !     t = t + 1.0d-5
    !   end do
    !   close (17)
    !   call abort(7)
    ! endif

    ! ! this is the OLD CODE. It is actually wrong, since it computes the time derivative
    ! ! in the moving refrence frame of the wing, which is not the actual angular acceleration
    ! Insect%rot_dt_wing_r_w = (Insect2%rot_rel_wing_r_w - Insect%rot_rel_wing_r_w)/dt
    ! Insect%rot_dt_wing_l_w = (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt


    if (Insect%wing_fsi == "yes") then
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! overwrite the left wings acceleration. the right hand side 18,19,20 is the
      ! time derivative of the angular velocity, so the acceleration is readily available
      Insect%rot_dt_wing_l_w = Insect%RHS_THIS(18:20)
    endif

  end subroutine wing_angular_accel


  subroutine delete_old_mask( time, mask, mask_color, us, Insect )
    use vars
    implicit none

    real(kind=pr), intent(in) :: time
    type(diptera),intent(inout) :: Insect
    real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    integer(kind=2) :: color_body, color_l, color_r

    ! colors for Diptera (one body, two wings)
    color_body = 1
    color_l = 2
    color_r = 3

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
      ! (and any other) call, and therefore we can safely reset the entire mask to zeros.
      mask = 0.d0
      mask_color = 0
      us = 0.d0
    endif

  end subroutine delete_old_mask


  !-----------------------------------------------------------------------------
  ! return the body rotation matrix
  !-----------------------------------------------------------------------------
  subroutine body_rotation_matrix( Insect, M_body )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=pr),intent(out) :: M_body(1:3,1:3)
    real(kind=pr), dimension(1:3,1:3) :: M1_b, M2_b, M3_b

    if (Insect%BodyMotion=="free_flight") then
      ! entries 7,8,9,10 of the Insect%STATE vector are the body quaternion
      call rotation_matrix_from_quaternion( Insect%STATE(7:10), M_body)
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
  end subroutine body_rotation_matrix

  !-----------------------------------------------------------------------------
  ! return the rotation matrix for the right wing
  !-----------------------------------------------------------------------------
  subroutine wing_right_rotation_matrix( Insect, M_wing_r )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=pr),intent(out) :: M_wing_r(1:3,1:3)
    real(kind=pr), dimension(1:3,1:3) :: M1, M2, M3, M_stroke_r


    call Rx(M1,pi)
    call Ry(M2,Insect%eta_stroke)
    M_stroke_r = matmul(M1,M2)

    ! note the coordinate system is rotated so we don't need to inverse the sign
    ! of theta, and the wings still rotate in opposite direction
    call Ry(M1,-Insect%alpha_r)
    call Rz(M2, Insect%theta_r)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,-Insect%phi_r)
    M_wing_r = matmul(M1,matmul(M2,matmul(M3,M_stroke_r)))
  end subroutine wing_right_rotation_matrix


  !-----------------------------------------------------------------------------
  ! return the rotation matrix for the left wing
  !-----------------------------------------------------------------------------
  subroutine wing_left_rotation_matrix( Insect, M_wing_l )
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=pr),intent(out) :: M_wing_l(1:3,1:3)
    real(kind=pr),dimension(1:3,1:3) :: M1, M2, M3, M_stroke_l

    if ( Insect%wing_fsi /= "yes" ) then
      ! we're not using the wing fsi solver, so the wings follow a prescribed
      ! motion and we can compute the rotation matrix from the angles
      call Ry(M1,Insect%eta_stroke)
      M_stroke_l = M1

      call Ry(M1,Insect%alpha_l)
      call Rz(M2,Insect%theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3,Insect%phi_l)
      M_wing_l = matmul(M1,matmul(M2,matmul(M3,M_stroke_l)))
    else
      !**********************************
      !** Wing fsi model               **
      !**********************************
      ! in the wing FSI case, a quaternion-based formulation is used to get the
      ! rotation matrix from the wing quaterion. note the wing quaternion are the
      ! entries 14,15,16,17 of the Insect%STATE vector
      ! entries 18,19,20 are the angular VELOCITY of the wing
      call rotation_matrix_from_quaternion( Insect%STATE(14:17), M_wing_l)
    endif
  end subroutine wing_left_rotation_matrix

  !-----------------------------------------------------------------------------
  ! Compute and print a couple of important numbers for insects
  !-----------------------------------------------------------------------------
  subroutine  print_insect_reynolds_numbers( Insect )
    implicit none
    type(diptera),intent(inout) :: Insect
    type(diptera) :: Insect_copy
    real(kind=pr) :: area, Re_f, Re
    real(kind=pr) :: time, dt
    real(kind=pr) :: phil_min, phil_max, phir_min, phir_max
    logical, save :: first_call = .true.

    ! the second call is just a return statement
    if ( first_call .eqv. .false.) return

    ! only root does this...
    if (root) then
      ! we need the wing area to compute the mean wing chord
      call compute_wing_surface(Insect, area)
      write(*,'(50("~"))')
      write(*,'("Wing area is A=",g15.8)') area
      write(*,'("Mean chord length is c_m=",g15.8)') area/1.d0 ! note c_m = A/R but R=1

      if (Insect%wing_fsi /= 'yes') then
        ! first we computethe stroke amplitude of the positional angle phi (for
        ! both wings). for safety, we make a copy of the insect, since the routines
        ! for the flapping motion write to this object, and we want to prevent any
        ! unwanted side effects
        Insect_copy = Insect
        time = 0.d0
        dt = 1.0d-3
        phil_min = 0.d0
        phil_max = 0.d0
        phir_min = 0.d0
        phir_max = 0.d0
        ! we use only one stroke ( the first one )
        do while (time < 1.d0)
          call FlappingMotion_left ( time, Insect_copy )
          call FlappingMotion_right ( time, Insect_copy )
          phil_min = min( phil_min, Insect_copy%phi_l )
          phil_max = max( phil_max, Insect_copy%phi_l )
          phir_min = min( phir_min, Insect_copy%phi_r )
          phir_max = max( phir_max, Insect_copy%phi_r )
          time = time + dt
        end do
        write(*,'("All following quantities are based on the first stroke 0.0 <= t <= 1.0")')
        write(*,'("Stroke amplitude is PHI_L=",g15.8)') abs(phil_max) + abs(phil_min)
        write(*,'("Stroke amplitude is PHI_R=",g15.8)') abs(phir_max) + abs(phir_min)
        write(*,'("Re_left  = 2*phi*R*f*c_m / nu =",g15.8)') 2.d0*(abs(phil_max)+abs(phil_min))*1.d0*1.d0*area/nu
        write(*,'("Re_right = 2*phi*R*f*c_m / nu =",g15.8)') 2.d0*(abs(phir_max)+abs(phir_min))*1.d0*1.d0*area/nu
        write(*,'("Re_f = R*R*f / nu =",g15.8)') 1.d0/nu
      else
        write(*,*) "In the case of wing_fsi problems, we do not know the wing"
        write(*,*) "kinematics before the computation. Therefore, we cannot tell"
        write(*,*) "the Reynolds number in advance."
        write(*,'("Re_f = R*R*f / nu =",g15.8)') 1.d0/nu
      endif
      write(*,'(50("~"))')
    endif

    first_call = .false.
  end subroutine print_insect_reynolds_numbers


  subroutine insect_clean(Insect)
    use vars
    implicit none
    type(diptera),intent(inout)::Insect

    if (allocated(particle_points)) deallocate ( particle_points )
    call load_kine_clean( Insect%kine_wing_l )
    call load_kine_clean( Insect%kine_wing_r )
  end subroutine insect_clean


  !-----------------------------------------------------------------------------
  ! In some cases, we need to reconstruct the mask or the body system in postprocessing
  ! if the free_flight solver was used, the body system cannnot be simply evaluated from closed-from
  ! expressions.
  ! In these case, we have to read the rigidsolidsolver.t file, which contains the Insect%STATE
  ! and from this the body orientation, rotation matrix, etc can be computed. So here we read this file
  ! and return Insect%STATE at the desired time (linear interpolation is used)
  !-----------------------------------------------------------------------------
  subroutine read_insect_STATE_from_file(time, Insect)
    use vars
    use ini_files_parser_mpi
    implicit none
    real(kind=pr), intent(in) :: time
    type(diptera),intent(inout) :: Insect
    integer :: num_lines, n_header = 1, i
    character(len=maxcolumns) :: dummy
    real(kind=pr), allocatable :: data(:,:)

    ! read rigidsolidsolver.t file
    ! skip header, count lines, read
    call count_lines_in_ascii_file_mpi('rigidsolidsolver.t', num_lines, n_header)
    ! read contents of file
    allocate( data(1:num_lines,1:14))
    call read_array_from_ascii_file_mpi('rigidsolidsolver.t', data , n_header)

    ! interpolate in time
    i=1
    do while (data(i,1) <= time .and. i<num_lines)
      i=i+1
    enddo
    ! we now have data(i-1,1) <= time < data(i,1)
    ! use linear interpolation
    Insect%STATE = 0.d0
    Insect%STATE(1:13) = data(i-1,2:14) + (time - data(i-1,1)) * (data(i,2:14)-data(i-1,2:14)) / (data(i,1)-data(i-1,1))

    if (root) write(*,*) "The extracted Insect%STATE vector is:"
    if (root) write(*,'(21(es12.4,1x))') time,Insect%STATE
    deallocate (data)
  end subroutine

end module
