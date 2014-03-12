!-------------------------------------------------------------------------------
! 2nd generation insect module
!-------------------------------------------------------------------------------
! contains now motion of the body (translation + yaw,pitch,roll
! contains a stroke plane angle for direct implementation of literature results
!-------------------------------------------------------------------------------


! Main routine for drawing insects. Loops over the entire domain, computes
! coordinates in various systems (global-, body-, stroke-, wing-) and calls
! subroutines doing the actual job of defining the mask. Note all surfaces are
! smoothed.
subroutine Draw_Insect ( time )
  use fsi_vars
  use mpi
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr) :: x(1:3), x_body(1:3), x_wing_l(1:3), x_wing_r(1:3)
  real(kind=pr) :: x_eye_r(1:3), x_eye_l(1:3)
  real(kind=pr) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt
  real(kind=pr) :: xc_body(1:3), alpha_l, phi_l, phi_r, alpha_r, eta_l, eta_r
  real(kind=pr) :: alpha_dt_l, alpha_dt_r, phi_dt_l, phi_dt_r
  real(kind=pr) :: theta_dt_l, theta_dt_r, eta_stroke, theta_r, theta_l
  real(kind=pr), dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, &
  M1,M2,M3, M_stroke_l, M_stroke_r
  real(kind=pr), dimension(1:3)::rot_l, rot_r,rot_body,xc_head,xc_eye_l,&
  xc_eye_r, xc_pivot_r,xc_pivot_l, x_head, vc_body, v_tmp
  integer :: ix, iy, iz
  
  Insect%safety = 2.d0*dz
  Insect%smooth = 1.d0*dz
  
  ! some checks
  if ((mpirank==0).and.((iMoving.ne.1).or.(iPenalization.ne.1))) then
    write (*,*) "insects.f90::DrawInsect: the parameters iMoving or iPenalization are wrong."
    stop
  endif
  
  !------------------------------------
  ! this is the relative coordinates (in the body system)
  ! of some interesting points on the insect (head, eyes, hinges)
  !------------------------------------
  xc_head = Insect%x_head
  xc_eye_l = Insect%x_eye_l
  xc_eye_r = Insect%x_eye_r
  xc_pivot_r = Insect%x_pivot_r
  xc_pivot_l = Insect%x_pivot_l
  
  
  !--------------------------
  ! fetch current motion state
  !--------------------------
  call BodyMotion ( time, psi, beta, gamma, psi_dt, beta_dt, gamma_dt, xc_body, vc_body )
  call FlappingMotion_right(time, phi_r, alpha_r, theta_r, phi_dt_r, alpha_dt_r, theta_dt_r )
  call FlappingMotion_left (time, phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l )  
  call StrokePlane (time, eta_stroke)

  !-------------------------------------------------------
  ! write kinematics to disk (Dmitry, 28 Oct 2013)
  !-------------------------------------------------------     
  if(mpirank == 0) then
    open(17,file='kinematics.t',status='unknown',position='append')
    write (17,'(14(e12.5,1x))') time, xc_body, psi, beta, gamma, eta_stroke, &
    alpha_l, phi_l, theta_l, alpha_r, phi_r, theta_r
    close(17)
  endif

  !-------------------------------
  ! define the rotation matrices to change between coordinate systems
  !-------------------------------
  call Rx(M1,psi)
  call Ry(M2,beta)
  call Rz(M3,gamma)
  M_body = matmul(M1,matmul(M2,M3))
  
  call Ry(M1,eta_stroke)
  M_stroke_l = M1
  
  call Rx(M1,pi)
  call Ry(M2,eta_stroke)  
  M_stroke_r = matmul(M1,M2)

  call Ry(M1,alpha_l)
  call Rz(M2,theta_l)   ! Order changed (Dmitry, 7 Nov 2013)
  call Rx(M3,phi_l)
  M_wing_l = matmul(M1,matmul(M2,matmul(M3,M_stroke_l)))

  ! note the coordinate system is rotated so we don't need to inverse the sign
  ! of theta, and the wings still rotate in opposite direction
  call Ry(M1,-alpha_r)
  call Rz(M2,theta_r)   ! Order changed (Dmitry, 7 Nov 2013) 
  call Rx(M3,-phi_r)
  M_wing_r = matmul(M1,matmul(M2,matmul(M3,M_stroke_r)))

  !------------------------------------
  ! angular velocity vectors
  !------------------------------------
  rot_l = (/ phi_dt_l, alpha_dt_l, theta_dt_l /) ! in the wing reference frame
  rot_r = (/-phi_dt_r,-alpha_dt_r, theta_dt_r /) ! no need to inverse theta_dt sign
  rot_body = (/psi_dt, beta_dt, gamma_dt /)

  
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !--------------------------------------
           ! define the various coordinate systems
           ! we are going to use
           !--------------------------------------
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           x_head   = x_body - xc_head
           x_eye_l  = x_body - xc_eye_l
           x_eye_r  = x_body - xc_eye_r
           x_wing_l = matmul(M_wing_l,x_body-xc_pivot_l)
           x_wing_r = matmul(M_wing_r,x_body-xc_pivot_r)
           
           !--------------------------------------
           ! call body subroutines
           !--------------------------------------
           call DrawBody(ix,iy,iz,x_body)
           call DrawHead(ix,iy,iz,x_head)
           call DrawEye(ix,iy,iz,x_eye_r)
           call DrawEye(ix,iy,iz,x_eye_l)
           
           !--------------------------------------
           ! wings
           !--------------------------------------
           call DrawWing(ix,iy,iz,x_wing_l,M_wing_l,rot_l)
           call DrawWing(ix,iy,iz,x_wing_r,M_wing_r,rot_r)
           
           !--------------------------------------
           ! add solid body rotation in the body-reference frame
           !--------------------------------------
           if (mask(ix,iy,iz) > 0.d0) then
            ! add solid body rotation to the translational velocity field
            v_tmp(1) = vc_body(1)+rot_body(2)*x_body(3)-rot_body(3)*x_body(2)
            v_tmp(2) = vc_body(2)+rot_body(3)*x_body(1)-rot_body(1)*x_body(3)
            v_tmp(3) = vc_body(3)+rot_body(1)*x_body(2)-rot_body(2)*x_body(1)
            
            us(ix,iy,iz,1:3)=matmul(transpose(M_body),us(ix,iy,iz,1:3)+v_tmp)
           endif
        enddo
    enddo
  enddo
end subroutine Draw_Insect


!-------------------------------------------------------
! short for the smooth step function.
! the smooting is defined in Insect%smooth, here we need only x, and the 
! thickness (i.e., in the limit, steps=1 if x<t and steps=0 if x>t
!-------------------------------------------------------
real(kind=pr) function steps(x,t)
  use fsi_vars
  use mpi
  implicit none
  real(kind=pr) :: f,x,t
  call smoothstep(f,x,t,Insect%smooth)
  steps=f
end function


!-------------------------------------------------------
! Compute angle from coefficients provided by Maeda
!-------------------------------------------------------
subroutine get_dangle( angles, F, a, b, shift_phase, initial_phase, dangle, dangle_dt )
  use fsi_vars
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



! Insect free flight dynamics.
! RHS of the ODE system.
subroutine dynamics_insect(time,it)
  use mpi
  use fsi_vars
  implicit none

  integer, intent(in) :: it
  integer :: ilegs  ! legs force: 0=off; 1=on
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: accx,accz,mass_solid,mass_fluid,gravity
  real (kind=pr) :: time_ref,length_ref,density_ref
  real (kind=pr) :: anglegs,kzlegsmax,dzlegsmax,t0,tmaxlegs,kzlegs0,anglegsend
  real (kind=pr) :: kxlegs,kzlegs,fxlegs,fzlegs,fxaero,fzaero,displacement_z

  ! reference values
  time_ref = 3.664845d-3 ! in s
  length_ref = 2.39d-3   ! in m
  density_ref = 1.225d0  ! in kg/m3

  ! mass
  mass_solid = 0.91d-6 / (density_ref*length_ref**3)
!  mass_solid = 2.0*0.91d-6 / (density_ref*length_ref**3)  ! heavy, x2
  mass_fluid = 0.0d0  ! TODO

  ! gravity acceleration
  gravity = -9.81d0*time_ref*time_ref/length_ref

  select case (Insect%BodyMotion)
  case ("takeoff")
    if (Insect%KineFromFile=="simplified_dynamic") then

      ! Current body center displacement
      displacement_z = SolidDyn%var_new(2)

      ! Legs model parameters
      ilegs = 1
!       anglegsend = 0.25d0*pi ! case 1 and 3
      anglegsend = 0.5d0*pi ! case2
!      kzlegsmax = 0.057d0 / (density_ref*length_ref**3/time_ref**2) ! case1
      kzlegsmax = 0.024d0 / (density_ref*length_ref**3/time_ref**2) ! case2
!      kzlegsmax = 0.08d0 / (density_ref*length_ref**3/time_ref**2) ! case3
      dzlegsmax = 0.00065d0 / length_ref
      t0 = 0.0005d0 / time_ref
      tmaxlegs = t0 + 0.0013d0 / time_ref
      kzlegs0 = (mass_solid-mass_fluid)*(-gravity) / dzlegsmax  

      ! Legs force
      fxlegs = 0.0d0
      fzlegs = 0.0d0
      if (ilegs>0) then
        ! If legs touch the ground
        if (displacement_z<dzlegsmax) then
          ! Compute torsion spring stiffness and angle
          if (time<t0) then
            kzlegs = kzlegs0
            anglegs = 0.5d0*pi
          elseif (time<tmaxlegs) then
            kzlegs = kzlegs0 + (kzlegsmax-kzlegs0) * (time-t0)/(tmaxlegs-t0)
            anglegs = 0.5d0*pi + (anglegsend-0.5d0*pi) * (time-t0)/(tmaxlegs-t0)
          else
            kzlegs = kzlegsmax
            anglegs = anglegsend
          endif
          if (tan(anglegs)<1.0d-8) then
            kxlegs = 0.0d0
          else
            kxlegs = kzlegs/tan(anglegs)
          endif
          ! Compute legs forces
          fxlegs = kxlegs*(dzlegsmax-displacement_z)
          fzlegs = kzlegs*(dzlegsmax-displacement_z)
        endif
      endif

      ! Fluid force. Unsteady correction treated implicitly
      fxaero = GlobalIntegrals%Force(1)
      fzaero = GlobalIntegrals%Force(3)

      ! Accelerations
      accx = (fxaero+fxlegs)/(mass_solid-mass_fluid)
      accz = (fzaero+fzlegs)/(mass_solid-mass_fluid) + gravity 

      ! RHS of the solid ODEs
      SolidDyn%rhs_this(1) = SolidDyn%var_this(3) ! horizontal velocity
      SolidDyn%rhs_this(2) = SolidDyn%var_this(4) ! vertical velocity
      SolidDyn%rhs_this(3) = accx ! horizontal acceleration
      SolidDyn%rhs_this(4) = accz ! vertical acceleration

      ! Write legs forces to file
      if(mpirank == 0) then
        open(14,file='legs.t',status='unknown',position='append')
        write (14,'(3(e12.5,1x))') time,fxlegs,fzlegs
        close(14)
      endif

    endif
  case default
    if (mpirank==0) then
    write (*,*) &
    "insects.f90::dynamics_insects case not defined"
    stop
    endif
  end select

end subroutine dynamics_insect


! Initialize insect free flight dynamics solver
subroutine dynamics_insect_init(idynamics)
  use mpi
  use fsi_vars
  implicit none

  integer, intent(out) :: idynamics

  ! dynamics solver inactive by default
  idynamics = 0

  select case (Insect%BodyMotion)
  case ("takeoff")
    if (Insect%KineFromFile=="simplified_dynamic") then
      idynamics = 1
      SolidDyn%var_new(1) = 0.0d0
      SolidDyn%var_new(2) = 0.0d0
      SolidDyn%var_new(3) = 0.0d0
      SolidDyn%var_new(4) = 0.0d0
    endif
  end select

end subroutine dynamics_insect_init



