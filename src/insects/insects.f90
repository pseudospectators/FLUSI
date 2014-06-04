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
  real(kind=pr) :: xc_body(1:3), alpha_l, phi_l, phi_r, alpha_r
  real(kind=pr) :: alpha_dt_l, alpha_dt_r, phi_dt_l, phi_dt_r, t1
  real(kind=pr) :: theta_dt_l, theta_dt_r, eta_stroke, theta_r, theta_l
  real(kind=pr), dimension(1:3,1:3) :: M_body, M_wing_l, M_wing_r, &
  M1, M2, M3, M_stroke_l, M_stroke_r, M_body_inv, M_wing_l_inv, M_wing_r_inv  
  real(kind=pr), dimension(1:3)::rot_l, rot_r,rot_body,xc_head,xc_eye_l,&
  xc_eye_r, xc_pivot_r,xc_pivot_l, x_head, vc_body, v_tmp
  integer :: ix, iy, iz
  integer :: color_body, color_l, color_r
  ! tell the code what type of subroutine to call for the wings: fourier or simple
  logical :: fourier_wing = .true. ! almost always we have this
  logical :: HasEye = .false., HasHead=.true.
  
  !-- decide what wing routine to call (call simplified wing 
  ! routines that don't use fourier)
  if (Insect%WingShape=='TwoEllipses') fourier_wing = .false.
  if (Insect%WingShape=='rectangular') fourier_wing = .false.
  !-- decide if to call draw_eye (logicals are much much faster)
  if (Insect%HasEye=="yes")  HasEye = .true.
  !-- decide if to call draw_head (logicals are much much faster)
  if (Insect%HasHead=="yes") HasHead=.true.
    
  !-- define the wings fourier coeffients, but only once  
  if (fourier_wing) call Setup_Wing_Fourier_coefficients()  
    
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

  Insect%vc_body = vc_body   ! This is required for aerodynamic power

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

  !-------------------------------------------------------
  ! inverse of the rotation matrices
  !-------------------------------------------------------     
  M_body_inv = transpose(M_body)
  M_wing_l_inv = transpose(M_wing_l)
  M_wing_r_inv = transpose(M_wing_r)

  !-------------------------------------------------------
  ! angular velocity in the global reference frame
  !-------------------------------------------------------
  Insect%rot_body_glob = matmul(M_body_inv, rot_body)
  Insect%rot_l_glob = matmul(M_body_inv, matmul(M_wing_l_inv, rot_l) + rot_body)
  Insect%rot_r_glob = matmul(M_body_inv, matmul(M_wing_r_inv, rot_r) + rot_body) 
 
  !-------------------------------------------------------
  ! vector from body centre to left/right pivot point in global reference frame, 
  ! for aerodynamic power
  !-------------------------------------------------------
  Insect%x_pivot_l_glob = matmul(M_body_inv, xc_pivot_l)
  Insect%x_pivot_r_glob = matmul(M_body_inv, xc_pivot_r) 

  !-------------------------------------------------------
  ! colors for Diptera (one body, two wings)
  !-------------------------------------------------------  
  color_body = 1
  color_l = 2
  color_r = 3

  !-------------------------------------------------------
  ! Draw indivudual parts of the Diptera. Separate loops are faster
  ! since the compiler can optimize them better
  !-------------------------------------------------------
  ! BODY
  !-------------------------------------------------------
  t1 = MPI_wtime()
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !-- define the various coordinate systems we are going to use
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           
           !-- call body subroutines
           call DrawBody(ix,iy,iz,x_body,color_body)
        enddo
     enddo
  enddo
  time_insect_body = time_insect_body + MPI_wtime() - t1
  
  !-------------------------------------------------------
  ! Eyes (not always present)
  !-------------------------------------------------------
  t1 = MPI_wtime()
  if (HasEye) then           
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !-- define the various coordinate systems we are going to use
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           x_eye_l  = x_body - xc_eye_l
           x_eye_r  = x_body - xc_eye_r
           
           !-- call eye subroutines
           call DrawEye(ix,iy,iz,x_eye_r,color_body)
           call DrawEye(ix,iy,iz,x_eye_l,color_body)
        enddo
     enddo
  enddo  
  endif
  time_insect_eye = time_insect_eye + MPI_wtime() - t1
  
  !-------------------------------------------------------
  ! Head (not always present)
  !-------------------------------------------------------
  t1 = MPI_wtime()
  if (HasHead) then           
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !-- define the various coordinate systems we are going to use
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           x_head   = x_body - xc_head
           !-- call body subroutines
           call DrawHead(ix,iy,iz,x_head,color_body)
        enddo
     enddo
  enddo  
  endif
  time_insect_head = time_insect_head + MPI_wtime() - t1
  
  !-------------------------------------------------------
  ! Wings (two subfunctions, for simple wings and those described
  ! by Fourier series)
  !-------------------------------------------------------
  t1 = MPI_wtime()
  if (fourier_wing) then
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !-- define the various coordinate systems we are going to use
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           x_wing_l = matmul(M_wing_l,x_body-xc_pivot_l)
           x_wing_r = matmul(M_wing_r,x_body-xc_pivot_r)
           
           !-- call wing subroutines
           call DrawWing_Fourier(ix,iy,iz,x_wing_l,M_wing_l,rot_l,color_l)
           call DrawWing_Fourier(ix,iy,iz,x_wing_r,M_wing_r,rot_r,color_r)
        enddo
     enddo
  enddo  
  
  else
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)
           !-- define the various coordinate systems we are going to use
           x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           x_body   = matmul(M_body,x-xc_body)
           x_wing_l = matmul(M_wing_l,x_body-xc_pivot_l)
           x_wing_r = matmul(M_wing_r,x_body-xc_pivot_r)
           
           !-- call wing subroutines
           call DrawWing_simple(ix,iy,iz,x_wing_l,M_wing_l,rot_l,color_l)
           call DrawWing_simple(ix,iy,iz,x_wing_r,M_wing_r,rot_r,color_r)
        enddo
     enddo
  enddo   
  endif
  time_insect_wings = time_insect_wings + MPI_wtime() - t1
  
  !-------------------------------------------------------
  ! Add solid body rotation (i.e. the velocity field that originates
  ! from the body rotation and translation. Until now, the wing velocities
  ! were the only ones set plus they are in the body reference frame
  !-------------------------------------------------------
  t1 = MPI_wtime()
  do ix = ra(1), rb(1)
     do iy = ra(2), rb(2)
        do iz = ra(3), rb(3)           
           ! add solid body rotation in the body-reference frame
           if (mask(ix,iy,iz) > 0.d0) then
            ! add solid body rotation to the translational velocity field
            v_tmp(1) = vc_body(1)+rot_body(2)*x_body(3)-rot_body(3)*x_body(2)
            v_tmp(2) = vc_body(2)+rot_body(3)*x_body(1)-rot_body(1)*x_body(3)
            v_tmp(3) = vc_body(3)+rot_body(1)*x_body(2)-rot_body(2)*x_body(1)
            
            us(ix,iy,iz,1:3)=matmul(M_body_inv,us(ix,iy,iz,1:3)+v_tmp)
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
  real (kind=pr), intent (in) :: time
  real (kind=pr) :: accx,accz,mass_solid,mass_fluid,gravity
  real (kind=pr) :: anglegs,kzlegsmax,dzlegsmax,t0,tmaxlegs,kzlegs0,anglegsend
  real (kind=pr) :: kxlegs,kzlegs,fxlegs,fzlegs,fxaero,fzaero,displacement_z

  ! mass
  mass_solid = Insect%mass_solid
  mass_fluid = 0.0d0  ! TODO

  ! gravity acceleration
  gravity = Insect%gravity

  select case (Insect%BodyMotion)
  case ("takeoff")
    if (Insect%KineFromFile=="simplified_dynamic") then

      ! Current body center displacement
      displacement_z = SolidDyn%var_new(2)

      anglegsend = Insect%anglegsend
      kzlegsmax = Insect%kzlegsmax 
      dzlegsmax = Insect%dzlegsmax
      t0 = Insect%t0legs
      tmaxlegs = t0 + Insect%tlinlegs
      kzlegs0 = (mass_solid-mass_fluid)*(-gravity) / dzlegsmax  

      ! Legs force
      fxlegs = 0.0d0
      fzlegs = 0.0d0
      if (Insect%ilegs>0) then
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
  integer :: mpicode  
  
  ! dynamics solver inactive by default
  idynamics = 0
  
  
  select case (Insect%BodyMotion)
  case ("takeoff")
    if (Insect%KineFromFile=="simplified_dynamic") then
      if (inicond(1:8).ne."backup::") then
        !--  we are not resuming a backup
        idynamics = 1
        SolidDyn%var_new(1) = 0.0d0
        SolidDyn%var_new(2) = 0.0d0
        SolidDyn%var_new(3) = 0.0d0
        SolidDyn%var_new(4) = 0.0d0
      else
        !-- we are resuming a backup
        idynamics = 1
        !-- root rank reads in backup file
        if (mpirank==0) then
          !-- backup files are called "runtime_backup0.h5.rigidsolver"
          write (*,*) "------"
          write (*,*) "Insect solver is resuming from file="//inicond(9:len_trim(inicond))//".rigidsolver"
          !-- open file
          open(10, file=inicond(9:len_trim(inicond))//".rigidsolver", form='formatted', status='old') 
          read(10, *) SolidDyn%var_new, SolidDyn%var_this, SolidDyn%rhs_this, SolidDyn%rhs_old
          write (*,*) SolidDyn%var_new
          write (*,*) SolidDyn%var_this
          write (*,*) SolidDyn%rhs_this
          write (*,*) SolidDyn%rhs_old
          !-- close file            
          close(10)
          write (*,*) "------"
        endif
        call MPI_BCAST( SolidDyn%var_new,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
        call MPI_BCAST( SolidDyn%var_this,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
        call MPI_BCAST( SolidDyn%rhs_this,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
        call MPI_BCAST( SolidDyn%rhs_old,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )        
      endif
    endif
  end select

end subroutine dynamics_insect_init


! Compute aerodynamic power
subroutine aero_power(apowtotal)
  use fsi_vars
  use mpi
  implicit none

  integer :: color_body, color_l, color_r
  real(kind=pr), dimension(1:3) :: omrel, momrel
  real(kind=pr), intent(out) :: apowtotal

  ! colors for Diptera (one body, two wings)
  color_body = 1
  color_l = 2
  color_r = 3

  ! body is not driven directly, therefore the power is set to zero
  Insect%PartIntegrals(color_body)%APow = 0.0d0

  ! wings
  ! all torques are computed with respect to body centre, therefore,
  ! an extra term appears (force cross position of the pivot point)

  ! left wing
  ! relative angular velocity
  omrel = Insect%rot_l_glob - Insect%rot_body_glob

  ! compute moment with respect to the pivot point
  ! initialize it as the moment with respect to insect's centre point
  momrel = Insect%PartIntegrals(color_l)%Torque + &
           Insect%PartIntegrals(color_l)%Torque_unst

  ! add correction
  momrel(1) = momrel(1) + &
           Insect%PartIntegrals(color_l)%Force(3) * Insect%x_pivot_l_glob(2) - &
           Insect%PartIntegrals(color_l)%Force(2) * Insect%x_pivot_l_glob(3)

  momrel(2) = momrel(2) + &
           Insect%PartIntegrals(color_l)%Force(1) * Insect%x_pivot_l_glob(3) - &
           Insect%PartIntegrals(color_l)%Force(3) * Insect%x_pivot_l_glob(1)

  momrel(3) = momrel(3) + &
           Insect%PartIntegrals(color_l)%Force(2) * Insect%x_pivot_l_glob(1) - &
           Insect%PartIntegrals(color_l)%Force(1) * Insect%x_pivot_l_glob(2)

  ! aerodynamic power
  Insect%PartIntegrals(color_l)%APow = - sum( momrel * omrel )

  ! right wing
  ! relative angular velocity
  omrel = Insect%rot_r_glob - Insect%rot_body_glob

  ! compute moment with respect to the pivot point
  ! initialize it as the moment with respect to insect's centre point
  momrel = Insect%PartIntegrals(color_r)%Torque + &
           Insect%PartIntegrals(color_r)%Torque_unst

  ! add correction
  momrel(1) = momrel(1) + &
           Insect%PartIntegrals(color_r)%Force(3) * Insect%x_pivot_r_glob(2) - &
           Insect%PartIntegrals(color_r)%Force(2) * Insect%x_pivot_r_glob(3)

  momrel(2) = momrel(2) + &
           Insect%PartIntegrals(color_r)%Force(1) * Insect%x_pivot_r_glob(3) - &
           Insect%PartIntegrals(color_r)%Force(3) * Insect%x_pivot_r_glob(1)

  momrel(3) = momrel(3) + &
           Insect%PartIntegrals(color_r)%Force(2) * Insect%x_pivot_r_glob(1) - &
           Insect%PartIntegrals(color_r)%Force(1) * Insect%x_pivot_r_glob(2)

  ! aerodynamic power
  Insect%PartIntegrals(color_r)%APow = - sum( momrel * omrel )

  ! Total aerodynamic power
  apowtotal = Insect%PartIntegrals(color_body)%APow + &
   Insect%PartIntegrals(color_l)%APow + Insect%PartIntegrals(color_r)%APow

end subroutine aero_power

