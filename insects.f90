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
  use mpi_header
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
  
  Insect%safety = 2.d0*dx
  
  ! some checks
  if ((mpirank==0).and.((iMoving.ne.1).or.(iPenalization.ne.1))) then
    write (*,*) "insects.f90::DrawInsect: the parameters iMoving or iPenalization are wrong."
    stop
  endif
  
  ! initialize mask and solid velocity as zero
  mask = 0.d0
  us = 0.d0
  
  !------------------------------------
  ! this is the relative coordinates (in the body system)
  ! of some interesting points on the beam (head, eyes, hinges)
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
  call Rx(M2,phi_l)
  call Rz(M3,theta_l)  
  M_wing_l = matmul(M1,matmul(M2,matmul(M3,M_stroke_l)))

  ! note the coordinate system is rotated so we don't need to inverse the sign
  ! of theta, and the wings still rotate in opposite direction
  call Ry(M1,-alpha_r)
  call Rx(M2,-phi_r)
  call Rz(M3,theta_r)  
  M_wing_r = matmul(M1,matmul(M2,matmul(M3,M_stroke_r)))
  
  !------------------------------------
  ! angular velocity vectors
  !------------------------------------
  rot_l = (/ phi_dt_l, alpha_dt_l, theta_dt_l /)
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
            
            us(ix,iy,iz,1:3)=us(ix,iy,iz,1:3)+matmul(transpose(M_body),v_tmp)
           endif
        enddo
    enddo
  enddo
end subroutine Draw_Insect




!------------------------------------------------------------------------------
! Draws a wing
! here, a wing is a rigid plate of constant thickness that differs from
! a rectangular plate only in the x-direction
! 
! note to save a bit of computing time, we first check the easy
! conditions (thickness and spanwise length) and then the shape
! function since this saves many evaluations of the shape.
subroutine DrawWing(ix,iy,iz,x_wing,M,rot)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr) :: a_body, R, R0, steps, x_top, x_bot
  real(kind=pr) :: y_tmp, x_tmp, z_tmp
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x_wing(1:3), rot(1:3), M(1:3,1:3)

  ! spanwise length:
  if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
  ! thickness: (note left and right wing have a different orientation of the z-axis
  ! but this does not matter since this is the same.
  if (abs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
      
      
  ! wing shape (determine between which x-values (x_bot, x_top) the wing is
  ! these values depend on the spanwise direction (which is y)
  select case(Insect%WingShape)
    case ('TwoEllipses')
      a_body = 0.5d0 * Insect%L_span
      if ((1.d0 - ((x_wing(2)-a_body)**2)/(a_body**2)) >= 0.d0) then
      x_top =  dsqrt((Insect%b_top**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
      x_bot = -dsqrt((Insect%b_bot**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
      else
      x_top = 0.d0
      x_bot = 0.d0
      endif
    case ('rectangular')
      x_top = Insect%b_top
      x_bot =-Insect%b_bot
    case default
      write (*,*) 'DrawWing:: unknown type of wing shape..'
      stop
  end select
      
  
  ! in the x-direction, the actual wing shape plays.    
  if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then        
      
    ! smooth length
    if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
      y_tmp = steps(-x_wing(2),0.d0,dz)
    else
      y_tmp = steps( x_wing(2),Insect%L_span,dz)
    endif

    ! smooth height
    z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness,dz) ! thickness       

    ! smooth shape
    if (x_wing(1)<0.d0) then
      x_tmp = steps(-x_wing(1),-x_bot,dz)
    else
      x_tmp = steps( x_wing(1), x_top,dz)
    endif
    
    mask(ix,iy,iz) = max(z_tmp*y_tmp*x_tmp, mask(ix,iy,iz))
    us(ix,iy,iz,1:3) = matmul(transpose(M), &
                      (/ rot(2)*x_wing(3)-rot(3)*x_wing(2), &
                         rot(3)*x_wing(1)-rot(1)*x_wing(3), &
                         rot(1)*x_wing(2)-rot(2)*x_wing(1)/) ) 
  endif
  
  
  
  endif
  endif
  
end subroutine DrawWing





!------------------------------------------------------------------------------
! Draws an insect's body, several options available.
! the body is, in the local coordinate system, always aligned with the
! x-axis. also, we currently use only rotational symmetric bodies.
subroutine DrawBody(ix,iy,iz,x_body)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr) :: a_body, R, R0, steps
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x_body(1:3)
  
  select case (Insect%BodyType)
  case ('ellipsoid')
    ! ellipsoid body    
    a_body = Insect%L_body / 2.d0
    
    ! check if inside the surrounding box (save comput. time)
    if ( dabs(x_body(2)) <= Insect%b_body + Insect%safety ) then
    if ( dabs(x_body(3)) <= Insect%b_body + Insect%safety ) then
    
    ! check for length inside ellipsoid:
    if ( dabs(x_body(1) ) < Insect%L_body/2 + Insect%safety ) then
        R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
        ! this gives the R(x) shape
        if ( (x_body(1)/a_body)**2 <= 1.d0) then
        R0 = dsqrt( Insect%b_body**2 *(1.d0- (x_body(1)/a_body)**2 ) )

        if ( R < R0 + Insect%safety ) then
            mask(ix,iy,iz)= max(steps(R,R0,dz),mask(ix,iy,iz))
        endif
        endif
    endif
    endif
    endif
  case ('nobody')
    ! doesn't do anything
  case default
    if (mpirank==0) then
    write (*,*) "In DrawBody: unkown body type."
    stop
    endif
  end select
  
end subroutine



!------------------------------------------------------------------------------
! Draws a sphere with radius R, as we need for the head and the eyes, if present
subroutine DrawSphere(ix,iy,iz,x,R0)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr), intent(in) :: R0
  real(kind=pr) :: R, steps
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  
  if (abs(x(1))<R0+Insect%safety) then
  if (abs(x(2))<R0+Insect%safety) then
  if (abs(x(3))<R0+Insect%safety) then
      R = sqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
      if ( R <= R0+Insect%safety ) then
        mask(ix,iy,iz) = max(steps(R,R0,dz),mask(ix,iy,iz))
      endif
  endif
  endif
  endif  
end subroutine

! as long as we have only spherical eyes, this is just a wrapper
subroutine DrawEye(ix,iy,iz,x)
  use fsi_vars
  use mpi_header
  implicit none
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  if (Insect%HasEye=="yes") then
  call DrawSphere(ix,iy,iz,x,Insect%R_eye)
  endif
end subroutine

! as long as we have only spherical heads, this is just a wrapper
subroutine DrawHead(ix,iy,iz,x)
  use fsi_vars
  use mpi_header
  implicit none
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  if (Insect%HasHead=="yes") then
  call DrawSphere(ix,iy,iz,x,Insect%R_head)
  endif
end subroutine



!-------------------------------------------------------------------------------
! Body motion protocoll, different choices.
! Input: 
!       time (self explanatory)
! Output:
!       psi:      roll angle
!       beta:     pitch angle
!       gamma:    yaw angle
!       psi_dt:   roll angular velocity
!       beta_dt:  pitch angular velocity
!       gamma_dt: yaw angular velocity
!       xc:       center of gravity coordinate
!       vc:       translational velocity of the body
! The actual motion depends on the choices in the parameter file, namely
! Insect%BodyMotion, and sub-parameters that may further precise a given motion 
! protocoll
subroutine BodyMotion(time, psi, beta, gamma, psi_dt, beta_dt, gamma_dt, xc, vc)
  use fsi_vars
  use mpi_header
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt
  real(kind=pr), intent(out) :: xc(1:3), vc(1:3)
  real(kind=pr) :: f,T, R
  
  select case (Insect%BodyMotion)
  case ("fixed")
    psi      = 0.d0
    beta     = 0.d0
    gamma    = 0.d0
    psi_dt   = 0.d0
    beta_dt  = 0.d0
    gamma_dt = 0.d0
    xc = (/0.5*xl, 0.5*yl,0.5*zl/)
    vc = (/0.0, 0.0, 0.0/)
    
  case ("wheeling")
    T = 20.0 ! time to do one turn    
    R = 1.5  ! circle radius
    
    psi      = deg2rad(-30.d0)
    beta     = 0.0
    gamma    = (2.d0*pi/T )*time
    psi_dt   = 0.0
    beta_dt  = 0.0
    gamma_dt = 2.d0*pi/T  
    
    xc = (/R*cos(1.5*pi+gamma)+0.5*xl, R*sin(1.5*pi+gamma)+0.5*yl, 0.5*zl/)
    vc = (/-R*sin(1.5*pi+gamma)*gamma_dt, R*cos(1.5*pi+gamma)*gamma_dt,0.d0/)
  case ("hovering")
    psi      = 0.0
    beta     = 0.0
    gamma    = 0.0
    psi_dt   = 0.0
    beta_dt  = 0.0
    gamma_dt = 0.0  
    xc = (/0.0, 0.0, 0.0/)
    vc = (/0.0, 0.0, 0.0/)
    
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::BodyMotion: motion case (Insect%BodyMotion) undefined"
    stop
    endif
  end select
  
end subroutine BodyMotion




!-------------------------------------------------------------------------------
! Flapping wing motion protocoll, different choices.
! Input: 
!       time (self explanatory)
!       protocoll: string containing what motion you want to have (may be 
!                  different for both wings)
! Output:
!       phi: flapping/positonal angle
!       alpha: feathering angle / angle of attack
!       theta: deviation angle / Out-of-stroke-plane
!       phi_dt: flapping angular velocity
!       alpha_dt: feathering angular velocity
!       theta_dt: deviatory angular velocity
! The actual motion depends on the choices in the parameter file, namely
! Insect%WingMotion, and sub-parameters that may further precise a given motion 
! protocoll. Note we allow both wings to follow a differen motion, but they both 
! call this routine here.
subroutine FlappingMotion(time, protocoll, phi, alpha, theta, phi_dt, alpha_dt, theta_dt)
  use fsi_vars
  use mpi_header
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  character (len=*), intent(in) :: protocoll
  real(kind=pr) :: phi_max,alpha_max, phase,f
  
  select case ( protocoll )
  case ("simplified")
    phi_max     = 60.d0*pi/180.d0  ! phi is up/down angle (flapping)
    alpha_max   = 45.d0*pi/180.d0  ! alpha is tethering
    phase       = 10.d0*pi/180.d0  ! phase shift between flapping and tethering
    f = 1.d0*2.0*pi
    phi      = phi_max  *dcos(f*time)
    alpha    = alpha_max*dsin(f*(time+phase))
    theta    = 0.0
    phi_dt   =-phi_max *f *dsin(f*time)
    alpha_dt = alpha_max*f*dcos(f*(time+phase))
    theta_dt = 0.0
  case ("debug")
    phi      = deg2rad(45.d0)   
    alpha    = deg2rad(0.d0)
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::FlappingMotion: motion case (protocoll) undefined"
    stop
    endif    
  end select
  
end subroutine FlappingMotion



!-------------------------------------------------------------------------------
! Stroke plane
! Input: 
!       time 
! Output:
!       eta_stroke: stroke plane angle 
subroutine StrokePlane ( time, eta_stroke )
  use fsi_vars
  use mpi_header
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: eta_stroke

  eta_stroke = deg2rad(0.d0) !+ 0.5d0*pi
  
end subroutine StrokePlane




!-----------------------
! Motion protocoll wrapper left wing
!-----------------------
subroutine FlappingMotion_left ( time, phi, alpha, theta, phi_dt, alpha_dt, theta_dt )
  use fsi_vars
  use mpi_header
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  call FlappingMotion ( time, Insect%FlappingMotion_left, &
                        phi, alpha, theta, phi_dt, alpha_dt, theta_dt )  
end subroutine FlappingMotion_left


!-----------------------
! Motion protocoll wrapper right wing
!-----------------------
subroutine FlappingMotion_right ( time, phi, alpha, theta, phi_dt, alpha_dt, theta_dt )
  use fsi_vars
  use mpi_header
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  call FlappingMotion ( time, Insect%FlappingMotion_left, &
                        phi, alpha, theta, phi_dt, alpha_dt, theta_dt )  
end subroutine FlappingMotion_right









!-----------------------
! Rotation Matrices
!-----------------------
subroutine Rx (R,angle)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R(1:3,1:3)
  R(1,:) = (/ 1.d0, 0.d0, 0.d0/)
  R(2,:) = (/ 0.d0, cos(angle), sin(angle) /)
  R(3,:) = (/ 0.d0, -sin(angle), cos(angle) /)
end subroutine 


subroutine Ry (R,angle)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R(1:3,1:3)
  R(1,:) = (/ cos(angle), 0.d0, -sin(angle)/)
  R(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
  R(3,:) = (/ +sin(angle), 0.d0, cos(angle) /)
end subroutine 


subroutine Rz (R,angle)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr), intent (in) :: angle
  real(kind=pr),dimension(1:3,1:3), intent(out) :: R(1:3,1:3)
  R(1,:) = (/ cos(angle), +sin(angle), 0.d0/)
  R(2,:) = (/ -sin(angle), cos(angle), 0.d0/)
  R(3,:) = (/ 0.d0, 0.d0, 1.d0  /)  
end subroutine 


real(kind=pr) function steps(x,t,h)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr) :: f,x,t,h
  call smoothstep(f,x,t,h)
  steps=f
end function