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
  
  Insect%safety = 2.d0*dz
  Insect%smooth = 1.d0*dz
  
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
            
            us(ix,iy,iz,1:3)=matmul(transpose(M_body),us(ix,iy,iz,1:3)+v_tmp)
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
  real(kind=pr) :: a_body, R, R0, steps, x_top, x_bot, R_tmp
  real(kind=pr) :: y_tmp, x_tmp, z_tmp, xroot,yroot, f,xc,yc
  real(kind=pr) :: ai(1:40), bi(1:40), a0, theta
  real(kind=pr) :: v_tmp(1:3), mask_tmp
  integer, intent(in) :: ix,iy,iz
  integer :: i
  real(kind=pr),intent(in) :: x_wing(1:3), rot(1:3), M(1:3,1:3)

  select case (Insect%WingShape) 
  ! in these two cases, we have two given x_w(y_w) that delimit the wing
  case ('TwoEllipses','rectangular')  
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
    end select
        
    
    ! in the x-direction, the actual wing shape plays.    
    if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then        
        
      ! smooth length
      if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
        y_tmp = steps(-x_wing(2),0.d0)
      else
        y_tmp = steps( x_wing(2),Insect%L_span)
      endif

      ! smooth height
      z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness       

      ! smooth shape
      if (x_wing(1)<0.d0) then
        x_tmp = steps(-x_wing(1),-x_bot)
      else
        x_tmp = steps( x_wing(1), x_top)
      endif
      
      mask_tmp = z_tmp*y_tmp*x_tmp
      
      if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then 
        mask(ix,iy,iz) = mask_tmp
        !------------------------------------------------
        ! solid body rotation
        ! Attention: the Matrix transpose(M) brings us back to the body
        ! coordinate system, not to the inertial frame. this is done in 
        ! the main routine Draw_Insect
        !------------------------------------------------
        v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
        v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
        v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)
        
        ! note we set this only if it is a part of the wing
        us(ix,iy,iz,1:3) = matmul(transpose(M), v_tmp)
      endif
    endif  
    
    endif
    endif
  
  
  ! in this case, we have given the wing shape as a function R(theta) which is 
  ! given by some Fourier coefficients
  case ('drosophila')
    ! first, check if the point lies inside the rectanglee L_span x L_span
    ! here we assume that the chordlength is NOT greater than the span
    if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
    if ((x_wing(1)>=-(Insect%L_span+Insect%safety))&
        .and.&
        (x_wing(1)<=Insect%L_span+Insect%safety))&
    then
    if (abs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
      ! Fourier coefficients
      a0 = 0.5140278
      ai = (/0.1276258,-0.1189758,-0.0389458,0.0525938,0.0151538,-0.0247938,&
             -0.0039188,0.0104848,-0.0030638,-0.0064578,0.0042208,0.0043248,&
             -0.0026878,-0.0021458,0.0017688,0.0006398,-0.0013538,-0.0002038,&
             0.0009738,0.0002508,-0.0003548,-0.0003668,-0.0002798,0.0000568,&
             0.0003358,0.0001408,-0.0002208,0.0000028,0.0004348,0.0001218,&
             -0.0006458,-0.0003498,0.0007168,0.0003288,-0.0007078,-0.0001368,&
             0.0007828,0.0001458,-0.0007078,-0.0001358/) 
             
      bi = (/-0.1072518,-0.0449318,0.0296558,0.0265668,-0.0043988,-0.0113218,&
             -0.0003278,0.0075028,0.0013598,-0.0057338,-0.0021228,0.0036178,&
             0.0013328,-0.0024128,-0.0007688,0.0011478,0.0003158,-0.0005528,&
             0.0000458,0.0003768,0.0002558,0.0000168,-0.0006018,-0.0006338,&
             0.0001718,0.0007758,0.0001328,-0.0005888,-0.0001088,0.0006298,&
             0.0000318,-0.0008668,-0.0000478,0.0009048,0.0001198,-0.0008248,&
             -0.0000788,0.0007028,-0.0000118,-0.0006608/)
             
      ! wing root point        
      xroot =+0.1122
      yroot =-0.0157
      ! center of circle
      xc =-0.1206 + xroot
      yc = 0.3619 + yroot            
      ! normalized angle
      theta = atan2 (x_wing(2)-yc,x_wing(1)-xc )
      theta = ( theta + pi ) / (2.d0*pi)
      
      ! fourier series
      R0 = a0/2.0
      f = 2.d0*pi      
      do i = 1, 40
        R0 = R0 + ai(i)*dcos(f*dble(i)*theta) + bi(i)*dsin(f*dble(i)*theta)
      enddo

      R = sqrt ( (x_wing(1)-xc)**2 + (x_wing(2)-yc)**2 )
      R_tmp = steps(R,R0)
      
      z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness
      
      
      mask_tmp = z_tmp*R_tmp
      
      if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then 
        mask(ix,iy,iz) = mask_tmp
        !------------------------------------------------
        ! solid body rotation
        ! Attention: the Matrix transpose(M) brings us back to the body
        ! coordinate system, not to the inertial frame. this is done in 
        ! the main routine Draw_Insect
        !------------------------------------------------
        v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
        v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
        v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)
        ! note we set this only if it is a part of the wing
        us(ix,iy,iz,1:3) = matmul(transpose(M), v_tmp)
      endif
      
    endif
    endif
    endif
  end select
end subroutine DrawWing





!------------------------------------------------------------------------------
! Draws an insect's body, several options available.
! the body is, in the local coordinate system, always aligned with the
! x-axis. also, we currently use only rotational symmetric bodies.
subroutine DrawBody(ix,iy,iz,x_body)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr) :: a_body, R, R0, steps, x, x_tmp, R_tmp
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x_body(1:3)
  
  select case (Insect%BodyType)
  case ('ellipsoid')  
    ! ------------------------------------
    ! ellipsoid body (jerry)
    ! ------------------------------------
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
          mask(ix,iy,iz)= max(steps(R,R0),mask(ix,iy,iz))
        endif
        endif
    endif
    endif
    endif
    
    
  case ('drosophila')
    ! ------------------------------------
    ! two b-splines body (abdomen+thorax)
    ! ------------------------------------    
    x = x_body(1) + 0.8067 ! centers the thickest part of the thorax at the origin
    
    ! check if inside body bounds (in x-direction)
    if ( (x>=-Insect%safety) .and. (x<=1.2+Insect%safety) ) then    
      R0=0.0
      ! compute radius as a function of x (counting from the tail on)
      if (x < 0.6333) then
        ! we're in the ABDOMEN
        R0 = max( -1.2990*x**2 + 0.9490*x + 0.0267, 0.d0)
      elseif ((x >= 0.6333) .and. (x <=1.0 )) then
        ! we're in the THORAX 
        R0 = max( -2.1667*x**2 + 3.4661*x - 1.2194, 0.d0)
      elseif ((x >= 1.0) .and. (x <=1.2 )) then
        ! we're in the HEAD
        R0 = max( -12.68*x**2 + 27.4960*x - 14.7360, 0.d0)
      endif
    
      ! radius at this point
      R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
      
      ! smoothing in x-direction
      if (x<Insect%safety) then  ! xs is chordlength coordinate
        x_tmp = steps(-x, Insect%smooth)
      else
        x_tmp = steps( x,1.2-Insect%smooth)
      endif     
      
      
      if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
        R_tmp = steps(R,R0)        
        mask(ix,iy,iz)= max( R_tmp*x_tmp , mask(ix,iy,iz) )
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
        mask(ix,iy,iz) = max(steps(R,R0),mask(ix,iy,iz))
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
  
!     select case (Insect%BodyType)
!     case ('ellipsoid')  
      ! an ellipsoid body goes with a spherical head
      call DrawSphere(ix,iy,iz,x,Insect%R_head)
!     case ('drosophila')
      ! drosophilae have different heads.
!     end select
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
    beta     = deg2rad(-55.d0)
    gamma    = 0.0
    psi_dt   = 0.0
    beta_dt  = 0.0
    gamma_dt = 0.0  
    xc = (/0.5*xl, 0.5*yl,0.5*zl/)
    vc = (/0.0, 0.0, 0.0/)
    
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::BodyMotion: motion case (Insect%BodyMotion) undefined"
    stop
    endif
  end select
  
  ! for compability, we update the x0,y0,z0 also
  ! this is used e.g. for torque computation  
  x0 = xc(1)
  y0 = xc(2)
  z0 = xc(3)
  
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
  real(kind=pr) :: ai_phi(1:10), bi_phi(1:10), ai_theta(1:10), bi_theta(1:10)
  real(kind=pr) :: ai_alpha(1:10), bi_alpha(1:10)
  real(kind=pr) :: a0_alpha, a0_phi, a0_theta, s,c
  integer :: i
  
  select case ( protocoll )
  case ("Drosophila_hovering_fry")
    !---------------------------------------------------------------------------
    ! motion protocoll digitalized from Fry et al JEB 208, 2303-2318 (2005)
    !
    ! fourier coefficients analyzed with matlab
    !---------------------------------------------------------------------------
    a0_phi   =25.4649398
    a0_alpha =-0.3056968
    a0_theta =17.8244658
    ai_phi   =(/71.1061858,2.1685448,-0.1986978,0.6095268,-0.0311298,&
               -0.1255648,-0.0867778,0.0543518,0.0,0.0/)
    bi_phi   =(/5.4547058,-3.5461688,0.6260698,0.1573728,-0.0360498,-0.0205348,&
               -0.0083818,-0.0076848,0.0,0.0/)
    ai_alpha =(/3.3288788,0.6303878,-10.9780518,2.1123398,-3.2301198,&
               -1.4473158,0.6141758,-0.3071608,0.1458498,0.0848308/)
    bi_alpha =(/67.5430838,0.6566888,9.9226018,3.9183988,-2.6882828,0.6433518,&
                -0.8792398,-0.4817838,0.0300078,-0.1015118/)
    ai_theta =(/3.9750378,8.2808998,-0.0611208,-0.3906598,0.4488778,-0.120087,&
               -0.0717048,0.0699578,0.0,0.0/)
    bi_theta =(/2.2839398,3.5213068,-1.9296668,1.0832488,0.3011748,-0.1786648,&
                0.1228608,-0.0004808,0.0,0.0/)
    
    ! mean values
    phi = a0_phi/2.0
    alpha = a0_alpha/2.0
    theta = a0_theta/2.0
    
    phi_dt = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
    
    ! frequency
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,10
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi(i)   * c + bi_phi(i)   * s
      theta = theta + ai_theta(i) * c + bi_theta(i) * s
      alpha = alpha + ai_alpha(i) * c + bi_alpha(i) * s
      
      ! you checked this in matlab, it is correct.
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi(i)   * s + bi_phi(i)   * c)
      theta_dt = theta_dt + f*dble(i)*(-ai_theta(i) * s + bi_theta(i) * c)
      alpha_dt = alpha_dt + f*dble(i)*(-ai_alpha(i) * s + bi_alpha(i) * c)
    enddo
    
    if(mpirank == 0) then
    open(14,file='motion.time',status='unknown',position='append')
    write (14,'(7(e12.5,1x))') time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt
    close(14)
    endif
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = -deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = -deg2rad(theta_dt)
    
  case ("simplified")
    !---------------------------------------------------------------------------
    ! simplified motion protocoll
    !
    ! J. comput. Phys. 231 (2012) 1822-1847 "A fluid–structure interaction 
    ! model of insect flight with flexible wings"
    !
    ! the pase shift "phase" was my idea
    !---------------------------------------------------------------------------
    phi_max     = 60.d0*pi/180.d0  ! phi is up/down angle (flapping)
    alpha_max   = 0.d0!45.d0*pi/180.d0  ! alpha is tethering
    phase       = 0.d0! 10.d0*pi/180.d0  ! phase shift between flapping and tethering
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
  case ("none")
    phi      = 0.0
    alpha    = 0.0
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

  select case (Insect%BodyMotion)
  case ("fixed")
    eta_stroke = deg2rad(0.d0)    
  case ("wheeling")
    eta_stroke = deg2rad(0.d0)
  case ("hovering")
    eta_stroke = deg2rad(-35.d0)
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::StrokePlane: motion case (Insect%BodyMotion) undefined"
    stop
    endif
  end select
  
  
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
  call FlappingMotion ( time, Insect%FlappingMotion_right, &
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


! short for the smooth step function.
! the smooting is defined in Insect%smooth, here we need only x, and the 
! thickness (i.e., in the limit, steps=1 if x<t and steps=0 if x>t
real(kind=pr) function steps(x,t)
  use fsi_vars
  use mpi_header
  implicit none
  real(kind=pr) :: f,x,t
  call smoothstep(f,x,t,Insect%smooth)
  steps=f
end function