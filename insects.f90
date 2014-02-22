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
  mask = 0.0d0
  Insect%maskpart = 0.d0
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
  real(kind=pr) :: y_tmp, x_tmp, z_tmp, xroot,yroot, f,xc,yc, a0
  real(kind=pr), dimension(:), allocatable :: ai, bi
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta
  integer :: n_fft
  integer, intent(in) :: ix,iy,iz
  integer :: i
  real(kind=pr),intent(in) :: x_wing(1:3), rot(1:3), M(1:3,1:3)

  select case (Insect%WingShape) 
  
  !*****************************************************************************
  ! in these two cases, we have two given x_w(y_w) that delimit the wing
  !*****************************************************************************
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
        Insect%maskpart(ix,iy,iz,1) = mask_tmp ! For wing/body forces
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
    
  
  !*****************************************************************************
  ! in this case, we have given the wing shape as a function R(theta) which is 
  ! given by some Fourier coefficients
  !*****************************************************************************
  case ('drosophila','drosophila_mutated','drosophila_sandberg',&
        'drosophila_maeda','flapper_sane')
    ! first, check if the point lies inside the rectanglee L_span x L_span
    ! here we assume that the chordlength is NOT greater than the span
    if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
    if ((x_wing(1)>=-(Insect%L_span+Insect%safety)).and.(x_wing(1)<=Insect%L_span+Insect%safety)) then
    if (abs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
    
      !-----------------------------------------
      ! hard-coded Fourier coefficients for R(theta)
      !-----------------------------------------
      if (Insect%WingShape == 'drosophila') then
        !********************************************
        ! Drosophila wing from Jan Gruber's png file
        !********************************************
        n_fft = 40
        allocate ( ai(1:n_fft), bi(1:n_fft) )
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
      elseif (Insect%WingShape == 'drosophila_mutated') then
        !********************************************
        ! mutated Drosophila wing from Jan Gruber's png file
        !********************************************  
        n_fft = 70
        allocate ( ai(1:n_fft), bi(1:n_fft) )
        a0 = 0.4812548
        ai = (/0.1593968, -0.1056828, -0.0551518, 0.0508748, 0.0244538, -0.0264738,&
                -0.0080828, 0.0181228, 0.0023648, -0.0134578, -0.0037068, 0.0064508,&
                0.0028748, -0.0014258, -0.0006028, -0.0008898, -0.0020408, 0.0009218,&
                0.0029938, 0.0002768, -0.0026968, -0.0011518, 0.0017798, 0.0016538,&
                -0.0006098, -0.0012998, -0.0001918, 0.0003478, 0.0001408, 0.0003098,&
                0.0001078, -0.0005568, -0.0005998, 0.0006128, 0.0009078, -0.0003798,&
                -0.0009268, 0.0002128, 0.0009098, -0.0000598, -0.0010668, -0.0003428,&
                0.0009228, 0.0007688, -0.0003568, -0.0010458, -0.0004378, 0.0008738,&
                0.0009478, -0.0004108, -0.0012248, -0.0000638, 0.0013148, 0.0004978,&
                -0.0010638, -0.0007148, 0.0006338, 0.0007438, -0.0003278, -0.0006078,&
                0.0001838, 0.0003768, -0.0001698, -0.0002148, 0.0001318, 0.0001628,&
                -0.0000878, 0.0000068, 0.0001478, -0.0001128/) 
              
        bi = (/-0.1132588, -0.0556428, 0.0272098, 0.0221478, -0.0063798, -0.0059078,&
                  0.0043788, 0.0043208, -0.0003308, -0.0026598, -0.0013158, 0.0025178,&
                  0.0022438, -0.0023798, -0.0037048, 0.0001528, 0.0031218, 0.0022248,&
                  -0.0007428, -0.0027298, -0.0018298, 0.0014538, 0.0028888, 0.0000648,&
                  -0.0023508, -0.0009418, 0.0017848, 0.0016578, -0.0008058, -0.0017348,&
                  -0.0001368, 0.0011138, 0.0004218, -0.0005918, -0.0002798, 0.0002388,&
                  0.0002148, 0.0001408, 0.0000218, -0.0005138, -0.0003458, 0.0008208,&
                  0.0009888, -0.0007468, -0.0015298, 0.0002728, 0.0015588, 0.0002758,&
                  -0.0012498, -0.0006908,0.0008718, 0.0008848, -0.0003038, -0.0008048,&
                  -0.0001538, 0.0005418, 0.0003658, -0.0001988, -0.0003938, 0.0000048,&
                  0.0003008, 0.0000538, -0.0002748, -0.0000598, 0.0002898, 0.0001398,&
                  -0.0002108, -0.0001888, 0.0001838, 0.0001888 /)
              
        ! wing root point        
        xroot =+0.1122
        yroot =-0.0157
        ! center of circle
        xc =-0.1206 + xroot
        yc = 0.3619 + yroot        
      elseif (Insect%WingShape == 'drosophila_sandberg') then
        !********************************************
        !  Drosophila wing from Ramamurti & Sandberg ( JEB 210, 881-896, 2007)
        !********************************************        
        n_fft = 24 
        allocate ( ai(1:n_fft), bi(1:n_fft) )
        a0 = 0.4995578 
        ai = (/0.0164168,-0.1621518,0.0030938,0.0601108,-0.0083988,-0.0199988,&
        0.0049048,0.0047878,-0.0005648,-0.0001108,-0.0008638,-0.0006928,&
        0.0006608,0.0001978,0.0001558,0.0006878,-0.0007498,-0.0008018,&
        0.0003878,0.0007028,0.0000408,-0.0001108,-0.0001068,-0.0003958 &
        /)
        bi = (/-0.2083518,-0.0106488,0.0878308,-0.0018168,-0.0338278,0.0045768,&
        0.0113778,-0.0020678,-0.0026928,0.0002758,-0.0000838,-0.0001298,&
        0.0004118,0.0005638,-0.0001018,-0.0006918,-0.0002268,0.0005238,&
        0.0004008,-0.0001818,-0.0003038,-0.0000068,-0.0001218,0.0002008 &
        /)
        xc =-0.0235498  
        yc = 0.1531398 
      elseif (Insect%WingShape == 'drosophila_maeda') then
        !********************************************
        !  Drosophila wing from Maeda and Liu, similar to Liu and Aono, BB2009
        !********************************************        
        n_fft = 25
        allocate ( ai(1:n_fft), bi(1:n_fft) )
        a0 = 0.591294836514357
        ai = (/0.11389995408864588, -0.08814321795213981, -0.03495210456149335,&
        0.024972085605453047, 0.009422293191002384, -0.01680813499169695,&
        -0.006006435254421029, 0.012157932943676907, 0.00492283934032996,&
        -0.009882103857127606, -0.005421102356676356, 0.007230876076797827,&
        0.005272314598249222, -0.004519437431722127, -0.004658072133773225,&
        0.0030795046767766853, 0.003970792618725898, -0.0016315879319092456,&
        -0.002415442110272326, 0.0011118187761994598, 0.001811261693911865,&
        -2.6496695842951815E-4, -0.0012472769174353662, -1.7427507835680091E-4,&
        0.0010049640224536927/)
        bi = (/0.0961275426181888, 0.049085916171592914, -0.022051083533094627,&
        -0.014004783021121204, 0.012955446778711292, 0.006539648525493488,&
        -0.011873438993933363, -0.00691719567010525, 0.008479044683798266,&
        0.0045388280405204194, -0.008252172088956379, -0.005091347100627815,&
        0.004626409662755484, 0.004445034936616318, -0.0030708884306814804,&
        -0.004428808427471962, 0.0014113707529017868, 0.003061279043478891,&
        -8.658653756413232E-4, -0.002153349816945423, 3.317570161883452E-4,&
        0.001573518502682025, 2.14583094242007E-4, -0.0011299834277813852,&
        -5.172854674801216E-4/)
        !xc = 0.0 ! original mesh 
        xc = 0.0473 ! shifted towards t.e. to 1/4 of the root chord ("+" sign here)
        !xc = -0.0728 ! shifted towards l.e., to 0.2cmean from the l.e. (Liu and Aono BB 2009)
        yc = 0.7
      elseif (Insect%WingShape == 'flapper_sane') then
        !********************************************
        !  Mechanical model from Sane and Dickinson, JEB 205, 2002 
        !  'The aerodynamic effects...'
        !********************************************        
        n_fft = 25
        allocate ( ai(1:n_fft), bi(1:n_fft) )
        a0 = 0.5379588906565078
        ai = (/0.135338653455782,-0.06793162622123261,-0.0398235167675977,&
        0.006442194893963269,0.0012783260416583853,-0.007014398516674715,&
        0.0017710765408983137,0.006401601802033519,-2.970619204124993E-4,&
        -0.0038483478773981405,-6.180958756568494E-4,8.015784831786756E-4,&
        -6.957513357109226E-4,-1.4028929172227943E-4,0.0013484885717868547,&
        4.827827498543977E-4,-9.747844462919694E-4,-5.838504331939134E-4,&
        2.72834004831554E-4,2.8152492682871664E-5,-1.2802199282558645E-4,&
        4.117887216124469E-4,3.364169982438278E-4,-3.33258003686823E-4,&
        -3.5615733035757616E-4/)
        bi = (/2.686408368800394E-4,0.01649582345310688,0.01288513083639708,&
        0.004711436946785864,-0.0035725088809005073,-0.00898640397179334,&
        -0.003856509905612652,0.004536524572892801,0.004849677692836578,&
        2.9194421255236984E-4,-7.512780802871473E-4,7.12685261783966E-4,&
        -1.5519932673320404E-4,-0.0012695469974603026,2.2861692091158138E-4,&
        0.0016461316319681953,5.257476721137781E-4,-7.686482830046961E-4,&
        -3.108879176661735E-4,2.2437540206568518E-4,-2.578427217327782E-4,&
        -2.5120263516966855E-4,4.1693453021778877E-4,3.9290173948150096E-4,&
        -1.9762601237675826E-4/)
        xc = 0.0  
        yc = 0.6
      endif
      
      !-----------------------------------------
      ! get normalized angle (theta)
      !-----------------------------------------
      theta = atan2 (x_wing(2)-yc,x_wing(1)-xc )
      theta = ( theta + pi ) / (2.d0*pi)
      
      !-----------------------------------------
      ! construct R by evaluating the fourier series
      !-----------------------------------------
      R0 = a0/2.0
      f = 2.d0*pi    
      do i = 1, n_fft
        R0=R0 + ai(i)*dcos(f*dble(i)*theta) + bi(i)*dsin(f*dble(i)*theta)
      enddo
      deallocate (ai, bi)
      
      !-----------------------------------------
      ! get smooth (radial) step function
      !-----------------------------------------
      R = sqrt ( (x_wing(1)-xc)**2 + (x_wing(2)-yc)**2 )
      R_tmp = steps(R,R0)
      
      ! smooth also the thicknes
      z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness
      mask_tmp = z_tmp*R_tmp      
      
      !-----------------------------------------
      ! set new value for mask and velocity us
      !-----------------------------------------
      if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then 
        mask(ix,iy,iz) = mask_tmp
        Insect%maskpart(ix,iy,iz,1) = mask_tmp ! For wing/body forces
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
  real(kind=pr) :: a_body, R, R0, steps, x, y, z, s, s1, x1, x_tmp, R_tmp
  real(kind=pr) :: rbc, x0bc, z0bc, thbc1, thbc2, xcs, zcs
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
          Insect%maskpart(ix,iy,iz,2) = mask(ix,iy,iz) ! For wing/body forces
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
        Insect%maskpart(ix,iy,iz,2) = mask(ix,iy,iz) ! For wing/body forces
      endif      
    
    endif
    
  
  
  case ('drosophila_maeda')
    ! ------------------------------------
    ! approximation to mesh from Maeda
    ! similar to Aono et al.
    ! ------------------------------------    
    x = x_body(1)
    y = x_body(2)
    z = x_body(3)

    ! symmetry plane is xz
    ! +x direction is forward
    ! body centerline is an arc with center at x0bc,y0bc
    ! radius rbc and angles th counted from negative z
    rbc = 0.9464435146443515
    thbc1 = 112.0d0 *pi/180.0d0
    thbc2 = 53.0d0 *pi/180.0d0
    x0bc = -0.24476987447698745d0
    z0bc = -0.9301255230125524
  
    ! chordwise dimensionless coordinate, from head to abdomen
    s = (atan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1) 
    ! body center coordinates at s
    xcs = x0bc + (x-x0bc)*rbc/sqrt((x-x0bc)**2+(z-z0bc)**2)
    zcs = z0bc + (z-z0bc)*rbc/sqrt((x-x0bc)**2+(z-z0bc)**2)
    ! distance to the body center at s
    R = sqrt( (x-xcs)**2 + y**2 + (z-zcs)**2 )

    ! check if inside body bounds (in s-direction)
    if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then    
      R0 = 0.0d0
      ! distortion of s
      s1 = 1.0d0 - ( s + 0.08d0*tanh(30.0d0*s) ) / (1.0d0+0.08d0*tanh(30.0d0))
      s1 = ( s1 + 0.04d0*tanh(60.0d0*s1) ) / (1.0d0+0.04d0*tanh(60.0d0))
      s1 = (sin(1.2d0*s1)/sin(1.2d0))**1.25d0
      x1 = 1.075d0 * s1 
      ! compute radius as a function of x1 (counting from the tail on)
      ! same shape as 'drosophila'
      if (x1 < 0.6333d0) then
        ! we're in the ABDOMEN
        R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
      elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
        ! we're in the THORAX 
        R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
      endif
      ! distortion of R0
      R0 = 0.8158996d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0

      ! smoothing
      if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
        R_tmp = steps(R,R0)        
        mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
        Insect%maskpart(ix,iy,iz,2) = mask(ix,iy,iz) ! For wing/body forces
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
  real(kind=pr) :: x_head,z_head,dx_head,dz_head,R,R0,steps
  if (Insect%HasHead=="yes") then
  
     select case (Insect%BodyType)
     case ('ellipsoid')  
      ! an ellipsoid body goes with a spherical head
      call DrawSphere(ix,iy,iz,x,Insect%R_head)

     case ('drosophila')
      ! drosophilae have different heads.

     case ('drosophila_maeda')  
      ! ellipsoid head, assumes xc_head=0 in .ini file
      x_head = 0.17d0
      z_head = -0.1d0
      dx_head = 0.5d0*  0.185d0
      dz_head = 0.5d0*  0.27d0
      ! check if inside the surrounding box (save comput. time)
      if ( dabs(x(2)) <= dz_head + Insect%safety ) then
      if ( dabs(x(3)-z_head) <= dz_head + Insect%safety ) then
      ! check for length inside ellipsoid:
      if ( dabs(x(1)-x_head) < dx_head + Insect%safety ) then
        R  = dsqrt ( x(2)**2 + (x(3)-z_head)**2 )
        ! this gives the R(x) shape
        if ( ((x(1)-x_head)/dx_head)**2 <= 1.d0) then
        R0 = dz_head*dsqrt(1.d0- ((x(1)-x_head)/dx_head)**2 )
        if ( R < R0 + Insect%safety ) then
          mask(ix,iy,iz)= max(steps(R,R0),mask(ix,iy,iz))
          Insect%maskpart(ix,iy,iz,2) = mask(ix,iy,iz) ! For wing/body forces
        endif
        endif
      endif
      endif
      endif

     case default
      ! do nothing
     end select 
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
  use share_kine 
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(out) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt
  real(kind=pr), intent(out) :: xc(1:3), vc(1:3)
  real(kind=pr) :: f,T,R
  
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
    
    xc = (/R*dcos(1.5d0*pi+gamma)+0.5d0*xl, R*dsin(1.5d0*pi+gamma)+0.5d0*yl, 0.5d0*zl/)
    vc = (/-R*dsin(1.5d0*pi+gamma)*gamma_dt, R*dcos(1.5d0*pi+gamma)*gamma_dt,0.d0/)

  case ("hovering")
    psi      = 0.0
!    beta     = deg2rad(-55.d0)
    beta     = deg2rad(-45.d0)  ! Comparison with Maeda (Dmitry, 7 Nov 2013)
    gamma    = deg2rad(45.d0)
    psi_dt   = 0.0
    beta_dt  = 0.0
    gamma_dt = 0.0  

!    xc = (/0.5*xl, 0.5*yl, 0.5*zl/)  ! Dmitry, 26 Oct 2013
!    xc = (/0.5*xl, 0.5*yl, zl-1.0d0/)  ! Dmitry, 30 Oct 2013 -one wing length from top
    xc = (/0.5*xl, 0.5*yl, zl-1.3d0/)  ! Dmitry, 30 Oct 2013 -1.3 wing length from top
!    xc = (/0.5d0*xl, 0.5d0*yl, 0.8d0/)  ! Dmitry, 28 Oct 2013  - ground dist+0.3
    vc = (/0.0d0, 0.0d0, 0.0d0/)    

  case ("flapper")  ! Comparison with Dickinson et al. (Dmitry, 19 Nov 2013)
    psi      = 0.0
    beta     = deg2rad(-90.d0)
    gamma    = deg2rad(45.d0)
    psi_dt   = 0.0
    beta_dt  = 0.0
    gamma_dt = 0.0  

    xc = (/0.5*xl, 0.5*yl, zl-1.0d0/)  
    vc = (/0.0d0, 0.0d0, 0.0d0/)    

  case ("takeoff")  ! Takeoff kinematics read from file (Dmitry, 14 Nov 2013)
    if (Insect%KineFromFile=="yes") then
      call body_kine_interp(time,beta,xc(3),xc(1),beta_dt,vc(3),vc(1))
      ! takeoff velocity factor
      !xc(3) = xc(3) * 0.1d0
      !vc(3) = vc(3) * 0.1d0
      !xc(1) = xc(1) * 1.0d0
      !vc(1) = vc(1) * 1.0d0
      ! x coordinate
      xc(1) = xc(1)+ 2.0d0 !0.5d0*xl
      ! y coordinate
      xc(2) = 0.5d0*yl
      vc(2) = 0.0d0
      ! vertical position corrected
      xc(3) = xc(3) + 0.3d0 + 0.56d0 !(ground+legs)
!      xc(3) = xc(3) + 0.3d0 + 2.0d0 !(far from the ground)
      ! convert pitch angle to flusi conventions
      beta = -beta
      beta = deg2rad(beta)
      beta_dt = -beta_dt
      beta_dt = deg2rad(beta_dt)
      ! zero heading and yaw
      psi = 0.0d0
      psi_dt = 0.0d0
      gamma = 0.0d0
      gamma_dt = 0.0d0
    elseif (Insect%KineFromFile=="simplified_dynamic") then
      ! interpolate. xc(3),xc(1),vc(3),vc(1) are unused!
      call body_kine_interp(time,beta,xc(3),xc(1),beta_dt,vc(3),vc(1))
      ! y coordinate
      xc(2) = 0.5d0*yl
      vc(2) = 0.0d0
      ! convert pitch angle to flusi conventions
      beta = -beta
      beta = deg2rad(beta)
      beta_dt = -beta_dt
      beta_dt = deg2rad(beta_dt)
      ! zero heading and yaw
      psi = 0.0d0
      psi_dt = 0.0d0
      gamma = 0.0d0
      gamma_dt = 0.0d0
      ! Use data from flight dynamics solver
      xc(1) = SolidDyn%var_new(1) + 2.0d0
!      xc(3) = SolidDyn%var_new(2) + 0.3d0 + 0.56d0 !(ground+legs)
      xc(3) = SolidDyn%var_new(2) + 0.3d0 + 2.0d0 !(far from the ground)
      vc(1) = SolidDyn%var_new(3)
      vc(3) = SolidDyn%var_new(4)
    endif

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
  real(kind=pr) :: bi_alpha_flapper(1:29) ! For comparison with Sane&Dickinson
  real(kind=pr) :: ai_phi_flapper(1:31) ! For comparison with Sane&Dickinson
  real(kind=pr) :: tadv ! For comparison with Dickinson
  real(kind=pr) :: posi,elev,feth,posi_dt,elev_dt,feth_dt,angles ! Comp. w. Maeda
  real(kind=pr) :: dangle_posi,dangle_elev,dangle_feth ! Comp. w. Maeda
  real(kind=pr) :: dangle_posi_dt,dangle_elev_dt,dangle_feth_dt ! Comp. w. Maeda
  real(kind=pr) :: a_posi(1:4),b_posi(1:4),a_elev(1:4),b_elev(1:4),a_feth(1:4),b_feth(1:4)
  real(kind=pr) :: a0_alpha, a0_phi, a0_theta, s,c
  real(kind=pr) :: tau, phia, la, ta, dtt, t1, phic, phicdeg, ua
  real(kind=pr) :: alphac, alphacdeg, dtr, tr0
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
    a0_theta =-17.8244658  ! - sign (Dmitry, 10 Nov 2013)
    ai_phi   =(/71.1061858,2.1685448,-0.1986978,0.6095268,-0.0311298,&
               -0.1255648,-0.0867778,0.0543518,0.0,0.0/)
    bi_phi   =(/5.4547058,-3.5461688,0.6260698,0.1573728,-0.0360498,-0.0205348,&
               -0.0083818,-0.0076848,0.0,0.0/)
    ai_alpha =(/3.3288788,0.6303878,-10.9780518,2.1123398,-3.2301198,&
               -1.4473158,0.6141758,-0.3071608,0.1458498,0.0848308/)
    bi_alpha =(/67.5430838,0.6566888,9.9226018,3.9183988,-2.6882828,0.6433518,&
                -0.8792398,-0.4817838,0.0300078,-0.1015118/)
    ai_theta =(/-3.9750378,-8.2808998,0.0611208,0.3906598,-0.4488778,0.120087,&
               0.0717048,-0.0699578,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
    bi_theta =(/-2.2839398,-3.5213068,1.9296668,-1.0832488,-0.3011748,0.1786648,&
               -0.1228608,0.0004808,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
    
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
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)

!    if(mpirank == 0) then
!    open(14,file='motion.t',status='unknown',position='append')
!    write (14,'(7(e12.5,1x))') time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt
!    close(14)
!    endif

  case ("Drosophila_hovering_maeda")
    !---------------------------------------------------------------------------
    ! Drosophila hovering kinematics protocol 
    !
    ! Fourier coefficients provided by Maeda
    ! Diditized from Fry et al.
    !---------------------------------------------------------------------------
    a_posi = (/  0.22700d0,  1.24020d0,  0.03610d0, -0.00360d0/)
    b_posi = (/  0.00000d0,  0.08880d0, -0.07000d0,  0.01250d0/)
    a_elev = (/  0.16125d0,  0.06750d0,  0.14500d0,  0.00540d0/)
    b_elev = (/  0.00000d0,  0.03670d0,  0.06840d0, -0.03390d0/)
    a_feth = (/ -0.00864d0, -0.04890d0, -0.02056d0,  0.19649d0/)
    b_feth = (/  0.00000d0, -1.17586d0, -0.01216d0, -0.17590d0/)

    ! Initialize angles and velocities
    posi = 0.0d0
    elev = 0.0d0
    feth = 0.0d0
    posi_dt = 0.0d0
    elev_dt = 0.0d0
    feth_dt = 0.0d0

    do i=0,3 !! Fourier series
      !! time dependent angle
      angles  = 2.0d0*dble(i)*pi*time

      selectcase( i )
      case( 0 ) ! Fourier 0th order
        ! mean
        dangle_posi = a_posi(1) ! +shift_mean_posi_
        dangle_elev = a_elev(1) ! +shift_mean_elev_
        dangle_feth = a_feth(1) ! +shift_mean_feth_

        dangle_posi_dt = 0.0d0
        dangle_elev_dt = 0.0d0
        dangle_feth_dt = 0.0d0
  
      case default !! Fourier n-th orders
  
        call get_dangle( &
            & angles, &                !! intent(in)
            & i, &                     !! intent(in)
            & a_posi(i+1), & !! intent(in)
            & b_posi(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_posi, &      !! intent(out)
            & dangle_posi_dt &   !! intent(out)
        & )
  
        call get_dangle( &
            & angles, &                !! intent(in
            & i, &                     !! intent(in)
            & a_elev(i+1), & !! intent(in)
            & b_elev(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_elev, &      !! intent(out)
            & dangle_elev_dt &   !! intent(out)
        & )
  
        call get_dangle( &
            & angles, &                !! intent(in
            & i, &                     !! intent(in)
            & a_feth(i+1), & !! intent(in)
            & b_feth(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_feth, &      !! intent(out)
            & dangle_feth_dt &   !! intent(out)
        & )
  
      endselect
  
      posi = posi +dangle_posi
      elev = elev +dangle_elev
      feth = feth +dangle_feth

      posi_dt = posi_dt +dangle_posi_dt
      elev_dt = elev_dt +dangle_elev_dt
      feth_dt = feth_dt +dangle_feth_dt
    enddo
 
    ! Convert to FLUSI's variables
    phi = posi
    alpha = -feth
    theta = -elev
    
    phi_dt = posi_dt
    alpha_dt = -feth_dt
    theta_dt = -elev_dt
    
!    if(mpirank == 0) then
!    open(14,file='motion.t',status='unknown',position='append')
!    write (14,'(7(e12.5,1x))') time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt
!    close(14)
!    endif
    
  case ("flapper_sane")
    !---------------------------------------------------------------------------
    ! motion protocol from Sane and Dickinson, JEB 204, 2607-2626 (2001)
    !
    ! feathering: fourier coefficients analyzed with matlab, 2nd order
    !             Butterworth filter with cutoff at k=10
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 2 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in JEB 204, p. 2613 
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.807554373967804d0,&
     0.0d0,11.14661083909663d0,0.0d0,2.242734216805251d0,&
     0.0d0,-0.6141899985692184d0,0.0d0,-0.7426551158681146d0,&
     0.0d0,-0.2329560587573768d0,0.0d0,0.038749678276091284d0,&
     0.0d0,0.07083462320831221d0,0.0d0,0.028982501947490313d0,&
     0.0d0,-0.0025202918494477244d0,0.0d0,-0.010221019942802941d0,&
     0.0d0,-0.005614021318470698d0,0.0d0,1.1958884364596903d-6,&
     0.0d0,0.002186832241254999d0,0.0d0,0.0015347995090793172d0/)

    alpha = 0.0
    alpha_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of JEB 204 (eg alphedeg=90-50 for Fig 3D)
    alphacdeg = 90.0d0 - 00.0d0
    alpha = alphacdeg/40.0d0 * alpha
    alpha_dt = alphacdeg/40.0d0 * alpha_dt
   
    ! convert in radians 
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)
    
    ! *** II. position ***
    ai_phi_flapper =(/72.96795908179631d0,&
     0.0d0,8.064401876272864d0,0.0d0,2.769062401215844d0,&
     0.0d0,1.2200252377066352d0,0.0d0,0.5584689705779989d0,&
     0.0d0,0.2545617536476344d0,0.0d0,0.11829515180579572d0,&
     0.0d0,0.05754453975774996d0,0.0d0,0.02964141751269772d0,&
     0.0d0,0.016177705089515895d0,0.0d0,0.009315101869467001d0,&
     0.0d0,0.005625663922446026d0,0.0d0,0.0035424425357352385d0,&
     0.0d0,0.0023130422432356247d0,0.0d0,0.001558278163264511d0,&
     0.0d0,0.001078213692334021d0/)

    phi = 0.0
    phi_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo
    
    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0d0
    phi = phicdeg/180.0d0 * phi
    phi_dt = phicdeg/180.0d0 * phi_dt
   
    ! convert in radians 
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)
    
    ! *** III. elevation ***
    theta = 0.0d0
    theta_dt = 0.0d0


    !if(mpirank == 0) then
    !open(14,file='motion.t',status='unknown',position='append')
    !write (14,'(7(e12.5,1x))') time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt
    !close(14)
    !endif
    
  case ("flapper_dickinson")
    !---------------------------------------------------------------------------
    ! motion protocol from Dickinson, Lehmann and Sane, Science (1999)
    !
    ! feathering: fourier coefficients analyzed with matlab,
    !             Gaussian filter that fits fig 3D
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 5 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in Science
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.23094285611071d0,&
      0.0d0,10.224154661301371d0,0.0d0,2.1623763046726396d0,&
      0.0d0,0.05049394424178093d0,0.0d0,-0.17550942623071494d0,&
      0.0d0,-0.06634193748204852d0,0.0d0,-0.008925020495896451d0,&
      0.0d0,0.0011292567942149407d0,0.0d0,6.471071566666472d-4,&
      0.0d0,1.0018757795834964d-4,0.0d0,3.0105550216312524d-6,&
      0.0d0,-1.237567150768195d-6,0.0d0,-1.988004402010933d-7,&
      0.0d0,-1.10165545174181d-8,0.0d0,2.4135650975460306d-10/)

    ! Advanced rotation (+ sign) or delayed rotation (- sign)
    tadv = 0.08
    !tadv = - 0.08

    alpha = 0.0
    alpha_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+tadv)) 
      c = dcos(f*dble(i)*(time+tadv))
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of Science (eg alphedeg=90-40 for Fig 3)
    alphacdeg = 90.0d0 - 40.0d0
    alpha = alphacdeg/40.0d0 * alpha
    alpha_dt = alphacdeg/40.0d0 * alpha_dt
   
    ! convert in radians 
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)
    
    ! *** II. position ***
    ai_phi_flapper =(/63.24528806534019d0,&
      0.0d0,5.753991800610726d0,0.0d0,1.3887974015525626d0,&
      0.0d0,0.3889856512386744d0,0.0d0,0.10577402496901325d0,&
      0.0d0,0.026061339604144987d0,0.0d0,0.005623376646981709d0,&
      0.0d0,0.001042285996467963d0,0.0d0,1.639611509380189d-4,&
      0.0d0,2.1716252827442023d-5,0.0d0,2.408190194815521d-6,&
      0.0d0,2.2268710288534648d-7,0.0d0,1.7118916093759426d-8,&
      0.0d0,1.0914870312823793d-9,0.0d0,5.76135101855556d-11,&
      0.0d0,2.513944479978149d-12/)

    phi = 0.0
    phi_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo
    
    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0d0
    phi = phicdeg/180.0d0 * phi
    phi_dt = phicdeg/180.0d0 * phi_dt
   
    ! convert in radians 
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0d0
    theta_dt = 0.0d0

    !if(mpirank == 0) then
    !open(14,file='motion.t',status='unknown',position='append')
    !write (14,'(7(e12.5,1x))') time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt
    !close(14)
    !endif
   
  case ("takeoff")
    !--------------------------------------------------
    ! Fontaine et al. 
    !--------------------------------------------------
    if (Insect%KineFromFile/="no") then
      call wing_kine_interp(time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt)
      ! position angle
      phi = deg2rad(phi)
      phi_dt = deg2rad(phi_dt)
      ! feathering angle
      alpha = deg2rad(alpha)         
      alpha_dt = deg2rad(alpha_dt)    
      ! elevation angle in flusi coordinates
      theta = -theta
      theta_dt = - theta_dt
      theta = deg2rad(theta)
      theta_dt = deg2rad(theta_dt)
    endif 
 
  case ("simplified")
    !---------------------------------------------------------------------------
    ! simplified motion protocoll
    !
    ! J. comput. Phys. 231 (2012) 1822-1847 "A fluid-structure interaction 
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
!    eta_stroke = deg2rad(-35.d0)
    eta_stroke = deg2rad(-45.d0)  ! Comparison with Maeda (Dmitry, 7 Nov 2013)
  case ("flapper")   ! Comparison with Dickinson et al. (Dmitry, 19 Nov 2013)
    eta_stroke = deg2rad(0.d0)
  case ("takeoff")
    eta_stroke = deg2rad(-28.d0) ! 62-90, Fontaine et al., fig 13 (Dmitry, 14 Nov 2013)
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
  use mpi_header
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
      anglegsend = 0.25d0*pi ! case 1 and 3
!      anglegsend = 0.5d0*pi ! case2
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
  use mpi_header
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



