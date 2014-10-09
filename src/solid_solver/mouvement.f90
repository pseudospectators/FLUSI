!-------------------------------------------------------------------------------
! The subroutine describes the imposed motion of the beam 
! within the RELATIVE coordinate system.  The angle alpha and its time 
! derivatives are the roll angle around the z axis of a plate
!
! We note that in the 3D code, the offset of the beam (x0_plate|y0) is always set to zero
! so the beam is at the origin of the relative coordinate system
!
! We note further that because of this, even a heaving beam has always x0_plate|y0=0
! since it is the relative system that is translated. However, in that case,
! the leading edge velocity is non-zero
!-------------------------------------------------------------------------------
subroutine mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam)
  implicit none
  type(solid), intent(in) :: beam
  real(kind=pr), intent(out) :: alpha, alpha_t, alpha_tt 
  real(kind=pr), intent(in) :: time
  ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real(kind=pr), dimension(1:6), intent(out) :: LeadingEdge
  real(kind=pr) :: f,angle_max, R
  real(kind=pr) :: a,b,c,d,k,kt,ktt,y,yt,ytt, t
  
   LeadingEdge = 0.0

  select case (imposed_motion_leadingedge)
  case ("fixed_x0y0z0")
     !--------------------------------------------------------------------------
     ! fixed
     !--------------------------------------------------------------------------
     alpha    = deg2rad( beam%AngleBeam )
     alpha_t  = 0.0
     alpha_tt = 0.0!!
     
  case ("swimmer") 
      ! swimmer including the leading edge cylinder
      f = frequ !-- normalizaton -> f is unity
      angle_max = deg2rad(AngleBeam)
      alpha    = angle_max * sin(2.d0*pi*f*time)
      alpha_t  = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
      alpha_tt = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
  
      !-- leading edge acceleration
      LeadingEdge(5) = -R_cylinder*alpha_tt*dsin(alpha) + &
                        R_cylinder*(alpha_t**2)*(-dcos(alpha))
      LeadingEdge(6) =  R_cylinder*alpha_tt*dcos(alpha) + &
                        R_cylinder*(alpha_t**2)*(-dsin(alpha))
                        
  case ("swimmer_simplified","swimmer_simplified_2D") 
      ! simplified swimmer without the leading edge cylinder
      ! startup conditioner is applied on the first period
      t = time * frequ 
      if (t <= 1.0) then
        a = -20.d0; b= 70.d0; c=-84.d0; d=35.d0;
        k    = a*t**7 + b*t**6 + c*t**5 + d*t**4
        kt  = 7.d0*a*t**6 + 6.d0*b*t**5 + 5.d0*c*t**4 + 4.d0*d*t**3
        ktt = 42.d0*a*t**5 + 30.d0*b*t**4 + 20.d0*c*t**3 + 12.d0*d*t**2     
      else
        k = 1.d0; kt = 0.d0; ktt = 0.d0
      endif     
      !-- actual motion protocoll 
      angle_max = deg2rad(AngleBeam)  
      y    = angle_max * sin(2.d0*pi*frequ*time)
      yt   = angle_max * cos(2.d0*pi*frequ*time) * (2.d0*pi*frequ)
      ytt  = -1.d0 * angle_max * sin(2.d0*pi*frequ*time) * (2.d0*pi*frequ)**2
      
      !-- motion protocoll times startup conditioner
      alpha = k*y
      alpha_t = kt*y + yt*k
      alpha_tt = ktt*y + ytt*k + 2.d0*kt*yt 

  case ("flapper")
     R=1.d0
     LeadingEdge = 0.0 ! note that both x,y and u,v are zero (v0_plate contains the velocity)
     ! however, leading edge acceleration is not zero
     LeadingEdge(5) = 0.d0
     LeadingEdge(6) = R*deg2rad(45.d0)*dsin(time)
     alpha    = 0.0
     alpha_t  = 0.0
     alpha_tt = 0.0
     
  case ("turek")
     alpha    = 0.0
     alpha_t  = 0.0
     alpha_tt = 0.0    
     
  case default
      if (mpirank==0) write(*,*) "mouvement:: imposed_motion_leadingedge undefined"
      if (mpirank==0) write(*,*) imposed_motion_leadingedge
      call abort()
      
  end select

end subroutine mouvement

!-------------------------------------------------------------------------------
! Relative coordinate system in which the plate lives. The deflection line is 
! always in the x-y plane. The vector x0_plate is the origin of the system, and the 
! angles psi, beta, gamma are the rotation angles (same notation as for insects
! psi = roll
! beta = pitch
! gamma = yaw
!-------------------------------------------------------------------------------
subroutine plate_coordinate_system( time, x0_plate,v0_plate, psi, beta, gamma, &
               psi_dt, beta_dt, gamma_dt, M_plate )
  implicit none
  
  real(kind=pr),dimension(1:3),intent(out) :: x0_plate,v0_plate
  real(kind=pr),dimension(1:3,1:3),intent(out)::M_plate
  real(kind=pr),intent(out)::psi,beta,gamma,psi_dt,beta_dt,gamma_dt
  real(kind=pr),intent(in)::time
  real(kind=pr)::R, angle_max,f, alpha, alpha_t, alpha_tt
  real(kind=pr),dimension(1:3,1:3) :: M1,M2,M3
  
  select case (imposed_motion_leadingedge)
  case ("fixed_x0y0z0")  
      !-- beam is in the middle of the domain and bends in x-y direction
      !-- z direction is height
      x0_plate = (/ x0,y0,z0 /)
      v0_plate = 0.d0
      psi = 0.d0
      beta = 0.d0
      gamma = 0.d0
      psi_dt = 0.0
      beta_dt = 0.d0
      gamma_dt = 0.d0
      
  case ("turek")
      x0_plate = (/ 0.5d0*xl, 1.220287d0,  0.8291429d0 /)
      if (nx==1) x0_plate(1)=0.d0
      v0_plate = 0.d0
      psi = deg2rad(+90.d0)
      beta = 0.0
      gamma = deg2rad(+90.d0)
      psi_dt = 0.0
      beta_dt = 0.d0
      gamma_dt = 0.d0
      
  case ("swimmer") 
      !-- beam is in the middle of the domain and bends in x-y direction
      !-- z direction is height      
      f = frequ !-- normalizaton -> f is unity
      angle_max = deg2rad(AngleBeam)      
      alpha    = angle_max * sin(2.d0*pi*f*time)
      alpha_t  = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
      alpha_tt = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
      !-- note (/ x0,y0,z0 /) marks center of cylinder
      x0_plate = (/ x0,y0,z0 /) + R_cylinder*(/dcos(alpha),dsin(alpha),0.d0/)
      v0_plate = R_cylinder*alpha_t*(/-dsin(alpha),dcos(alpha),0.d0/)
      !-- no rotation of relative system in swimmer case
      psi = 0.d0
      beta = 0.d0
      gamma = 0.d0
      psi_dt = 0.d0
      beta_dt = 0.d0
      gamma_dt = 0.d0
      
  case ("swimmer_simplified") 
      x0_plate = (/ x0,y0,z0 /)
      v0_plate = 0.d0

      psi = 0.d0
      beta = 0.d0
      gamma = 0.d0
      psi_dt = 0.d0
      beta_dt = 0.d0
      gamma_dt = 0.d0

  case ("swimmer_simplified_2D") 
      x0_plate = (/ 0.d0,y0,z0 /)
      v0_plate = 0.d0

      psi = deg2rad(+90.d0)
      beta = 0.0
      gamma = deg2rad(+90.d0)
      psi_dt = 0.0
      beta_dt = 0.d0
      gamma_dt = 0.d0    
      
  case ("flapper")
      R = 1.0 ! you need to change that in the above routine as well
      psi = deg2rad(45.d0)*dsin(time)
      psi_dt = deg2rad(45.d0)*dcos(time)
      x0_plate = (/ 0.5*xl, 0.5*yl+R*dsin(-psi), 0.5*zl+R*dcos(-psi) /)
      v0_plate = (/0.d0, -R*psi_dt*dcos(-psi), +R*psi_dt*dsin(-psi)/)
      beta = 0.d0
      gamma = 0.d0
      beta_dt = 0.d0
      gamma_dt = 0.d0
          
  case default
      if (mpirank==0) write(*,*) "plate_coordinate_system:: imposed_motion_leadingedge undefined"
      if (mpirank==0) write(*,*) imposed_motion_leadingedge      
      call abort()
      
  end select  
  
  if (mpirank==0) then
   if (maxval(v0_plate)>0.d0) then
   write(*,*) "thomas, please be sure to check if the angular velocities are okay here so that you don't do the&
   & same mistake twice."
   call abort()
   endif
  endif
   
  
  
  !-- rotation matrices that take us to the relative system
  call Rx(M1,psi)
  call Ry(M2,beta)
  call Rz(M3,gamma)
  M_plate = matmul(M1,matmul(M2,M3))
  
end subroutine plate_coordinate_system
