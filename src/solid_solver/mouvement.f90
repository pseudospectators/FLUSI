!-------------------------------------------------------------------------------
! The subroutine describes the imposed motion of the beam 
! within the RELATIVE coordinate system.  The angle alpha and its time 
! derivatives are the roll angle around the z axis of a plate
!
! We note that in the 3D code, the offset of the beam (x0|y0) is always set to zero
! so the beam is at the origin of the relative coordinate system
!
! We note further that because of this, even a heaving beam has always x0|y0=0
! since it is the relative system that is translated. However, in that case,
! the leading edge velocity is non-zero
!-------------------------------------------------------------------------------
subroutine mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam)
  implicit none
  type(solid), intent(in) :: beam
  real (kind=pr), intent(out) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), intent(in) :: time
  ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real (kind=pr), dimension(1:6), intent(out) :: LeadingEdge
  
  real(kind=pr)::f,angle_max


  select case (imposed_motion_leadingedge)
  case ("fixed_middle")
     !--------------------------------------------------------------------------
     ! fixed
     !--------------------------------------------------------------------------
     LeadingEdge = 0.0
     LeadingEdge(1) = beam%x0
     LeadingEdge(2) = beam%y0
     alpha    = beam%AngleBeam*pi/180.0
     alpha_t  = 0.0
     alpha_tt = 0.0
     
  case ("pitching_middle") 
      LeadingEdge = 0.0
      LeadingEdge(1) = beam%x0
      LeadingEdge(2) = beam%y0
      f = 1.0
      angle_max = deg2rad(45.d0)
      
      alpha    = angle_max * sin(2.d0*pi*f*time)
      alpha_t  = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
      alpha_tt = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
      
  case default
      if (mpirank==0) write(*,*) "mouvement:: imposed_motion_leadingedge undefined"
      if (mpirank==0) write(*,*) imposed_motion_leadingedge
      stop
      
  end select

end subroutine mouvement

!-------------------------------------------------------------------------------
! Relative coordinate system in which the plate lives. The deflection line is 
! always in the x-y plane. The vector x0 is the origin of the system, and the 
! angles psi, beta, gamma are the rotation angles (same notation as for insects
! psi = roll
! beta = pitch
! gamma = yaw
!-------------------------------------------------------------------------------
subroutine plate_coordinate_system( time, x0,v0, psi, beta, gamma, psi_dt, beta_dt, gamma_dt )
  implicit none
  real(kind=pr),dimension(1:3),intent(out) :: x0,v0
  real(kind=pr),intent(out) :: psi,beta,gamma,psi_dt,beta_dt,gamma_dt
  real(kind=pr), intent(in) :: time
  
  select case (imposed_motion_leadingedge)
  case ("fixed_middle")  
      !-- beam is in the middle of the domain and bends in x-y direction
      !-- z direction is height
      x0 = (/ 0.5*xl,0.5*yl,0.5*zl /)
      v0 = 0.d0
      psi = deg2rad(45.d0)*sin(time)
      beta = 0.d0
      gamma = 0.d0
      psi_dt = deg2rad(45.d0)*cos(time)
      beta_dt = 0.d0
      gamma_dt = 0.d0
      
  case ("pitching_middle") 
      !-- beam is in the middle of the domain and bends in x-y direction
      !-- z direction is height
      x0 = (/ 0.5*xl,0.5*yl,0.5*zl /)
      v0 = 0.d0
      psi = 0.d0
      beta = 0.d0
      gamma = 0.d0
      psi_dt = 0.0
      beta_dt = 0.d0
      gamma_dt = 0.d0
      
  case default
      if (mpirank==0) write(*,*) "plate_coordinate_system:: imposed_motion_leadingedge undefined"
      if (mpirank==0) write(*,*) imposed_motion_leadingedge      
      stop
      
  end select  
end subroutine plate_coordinate_system