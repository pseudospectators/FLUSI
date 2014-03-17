module BeamForces
use motion
  implicit none
  contains

subroutine GetForces ( time, beams, press, u, timelevel )
  use share_vars
  use FieldExport
  implicit none
  type(solid), dimension(1:iBeam), intent(inout) :: beams
  real (kind=pr), intent (in)   :: time  
  real (kind=pr), intent (in)   :: u (0:nx-1, 0:ny-1,1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)    
  real (kind=pr)  :: pressure_beam (0:ns-1)
  real (kind=pr)  :: force_pressure (1:2)
  real (kind=pr)  :: beam (0:ns-1, 1:6)
  real (kind=pr)  :: tau_beam (0:ns-1)
  character(len=3), intent (in) :: timelevel
  integer :: i

  !-------------------------------------------------------------------
  ! This is a WRAPPER to call various interpolation routines for the beam
  ! so you can compare them
  !----------------------------------------------------------------------
  do i = 1, iBeam
      beam(:,1) = beams(i)%x
      beam(:,2) = beams(i)%y
      beam(:,3) = beams(i)%vx
      beam(:,4) = beams(i)%vy
      beam(:,5) = beams(i)%theta
      beam(:,6) = beams(i)%theta_dot
      
      call GetForcesInterp ( time, beam, pressure_beam, tau_beam,  press, u, force_pressure, 3 , i, beams(i) )
      
      if (timelevel=="new") then
        ! new time level
        beams(i)%pressure_new = pressure_beam
        beams(i)%tau_new      = tau_beam
        beams(i)%Force_press  = force_pressure
      elseif (timelevel=="old") then
        ! old time level
        beams(i)%pressure_old = pressure_beam
        beams(i)%tau_old      = tau_beam
        beams(i)%Force_press  = force_pressure
      else
        write(*,*) "error: timelevel unclear"
        stop
      endif      
  enddo

end subroutine GetForces





! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@





subroutine GetForcesInterp(time, beam, pressure_beam, tau_beam, press, u, force_pressure, SPL, beamnumber, beam_solid )
  ! tradional interpolation routine: set SPL=0 linear interpolation SPL=1 bciubic with FD derivatives (cheaper) SPL=2 bicubic with spectral derivatives
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer :: n
  integer :: i
  real (kind=pr), intent (in)   :: time
  integer       , intent (in)   :: SPL, beamnumber
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1), tau_beam(0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)
  real (kind=pr), intent (in)   :: u (0:nx-1, 0:ny-1,1:2)
  type (solid)                  :: beam_solid ! we need this for compatability
  real (kind=pr)                :: xu,yu,xb,yb,xu1,xb1,yu1,yb1
  real (kind=pr)                :: theta(0:ns-1), sigma_beam(0:ns-1)
  real (kind=pr)                :: alpha, alpha_t, alpha_tt, fx, fy,soft_startup,dx,dy,A1,B1,A2,B2,sigma_n1,sigma_n2,sigma_t2,sigma_t1
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  real (kind=pr)                :: press_dx (0:nx-1, 0:ny-1),press_dy (0:nx-1, 0:ny-1)
  real (kind=pr)                :: press_dxdy (0:nx-1, 0:ny-1),press_k (0:nx-1, 0:ny-1)
  real (kind=pr), dimension (0:nx-1, 0:ny-1)  :: temp, work1, work2, work3, stress_a, stress_b
  logical                       :: viscous_tensions
  character(len=1)              :: beamstr
  
  if (iViscous==1) viscous_tensions = .true.
  if (iViscous==0) viscous_tensions = .false.
  
  dx = xl/real(nx)
  dy = yl/real(ny)
  
  call mouvement ( time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  
  
  tau_beam      = 0.0
  pressure_beam = 0.0
  sigma_beam    = 0.0
  
  ! -----------------------------------------------------
  ! --- compute stress tensors
  ! -----------------------------------------------------
  if (viscous_tensions) then
    call coftxy  ( u(:,:,1), work1 )
    call cofdx   ( work1   , work2 ) 
    call cofitxy ( work2   , stress_a ) 
    stress_a = 2.0 * nu * stress_a
    !------------------------
    call coftxy  ( u(:,:,1), work1 )
    call cofdy   ( work1   , work2 ) 
    call cofitxy ( work2   , work3 )
    call coftxy  ( u(:,:,2), work1 )
    call cofdx   ( work1   , work2 ) 
    call cofitxy ( work2   , stress_b ) 
    stress_b = nu * ( work3 + stress_b )
    !   do ix=1,nx-2
    !   do iy=1,ny-2
    !     stress_a(ix,iy) = 2.0 * nu * (u(ix+1,iy,1)-u(ix-1,iy,1))/2.0/dx
    !     stress_b(ix,iy) = nu *(u(ix,iy+1,1)-u(ix,iy-1,1))/2.0/dy + nu *(u(ix+1,iy,2)-u(ix-1,iy,2))/2.0/dx
    !   enddo
    !   enddo
  endif
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (viscous_tensions) call BeamSurfaceIntegral (time, beam, u, stress_a, stress_b, press, beamnumber, beam_solid)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  

  !------------------------------------------------------------------------------------------------------------------------------
  !--           Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------

  if (SPL==2) then
    ! pressure
     call coftxy (press,press_k) 
     call cofdx  (press_k,temp)
     call cofitxy(temp,press_dx)     
     call cofdy  (press_k,temp)
     call cofitxy(temp,press_dy)     
     call cofdxdy(press_k,temp)
     call cofitxy(temp,press_dxdy)
  endif
  
  theta = beam(:,5) + alpha

  !$omp parallel do private(n,xu,yu,xb,yb,xu1,yu1,xb1,yb1,A1,A2,B1,B2,sigma_n1,sigma_t1,sigma_t2,sigma_n2) 
  do n = 0, ns-1
    !top point
    xu = beam(n,1)  - (t_beam+N_smooth*max(dy,dx)) * sin(theta(n))
    yu = beam(n,2)  + (t_beam+N_smooth*max(dy,dx)) * cos(theta(n))
    !bottom point
    xb = beam(n,1)  + (t_beam+N_smooth*max(dy,dx)) * sin(theta(n))
    yb = beam(n,2)  - (t_beam+N_smooth*max(dy,dx)) * cos(theta(n))
    !top point (pressure)
    xu1 = beam(n,1)  - (t_beam-0.0*N_smooth*max(dy,dx)) * sin(theta(n))
    yu1 = beam(n,2)  + (t_beam-0.0*N_smooth*max(dy,dx)) * cos(theta(n))
    !bottom point (pressure)
    xb1 = beam(n,1)  + (t_beam-0.0*N_smooth*max(dy,dx)) * sin(theta(n))
    yb1 = beam(n,2)  - (t_beam-0.0*N_smooth*max(dy,dx)) * cos(theta(n))    

    if (SPL==0) then            !--------------------------------------------------------------------------------------
    
    pressure_beam(n) = LinearInterpolation (xu1, yu1, press, 0.0,0.0,xl-dx,yl-dy) &
                     - LinearInterpolation (xb1, yb1, press, 0.0,0.0,xl-dx,yl-dy)  
                     
    elseif (SPL==1) then        !--------------------------------------------------------------------------------------
    
    pressure_beam(n) = BicubicInterpolationFD (xu1, yu1, press, 0.0,0.0,xl-dx,yl-dy) &
                     - BicubicInterpolationFD (xb1, yb1, press, 0.0,0.0,xl-dx,yl-dy) 
                  if (viscous_tensions) then   
                  A1 = BicubicInterpolationFD (xu, yu, stress_a, 0.0,0.0,xl-dx,yl-dy)
                  B1 = BicubicInterpolationFD (xu, yu, stress_b, 0.0,0.0,xl-dx,yl-dy)
                  A2 = BicubicInterpolationFD (xb, yb, stress_a, 0.0,0.0,xl-dx,yl-dy)
                  B2 = BicubicInterpolationFD (xb, yb, stress_b, 0.0,0.0,xl-dx,yl-dy) 
                  endif
    
    elseif (SPL==2) then        !--------------------------------------------------------------------------------------
    
    pressure_beam(n) = BicubicInterpolation (xu1, yu1, press,press_dx,press_dy,press_dxdy, 0.0, 0.0, xl-dx,yl-dy) &
                     - BicubicInterpolation (xb1, yb1, press,press_dx,press_dy,press_dxdy, 0.0, 0.0, xl-dx,yl-dy)   
                  if (viscous_tensions) then    
                  A1 = DeltaInterpolation (xu, yu, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  B1 = DeltaInterpolation (xu, yu, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  A2 = DeltaInterpolation (xb, yb, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  B2 = DeltaInterpolation (xb, yb, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")             
                  endif
                  
    elseif (SPL==3) then         !--------------------------------------------------------------------------------------
    
    pressure_beam(n) = DeltaInterpolation (xu1, yu1, press, 0.0, 0.0, xl-dx,yl-dy,"phi_4star") &
                     - DeltaInterpolation (xb1, yb1, press, 0.0, 0.0, xl-dx,yl-dy,"phi_4star")   
                  if (viscous_tensions) then    
                  A1 = DeltaInterpolation (xu, yu, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  B1 = DeltaInterpolation (xu, yu, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  A2 = DeltaInterpolation (xb, yb, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
                  B2 = DeltaInterpolation (xb, yb, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")             
                  endif           
    endif       !--------------------------------------------------------------------------------------         
    
    if (viscous_tensions) then 
      ! viscous stresses    
      sigma_n1 = - A1*cos(2.0*theta(n) ) - B1*sin(2.0*theta(n))
      sigma_t1 = - A1*sin(2.0*theta(n) ) + B1*cos(2.0*theta(n))       
      sigma_n2 = - A2*cos(2.0*theta(n) ) - B2*sin(2.0*theta(n))
      sigma_t2 = - A2*sin(2.0*theta(n) ) + B2*cos(2.0*theta(n))
      
      tau_beam(n)   = sigma_t1 - sigma_t2 
      sigma_beam(n) =-sigma_n1 + sigma_n2
    endif
  enddo   
  !$omp end parallel do
 
  if (viscous_tensions) then
    write (beamstr,'(i1)') beamnumber 
  
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'tau'//beamstr, status = 'unknown', access = 'append') ! Append output data file
    write (14,'(1x, 128(es15.8,1x))') time, (tau_beam(n), n=0,ns-1,4)
    close (14)
    
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'sigma'//beamstr, status = 'unknown', access = 'append') ! Append output data file
    write (14,'(1x, 128(es15.8,1x))') time, (sigma_beam(n), n=0,ns-1,4)
    close (14)
    
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'press'//beamstr, status = 'unknown', access = 'append') ! Append output data file
    write (14,'(1x, 128(es15.8,1x))') time, (pressure_beam(n), n=0,ns-1,4)
    close (14)
  endif
  
  
  if (viscous_tensions) then
    pressure_beam = pressure_beam + sigma_beam
  else
    pressure_beam = pressure_beam
    tau_beam = 0.0
  endif
  
  !------------------------------------------------------------------------------------------------------------------------------
  !--           startup (slow coupling or values from another simulation)
  !------------------------------------------------------------------------------------------------------------------------------ 
  if (time <= T_release) then
      soft_startup = 0.0
  elseif ( ( time >T_release ).and.(time<(T_release + tau)) ) then
      soft_startup =  ((time-T_release)**3)/(-0.5*tau**3)   + 3.*((time-T_release)**2)/tau**2
  else
      soft_startup = 1.0
  endif
  
  pressure_beam = pressure_beam*soft_startup 
  tau_beam      = tau_beam*soft_startup
    
  !------------------------------------------------------------------------------------------------------------------------------
  !--           Compute Drag and Lift resulting from the pressure on the beam (w/o viscous tensions)
  !------------------------------------------------------------------------------------------------------------------------------

  fx=0.0
  fy=0.0

  do i=0, ns-1
    fx = fx - (sin(theta(i))*pressure_beam(i) + cos(theta(i))*tau_beam(i))*ds
    fy = fy + (cos(theta(i))*pressure_beam(i) + sin(theta(i))*tau_beam(i))*ds
  enddo
  
  force_pressure(1) = fx
  force_pressure(2) = fy 
  
  !------------------------------------------------------------------------------------------------------------------------------
  !--           save energy contributions
  !------------------------------------------------------------------------------------------------------------------------------

  write (beamstr,'(i1)') beamnumber
  ! save force "integral" along top and bottom side. used for comparison with BeamSurfaceIntegral
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'int'//beamstr, status = 'unknown', access = 'append') ! Append output data file
  write (14,'(1x, 3(es15.8,1x))') time, fx,fy
  
  ! save energy of the external forces (that should be what's missing in the beam's energy)
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'external_power'//beamstr, status = 'unknown', access = 'append') ! Append output data file
  write (14,'(1x, 3(es15.8,1x))') time, &
  sum( pressure_beam * (-beam(:,3)*sin(theta)+beam(:,4)*cos(theta) ) )*ds,&
  sum( tau_beam      * ( beam(:,3)*cos(theta)+beam(:,4)*sin(theta) ) )*ds
  close (14)
  
end subroutine








!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@








subroutine BeamSurfaceIntegral (time, beam, u, stress_a, stress_b, press, beamnumber, beam_solid)
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer, parameter 		:: nc=128 ! number of points on a circle
  integer			:: i
  integer, intent(in)           :: beamnumber
  type (solid)                  :: beam_solid
  real (kind=pr), intent (in)   :: time
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (in)   :: stress_a (0:nx-1, 0:ny-1), stress_b(0:nx-1, 0:ny-1), u(0:nx-1, 0:ny-1,1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)  
  real (kind=pr)                :: alpha, alpha_t, alpha_tt, dx,dy, dphi, dL, beta, sigma_n, sigma_t, sigma_p, drag, lift, drag_surf, lift_surf, A, B, P
  real (kind=pr)		::  penal(0:nx-1, 0:ny-1,1:2)
  real (kind=pr), dimension(1:2)    :: Ft=0.0, Fn=0.0, Fp=0.0 !forces (tangential, normal, pressure)
  real (kind=pr), dimension(0:ns-1) :: xu,yu,xb,yb, theta, xu1,yu1,xb1,yb1 ! coordinates of the beam contours (up and down)
  real (kind=pr), dimension(0:nc-1) :: xc_trailing, yc_trailing, xc_leading, yc_leading, phi_trailing, phi_leading
  real (kind=pr), dimension(1:6)    :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  logical, save 	            :: warned=.false.
  character(len=1)                  :: beamstr
  
  dx = xl/real(nx)
  dy = yl/real(ny)
  
  if ((iCylinder ==1).and.(warned==.false.)) then
    write (33,*) "!!! Surface Integral cannot be computed with a cylinder at the leading edge..."    
    warned = .true.
  endif
  
  if (iCylinder==0) then
  
  
  
  Ft = 0.0
  Fn = 0.0
  Fp = 0.0
        
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  
 
  ! -----------------------------------------------------
  ! --- trailing edge half circle
  ! -----------------------------------------------------
  
  dphi = pi / real(nc-1)	! note that we include the first and last points
  beta = beam( ns-1, 5) + alpha ! angle for the cylinder
  
  do i=0,nc-1
    phi_trailing(i) = -( 0.5*pi-beta ) + real(i)*dphi
  enddo
  
  ! line element
  dL = dphi*(t_beam+N_smooth*max(dx,dy))
  ! coordinates on the circle
  xc_trailing = beam(ns-1,1) + (t_beam+N_smooth*max(dx,dy)) * cos(phi_trailing)
  yc_trailing = beam(ns-1,2) + (t_beam+N_smooth*max(dx,dy)) * sin(phi_trailing)
  
  do i=0, nc-1
!     A = BicubicInterpolationFD (xc_trailing(i), yc_trailing(i), stress_a, 0.0,0.0,xl-dx,yl-dy)
!     B = BicubicInterpolationFD (xc_trailing(i), yc_trailing(i), stress_b, 0.0,0.0,xl-dx,yl-dy)
!     P = BicubicInterpolationFD (xc_trailing(i), yc_trailing(i), press, 0.0,0.0,xl-dx,yl-dy)
    A = DeltaInterpolation (xc_trailing(i), yc_trailing(i), stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    B = DeltaInterpolation (xc_trailing(i), yc_trailing(i), stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    P = DeltaInterpolation (xc_trailing(i), yc_trailing(i), press, 0.0,0.0,xl-dx,yl-dy,"phi_4star")    
    
    
    sigma_n = + A*cos(2.0*phi_trailing(i) ) + B*sin(2.0*phi_trailing(i))
    sigma_t = - A*sin(2.0*phi_trailing(i) ) + B*cos(2.0*phi_trailing(i))
    sigma_p = - P
    
    Ft(1) = Ft(1) + dL*sigma_t * ( -sin(phi_trailing(i)) )
    Ft(2) = Ft(2) + dL*sigma_t * ( +cos(phi_trailing(i)) )
    Fn(1) = Fn(1) + dL*sigma_n * ( +cos(phi_trailing(i)) )
    Fn(2) = Fn(2) + dL*sigma_n * ( +sin(phi_trailing(i)) )    
    Fp(1) = Fp(1) + dL*sigma_p * ( +cos(phi_trailing(i)) )
    Fp(2) = Fp(2) + dL*sigma_p * ( +sin(phi_trailing(i)) )    
  enddo
  
  ! -----------------------------------------------------
  ! --- leading edge half circle
  ! -----------------------------------------------------
  
  dphi = pi / real(nc-1)	! note that we include the first and last points
  beta = beam( 1, 5) + alpha ! angle for the cylinder
  
  do i=0,nc-1
    phi_leading(i) = (0.5*pi+beta) + real(i)*dphi
  enddo
  
  ! line element
  dL = dphi*(t_beam+N_smooth*max(dx,dy))
  ! coordinates on the circle
  xc_leading = beam(0,1) + (t_beam+N_smooth*max(dx,dy)) * cos(phi_leading)
  yc_leading = beam(0,2) + (t_beam+N_smooth*max(dx,dy)) * sin(phi_leading)
  
  do i=0, nc-1
!     A = BicubicInterpolationFD (xc_leading(i), yc_leading(i), stress_a, 0.0,0.0,xl-dx,yl-dy)
!     B = BicubicInterpolationFD (xc_leading(i), yc_leading(i), stress_b, 0.0,0.0,xl-dx,yl-dy)
!     P = BicubicInterpolationFD (xc_leading(i), yc_leading(i), press, 0.0,0.0,xl-dx,yl-dy)
    A = DeltaInterpolation (xc_leading(i), yc_leading(i), stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    B = DeltaInterpolation (xc_leading(i), yc_leading(i), stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    P = DeltaInterpolation (xc_leading(i), yc_leading(i), press, 0.0,0.0,xl-dx,yl-dy,"phi_4star")    
    
    sigma_n = + A*cos(2.0*phi_leading(i) ) + B*sin(2.0*phi_leading(i))
    sigma_t = - A*sin(2.0*phi_leading(i) ) + B*cos(2.0*phi_leading(i))
    sigma_p = - P
    
    Ft(1) = Ft(1) + dL*sigma_t * ( -sin(phi_leading(i)) )
    Ft(2) = Ft(2) + dL*sigma_t * ( +cos(phi_leading(i)) )
    Fn(1) = Fn(1) + dL*sigma_n * ( +cos(phi_leading(i)) )
    Fn(2) = Fn(2) + dL*sigma_n * ( +sin(phi_leading(i)) )    
    Fp(1) = Fp(1) + dL*sigma_p * ( +cos(phi_leading(i)) )
    Fp(2) = Fp(2) + dL*sigma_p * ( +sin(phi_leading(i)) )    
  enddo  
  
  
  ! -----------------------------------------------------
  ! --- beam (top side)
  ! -----------------------------------------------------
  
  theta = beam(:,5) + alpha
  xu = beam(:,1) - (t_beam+N_smooth*max(dx,dy)) * sin( theta )
  yu = beam(:,2) + (t_beam+N_smooth*max(dx,dy)) * cos( theta )
  xu1 = beam(:,1) - t_beam * sin( theta )
  yu1 = beam(:,2) + t_beam * cos( theta )  
  
  do i=0,ns-1
!     A = BicubicInterpolationFD (xu(i), yu(i), stress_a, 0.0,0.0,xl-dx,yl-dy)
!     B = BicubicInterpolationFD (xu(i), yu(i), stress_b, 0.0,0.0,xl-dx,yl-dy)
!     P = BicubicInterpolationFD (xu(i), yu(i), press, 0.0,0.0,xl-dx,yl-dy)
    A = DeltaInterpolation (xu(i), yu(i), stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    B = DeltaInterpolation (xu(i), yu(i), stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    P = DeltaInterpolation (xu1(i), yu1(i), press, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    
    
    sigma_n = - A*cos(2.0*theta(i) ) - B*sin(2.0*theta(i))
    sigma_t = - A*sin(2.0*theta(i) ) + B*cos(2.0*theta(i))
    sigma_p = - P
    
    Ft(1) = Ft(1) + ds*sigma_t * ( +cos(theta(i)) )
    Ft(2) = Ft(2) + ds*sigma_t * ( +sin(theta(i)) )
    Fn(1) = Fn(1) + ds*sigma_n * ( -sin(theta(i)) )
    Fn(2) = Fn(2) + ds*sigma_n * ( +cos(theta(i)) )    
    Fp(1) = Fp(1) + ds*sigma_p * ( -sin(theta(i)) )
    Fp(2) = Fp(2) + ds*sigma_p * ( +cos(theta(i)) )    
  enddo
  
  ! -----------------------------------------------------
  ! --- beam (bottom side)
  ! -----------------------------------------------------
  xb = beam(:,1) + (t_beam+N_smooth*max(dx,dy)) * sin( theta )
  yb = beam(:,2) - (t_beam+N_smooth*max(dx,dy)) * cos( theta )
  xb1 = beam(:,1) + t_beam * sin( theta )
  yb1 = beam(:,2) - t_beam * cos( theta )  
  
  do i=0,ns-1
!     A = BicubicInterpolationFD (xb(i), yb(i), stress_a, 0.0,0.0,xl-dx,yl-dy)
!     B = BicubicInterpolationFD (xb(i), yb(i), stress_b, 0.0,0.0,xl-dx,yl-dy)
!     P = BicubicInterpolationFD (xb(i), yb(i), press, 0.0,0.0,xl-dx,yl-dy)
    A = DeltaInterpolation (xb(i), yb(i), stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    B = DeltaInterpolation (xb(i), yb(i), stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
    P = DeltaInterpolation (xb1(i), yb1(i), press, 0.0,0.0,xl-dx,yl-dy,"phi_4star")    
    
    sigma_n = - A*cos(2.0*theta(i) ) - B*sin(2.0*theta(i))
    sigma_t = - A*sin(2.0*theta(i) ) + B*cos(2.0*theta(i))
    sigma_p = - P
    
    Ft(1) = Ft(1) + ds*sigma_t * ( -cos(theta(i)) )
    Ft(2) = Ft(2) + ds*sigma_t * ( -sin(theta(i)) )
    Fn(1) = Fn(1) + ds*sigma_n * ( +sin(theta(i)) )
    Fn(2) = Fn(2) + ds*sigma_n * ( -cos(theta(i)) )    
    Fp(1) = Fp(1) + ds*sigma_p * ( +sin(theta(i)) )
    Fp(2) = Fp(2) + ds*sigma_p * ( -cos(theta(i)) )    
  enddo
  
  !-----------------------------------------------------------------------------
  penal(:,:,1) = mask * (u(:,:,1) - maskvx)  ! ux: x-component of the relative velocity
  penal(:,:,2) = mask * (u(:,:,2) - maskvy)  ! uy: y-component of the relative velocity
  drag = sum(penal(:,:,1))*dx*dy
  lift = sum(penal(:,:,2))*dx*dy
  
  drag_surf = Ft(1) + Fn(1) + Fp(1)
  lift_surf = Ft(2) + Fn(2) + Fp(2)
  
  write(beamstr,'(i1)') beamnumber
  
  
  if (iWalls==0) then
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'surface'//beamstr, status = 'unknown', access = 'append') ! Append output data file
    write (14,'(1x, 5(es15.8,1x),"  ",3(es15.8,1x),"  ", 3(es15.8,1x))') time, drag, lift, drag_surf, lift_surf, Ft(1),Fn(1),Fp(1),Ft(2),Fn(2),Fp(2)
    close (14)
  else
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'surface'//beamstr, status = 'unknown', access = 'append') ! Append output data file
    write (14,'(1x, 5(es15.8,1x),"  ",3(es15.8,1x),"  ", 3(es15.8,1x))') time, 0.0, 0.0, drag_surf, lift_surf, Ft(1),Fn(1),Fp(1),Ft(2),Fn(2),Fp(2)
    close (14)
  endif
  

  
  
  endif
  
end subroutine BeamSurfaceIntegral



















end module


