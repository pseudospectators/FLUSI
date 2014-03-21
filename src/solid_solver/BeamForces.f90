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

    call mouvement ( time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )

    tau_beam      = 0.0
    pressure_beam = 0.0
    sigma_beam    = 0.0

    theta = beam(:,5) + alpha

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

      pressure_beam(n) = DeltaInterpolation (xu1, yu1, press, 0.0, 0.0, xl-dx,yl-dy,"phi_4star") &
            - DeltaInterpolation (xb1, yb1, press, 0.0, 0.0, xl-dx,yl-dy,"phi_4star")   
      if (viscous_tensions) then    
          A1 = DeltaInterpolation (xu, yu, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
          B1 = DeltaInterpolation (xu, yu, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
          A2 = DeltaInterpolation (xb, yb, stress_a, 0.0,0.0,xl-dx,yl-dy,"phi_4star")
          B2 = DeltaInterpolation (xb, yb, stress_b, 0.0,0.0,xl-dx,yl-dy,"phi_4star")             
      endif
      
    enddo

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
    
  end subroutine GetForcesInterp

end module BeamForces


