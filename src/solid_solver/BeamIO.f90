! ================================================================================================================================================
! 							Beam Input/Output subroutines
! ================================================================================================================================================



subroutine InitBeamFiles()
  use share_vars
  use motion
  implicit none
  character(len=1) :: beamstr
  integer :: i

  
  if (iSaveBeam == 1) then
  
       do i=1, iBeam
          write(beamstr,'(i1)') i   
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_x'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_y'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vx'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vy'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14)
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta_dot'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          close (14) 
          
          open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data'//beamstr, status = 'replace')
          call WriteHeaderToFile(14)
          write (14,'(A)') "%  beam data file, all quantities for the end point "
          write (14,'(A)') "%         time   pressure(ns-1)  pressure(3ns/4) pressure(ns/2)  x-coordinate    y-coordinate    x-velocity      y-velocity      theta           theta_dot       drag            lift            drag_unst       lift_unst       F_press_x       F_press_y       E_kin_solid     E_elastic_sol   E_potential   "
          close (14)
      enddo
  
  elseif (iSaveBeam == 3) then
       do i=1, iBeam
        write(beamstr,'(i1)') i      
        open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data'//beamstr, status = 'replace')
        call WriteHeaderToFile(14)
        write (14,'(A)') "%  beam data file, all quantities for the end point "
        write (14,'(A)') "%         time   pressure(ns-1)  pressure(3ns/4) pressure(ns/2)  x-coordinate    y-coordinate    x-velocity      y-velocity      theta           theta_dot       drag            lift            drag_unst       lift_unst       F_press_x       F_press_y       F_INERT_x       F_INERT_y       E_kin_solid     E_elastic_sol   E_potential   "
        close (14)      
        open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_complete'//beamstr, status = 'replace')
        call WriteHeaderToFile(14)
        close (14)
      enddo 
  
  endif
  

  
  do i = 1, iBeam
    write (beamstr,'(i1)') i
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'mouvement'//beamstr, status = 'replace')
    call WriteHeaderToFile(14)
    write (14,'(A)') " %         time  alpha           alpha_t         alpha_tt        x-pos           y-pos           x-vel          y-vel            x-accel        y-accel"
    close (14)  
  enddo
  
  
  
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'energies', status = 'replace')
  call WriteHeaderToFile(14)
  write (14,'(A)') "%         time   dt              vor_rms         E_kinetic_f     dissipation     Penalty diss    Inflow Solid    Inflow Mean   "  
  close (14)  
  
  
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'u_max', status = 'replace')
  call WriteHeaderToFile(14)
  write (14,'(A)') "%  time         dt           u_max        mean_ux      mean_uy"  
  close (14) 
  

  
  call RewriteStopFile
  
  !----------------------------------------------------------------------------------------------------------------------------
  open  (10, file = trim(dir_name)//'/'//trim(simulation_name)//"performance_details", status = 'replace')
  write (10,'(A)') "%    time_left     runtime         time        dt             sec / dt    sum         time_NST         time_pressure    time_mask        time_solid"
  close (10)
  
end subroutine InitBeamFiles

! =================================================================================================================================================================================


subroutine SaveBeamData( time, beams, dt1 )
  use share_vars
  use SolidSolver
  use motion
  implicit none
  real (kind=pr), intent (in) :: time,dt1
  type (solid), dimension(1:iBeam), intent (in) :: beams
  character(len=16) :: format_ns1
  character(len=3)  :: ns1_string
  character(len=1)  :: beamstr
  real (kind=pr) :: alpha, alpha_t, alpha_tt , dx, dy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  integer :: n,step,i
  
  dx=xl/real(nx)
  dy=yl/real(ny)  

  
  ! set up formats
  write(ns1_string, '(I3)') ns+1
  format_ns1 = '('//ns1_string//'(es12.5,1x))'
  
  
  !--------------------------------------------------------------
  ! save the entire beam data ?
  !--------------------------------------------------------------
  if (iSaveBeam==1) then !save the entire beam
  
    do i=1, iBeam
      write (beamstr,'(i1)') i   
      
      step=1 ! save the entire beam? or only every 2..3 points?
      ! Save the pressure-jumps at the beampoint, after multiplying with startup conditioner
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%pressure_new(n), n=0,ns-1,step)
      close (14) 
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%theta(n), n=0,ns-1,step)
      close (14) 
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_x'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%x(n), n=0,ns-1,step)
      close (14) 
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_y'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%y(n), n=0,ns-1,step)
      close (14) 
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vx'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%vx(n), n=0,ns-1,step)
      close (14) 
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vy'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns1) time, (beams(i)%vy(n), n=0,ns-1,step)
      close (14) 
      
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data'//beamstr, status = 'unknown', access = 'append')
      write (14, '(21(es15.8,1x))') &
        time,&
        beams(i)%pressure_new(ns-1),&
        beams(i)%pressure_new(3*ns/4),&
        beams(i)%pressure_new(ns/2),&
        beams(i)%x(ns-1),&
        beams(i)%y(ns-1),&
        beams(i)%vx(ns-1),&
        beams(i)%vy(ns-1),&
        beams(i)%theta(ns-1),&
        beams(i)%theta_dot(ns-1),&
        beams(i)%Force(1),&
        beams(i)%Force(2),&
        beams(i)%Force_unst(1),&
        beams(i)%Force_unst(2),&
        beams(i)%Force_press(1),&
        beams(i)%Force_press(2),&
        beams(i)%Inertial_Force(1),&
        beams(i)%Inertial_Force(2),&        
        beams(i)%E_kinetic,&
        beams(i)%E_elastic,&
        beams(i)%E_pot
      close (14)
    enddo
      
      
      
  !--------------------------------------------------------------
  ! save the end-point beam data ?
  !--------------------------------------------------------------
  elseif (iSaveBeam == 3) then !save only the last point
  
    do i=1, iBeam ! loop over beams
      ! for naming files..
      write (beamstr,'(i1)') i        
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data'//beamstr, status = 'unknown', access = 'append')
      write (14, '(21(es15.8,1x))') &
        time,&
        beams(i)%pressure_new(ns-1),&
        beams(i)%pressure_new(3*ns/4),&
        beams(i)%pressure_new(ns/2),&
        beams(i)%x(ns-1),&
        beams(i)%y(ns-1),&
        beams(i)%vx(ns-1),&
        beams(i)%vy(ns-1),&
        beams(i)%theta(ns-1),&
        beams(i)%theta_dot(ns-1),&
        beams(i)%Force(1),&
        beams(i)%Force(2),&
        beams(i)%Force_unst(1),&
        beams(i)%Force_unst(2),&
        beams(i)%Force_press(1),&
        beams(i)%Force_press(2),&
        beams(i)%Inertial_Force(1),&
        beams(i)%Inertial_Force(2),&        
        beams(i)%E_kinetic,&
        beams(i)%E_elastic,&
        beams(i)%E_pot
      close (14)
    enddo 
    
  endif
  
  
  !*********************
  ! the drag / lift is still saved 
  ! in cases where no beam is present. 
  ! however, as there will be no 
  !*********************
  
  !--------------------------------------------------------------
  ! save leading edge motions
  !--------------------------------------------------------------
  do i=1, iBeam
    write (beamstr,'(i1)') i 
    call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )
    
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'mouvement'//beamstr, status = 'unknown', access = 'append')
    write (14, '(10(es15.8,1x))') time, alpha, alpha_t, alpha_tt, (LeadingEdge(n), n=1,6)
    close (14)
    
  enddo


end subroutine SaveBeamData




subroutine SaveDeflectionLine ( time, beams )
  ! --------------------------------------------------------------------------------------------------------
  ! just save deflection line (complete beam, full resolution)
  ! --------------------------------------------------------------------------------------------------------
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time  
  character(len=22) :: format_ns
  character(len=1) :: beamstr
  type(solid), dimension(1:iBeam), intent(in) :: beams
  integer :: n,i  
  
  do i = 1, iBeam  
    !--for naming files..
    write (beamstr,'(i1)') i        
    write (format_ns, '("(",(i4.4),"(es12.5,1x))")') 2*ns+1   
    open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_complete'//beamstr, status = 'unknown', access = 'append')
    write (14, format_ns) time, (beams(i)%x(n), n=0,ns-1), (beams(i)%y(n), n=0,ns-1)
    close (14)
  enddo

end subroutine SaveDeflectionLine




subroutine SaveBeamPositions ( time, beams )
  ! --------------------------------------------------------------------------------------------------------
  ! this subroutine saves a complete deflection line, but it uses only 100 points (if there are enough)
  ! higher resolutions are not really required (and there may be buffer overflows in that case)
  !
  ! corrected 05 apr 2012 now working correctly
  !
  ! --------------------------------------------------------------------------------------------------------
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  type(solid), dimension(1:iBeam), intent(in) :: beams
  character(len=16) :: format_ns
  character(len=3)  :: ns_string
  character(len=1)  :: beamstr
  integer :: n,i,nsave=100, step, points
  
  ! loop over beams
  do i = 1, iBeam  
  write(beamstr,'(i1)') i

  if ((iSaveBeam==3).and.(ns>=100)) then
  !more than 100 points on the beam
      step=ns/nsave ! integer division
      points = ns/step !integer division
      
      write(ns_string, '(I3)') 2*points+1+4 !take a bit more place in the file (time+safety)
      format_ns  = '('//ns_string//'(es12.5,1x))'  
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_complete'//beamstr, status = 'unknown', access = 'append')
      write (14, format_ns) time, (beams(i)%x(n), n=0,ns-1,step), (beams(i)%y(n), n=0,ns-1,step)
      close (14)
  elseif ((iSaveBeam==3).and.(ns<100)) then
  !less than 100 points
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_complete'//beamstr, status = 'unknown', access = 'append')
      write (14, '(201(es12.5,1x))') time, (beams(i)%x(n), n=0,ns-1), (beams(i)%y(n), n=0,ns-1)
      close (14)
  endif

  enddo
end subroutine SaveBeamPositions








subroutine CreateHeader()
  use share_vars
  !----------------------------------------------------------------
  ! This routine reads in the PARAMS file and stores it in a
  ! global array. This array can now be written in every
  ! output file, for disambiguation.
  !----------------------------------------------------------------
  implicit none
  integer, parameter :: Nchar = 256 ! remember to change in share_vars, too, if required.
  integer :: io_error=0, linenumber=0,i=0
  character (len=Nchar) :: dummy
  !-------------------------------------------------------------------------------------------------------------
  ! first, lets determine how many entrys the PARAMS file has
      open (14, file = 'PARAMS.m', status='old', action='read')
      do while (io_error==0)
        read (14,'(A)',iostat=io_error) dummy
        linenumber=linenumber+1
      enddo
      linenumber=linenumber-1 !counted one too far
      close (14)
  !-------------------------------------------------------------------------------------------------------------
  ! then allocate the table for the header and fill the array
      open (14, file = 'PARAMS.m', status='old', action='read')
      allocate (Params_Header(1:linenumber+2))
      Params_Header = ""
      Params_Header(1)            ="% ======================================PARAMS begin======================================================="
      Params_Header(linenumber+2) ="% ======================================PARAMS END========================================================="
      do i=2,linenumber+1
        read (14,'(A)',iostat=io_error) dummy
        Params_Header(i) = "% "//dummy
      enddo
      close (14)
  !-------------------------------------------------------------------------------------------------------------
end subroutine CreateHeader




!================================================================================================================================================




subroutine WriteHeaderToFile(file_identifier)
  !----------------------------------------------------------------
  ! This routine writes the PARAMS file as a backup in the specified file_identifier
  !----------------------------------------------------------------
  use share_vars
  implicit none
  integer, intent (in) :: file_identifier
  integer :: i
  !------------------------------------------------------------------
  do i=1,size(Params_Header)
    write(file_identifier,'(A)') Params_Header(i)
  enddo
end subroutine WriteHeaderToFile