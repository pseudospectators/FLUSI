subroutine init_beams ( beams )
  !---------------------------------------------------
  ! initializes an array of beams. the initial state is always
  ! straight lines, possible oriented with different angles, at rest.
  !---------------------------------------------------
  use share_vars
  use motion
  implicit none
  integer :: n, i
  type(solid), dimension (1:iBeam), intent (out) :: beams  
  real (kind=pr) :: alpha, alpha_t, alpha_tt, dx, dy  
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  dx = xl/real(nx)
  dy = yl/real(ny)
  
  write (*,'("*** initializing ",i1," beams")')  iBeam
  write (*,'("--- Beam width (2*t_beam) covers ",(f4.1)," points")') 2.0*t_beam/max(dy,dx)
  

  
  !-------------------------------------------
  ! allocate beam storage for each beam
  !-------------------------------------------
  do i = 1, iBeam
    allocate ( beams(i)%x(0:ns-1) )
    allocate ( beams(i)%y(0:ns-1) )
    allocate ( beams(i)%vx(0:ns-1) )
    allocate ( beams(i)%vy(0:ns-1) )
    allocate ( beams(i)%ax(0:ns-1) )
    allocate ( beams(i)%ay(0:ns-1) )
    allocate ( beams(i)%theta(0:ns-1) )
    allocate ( beams(i)%theta_dot(0:ns-1) )
    allocate ( beams(i)%pressure_old(0:ns-1) )
    allocate ( beams(i)%pressure_new(0:ns-1) )
    allocate ( beams(i)%tau_old(0:ns-1) )
    allocate ( beams(i)%tau_new(0:ns-1) )    
    allocate ( beams(i)%beam_oldold(0:ns-1,1:6) )    
  enddo 
  
  !---------------------------------------------
  ! define adjustable parameters for each beam 
  ! this is position and motion protocoll
  !--------------------------------------------
  beams(1)%x0 = x0
  beams(1)%y0 = y0
  beams(1)%iMouvement = iMotion
  beams(1)%AngleBeam = AngleBeam
  beams(1)%phase = 0.0
  
  if (iBeam>1) then
  beams(2)%x0 = x0 + 1.25
  beams(2)%y0 = y0
  beams(2)%iMouvement = iMotion
  beams(2)%AngleBeam = AngleBeam
  select case ( beams(2)%iMouvement )
    case (300)
      beams(2)%phase = 0.00*pi ! in-phase case
    case (301)
      beams(2)%phase = 0.25*pi
    case (302)
      beams(2)%phase = 0.50*pi
    case (303)
      beams(2)%phase = 0.75*pi
    case (304)
      beams(2)%phase = 1.00*pi ! out-of-phase        
    case (305)
      beams(2)%phase = 1.25*pi
    case (306)
      beams(2)%phase = 1.50*pi
    case (307)
      beams(2)%phase = 1.75*pi      
    case default
      write (*,*) "error...iMotion bad bad"
      stop
  end select
  write (*,'("*** (dragonfly) phase=",es12.4)') beams(2)%phase/pi
  endif
  
  
  !---------------------------------------------
  ! loop over beams and initialize them
  !---------------------------------------------
  do i = 1, iBeam
    ! fetch leading edge position 
    call mouvement(0.0, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )
    ! initialize as zero
    beams(i)%x = 0.0
    beams(i)%y = 0.0
    beams(i)%vx = 0.0
    beams(i)%vy = 0.0
    beams(i)%ax = 0.0
    beams(i)%ay = 0.0
    beams(i)%theta = 0.0
    beams(i)%theta_dot = 0.0
    beams(i)%pressure_old = 0.0
    beams(i)%pressure_new = 0.0
    beams(i)%tau_old = 0.0
    beams(i)%tau_new = 0.0   
    beams(i)%beam_oldold = 0.0
    beams(i)%Inertial_Force = 0.0
    ! this is used for unst correction computation:
    beams(i)%drag_unst_new = 0.0
    beams(i)%drag_unst_old = 0.0
    beams(i)%lift_unst_new = 0.0
    beams(i)%lift_unst_old = 0.0
    beams(i)%UnsteadyCorrectionsReady = .false.
    beams(i)%StartupStep = .true.
    beams(i)%dt_old = 0.0
    
    ! remember that what follows is described in relative coordinates 
    do n = 0, ns-1 
      beams(i)%x(n) = LeadingEdge(1) + real(n) *ds
      beams(i)%y(n) = LeadingEdge(2)
      beams(i)%vx(n)= 0.0
      beams(i)%vy(n)= 0.0
      beams(i)%ax(n)= 0.0
      beams(i)%ay(n)= 0.0
    enddo
    call integrate_position (0.0, beams(i)) ! to take static angle into account
  enddo
  
end subroutine init_beams
