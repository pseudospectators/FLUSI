subroutine SaveBeamData( time, beams, dt1 )
  use fsi_vars
  implicit none
  real (kind=pr), intent (in) :: time,dt1
  type (solid), dimension(1:nBeams), intent (in) :: beams
  character(len=16) :: format_ns1
  character(len=3)  :: ns1_string
  character(len=1)  :: beamstr
  real (kind=pr) :: alpha, alpha_t, alpha_tt
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  integer :: n,step,i
  

   
  ! set up formats
  write(ns1_string, '(I3)') ns+1
  format_ns1 = '('//ns1_string//'(es12.5,1x))'
  
  
  do i=1, nBeams ! loop over beams
    ! for naming files..
    write (beamstr,'(i1)') i        
    open  (14, file = 'beam_data'//beamstr, status = 'unknown',position='append')
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
  ! save leading edge motions
  !--------------------------------------------------------------
  do i=1, nBeams
    write (beamstr,'(i1)') i 
    call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )
    
    open (14, file = 'mouvement'//beamstr, status = 'unknown',position='append')
    write (14, '(10(es15.8,1x))') time, alpha, alpha_t, alpha_tt, (LeadingEdge(n), n=1,6)
    close (14)
    
  enddo


end subroutine SaveBeamData
