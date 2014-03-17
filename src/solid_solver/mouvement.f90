subroutine mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam)
  implicit none
  type(solid), intent(in) :: beam
  real (kind=pr), intent(out) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension(1:6), intent(out) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)


  select case (beam%iMouvement)
    case (0) ! fixed beam (no leading edge mouvement)    
      !--------------------------------------------------------------------------------------------------------
      ! fixed
      !--------------------------------------------------------------------------------------------------------    
      LeadingEdge = 0.0
      LeadingEdge(1) = beam%x0
      LeadingEdge(2) = beam%y0
      alpha    = beam%AngleBeam*pi/180.0
      alpha_t  = 0.0
      alpha_tt = 0.0
      
  end select

end subroutine