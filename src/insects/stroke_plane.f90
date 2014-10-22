!-------------------------------------------------------------------------------
! Stroke plane
! Input: 
!       time 
! Output:
!       eta_stroke: stroke plane angle 
subroutine StrokePlane ( time, Insect )
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: eta_stroke

  select case (Insect%BodyMotion)
  case ("fixed")
    eta_stroke = deg2rad(0.d0)   
  case ("fixed45")
    eta_stroke = deg2rad(-45.d0)
  case ("wheeling")
    eta_stroke = deg2rad(0.d0)
  case ("hovering")
!    eta_stroke = deg2rad(-35.d0)
    eta_stroke = deg2rad(-45.d0)  ! Comparison with Maeda (Dmitry, 7 Nov 2013)
  case ("x0y0z0")
    eta_stroke = deg2rad(-45.d0)  ! Comparison with Maeda (Dmitry, 7 Nov 2013)
  case ("flapper")   ! Comparison with Dickinson et al. (Dmitry, 19 Nov 2013)
    eta_stroke = deg2rad(0.d0)
  case ("takeoff")
!    eta_stroke = deg2rad(-28.d0) ! 62-90, Fontaine et al., fig 13 (Dmitry, 14 Nov 2013)
    eta_stroke = Insect%eta_stroke ! read from file
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::StrokePlane: motion case (Insect%BodyMotion) undefined"
    call abort()
    endif
  end select
  
  ! save it in the insect
  Insect%eta_stroke = eta_stroke
  
end subroutine StrokePlane
