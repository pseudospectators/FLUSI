! Wrapper for different (possibly time-dependend) mask functions
subroutine create_mask(time,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  real(kind=pr), intent(in) :: time  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nbeams), intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  
  real(kind=pr) :: eps_inv
  real(kind=pr) :: t1
  t1 = MPI_wtime() 
  
  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  !-------------------------------------------------------------
  ! create obstacle mask
  !-------------------------------------------------------------  
  ! do not create any mask when not using penalization
  if (iPenalization==1) then
    ! Actual mask functions:
    select case (iMask)
    case ("sphere","Sphere")    
      call Draw_Sphere(mask, mask_color, us)
    case ("cylinder","cylinder_x")    
      call Draw_cylinder_x(mask, mask_color, us)
    case ("romain_open_cavity")    
      call romain_open_cavity(mask, mask_color, us)
    case ("Flapper")    
      call Flapper (time, mask, mask_color, us)    
    case ("Insect")
      call Draw_Insect (time, Insect, mask, mask_color, us)
    case("Flexibility")      
      call Draw_flexible_plate(time, beams(1), mask, mask_color, us)
    case ("plate","Plate")
      call Draw_Plate (time, mask, mask_color, us ) ! 2d plate, etc (Dmitry, 25 Oct 2013)
    case ("noncircular_cylinder")
      call noncircular_cylinder( mask, mask_color, us )
    case ("couette")
      call taylor_couette(mask, mask_color, us)
    case("none")
      mask = 0.d0
    case default    
      write (*,*) "iMask="//iMask//" not properly set; stopping."
      call abort()
    end select
  endif

  !-------------------------------------------------------------
  ! add cavity / channel mask 
  !-------------------------------------------------------------  
  ! if desired, add cavity mask surrounding the domain
  if ((iCavity/="no").and.(iPenalization==1)) then
    call Add_Cavity (mask, mask_color, us)
  endif
  
  ! if desired, add channel mask 
  if ((iChannel/="no").and.(iPenalization==1)) then
    call Add_Channel (mask, mask_color, us)
  endif

  ! Attention: division by eps is done here, not in subroutines.
  eps_inv = 1.d0/eps
  mask = mask*eps_inv  

  ! -- for global timing.
  time_mask = time_mask + MPI_wtime() - t1
end subroutine create_mask


!-------------------------------------------------------------------------------
! This subroutine returns the value f of a smooth step function
! The sharp step function would be 1 if x<=t and 0 if x>t
! h is the semi-size of the smoothing area, so
! f is 1 if x<=t-h
! f is 0 if x>t+h
! f is variable (smooth) in between
!-------------------------------------------------------------------------------
subroutine smoothstep(f,x,t,h)
  use vars
  implicit none
  real (kind=pr), intent (out) :: f
  real (kind=pr), intent (in)  :: x,t,h

  !-------------------------------------------------
  ! cos shaped smoothing (compact in phys.space)
  !-------------------------------------------------
  if (x<=t-h) then
    f = 1.d0
  elseif (((t-h)<x).and.(x<(t+h))) then
    f = 0.5d0*(1.d0+dcos((x-t+h)*pi/(2.d0*h)) )
  else
    f = 0.d0
  endif
  
end subroutine smoothstep