! Wrapper for different (possibly time-dependend) mask functions
subroutine create_mask(time,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  real(kind=pr), intent(in) :: time  
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  type(solid), dimension(1:nbeams), intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  
  real(kind=pr) :: eps_inv
  real(kind=pr) :: t1
  t1 = MPI_wtime() 
  
  ! Attention: mask is reset here (and not in subroutines)
  mask = 0.d0

  ! Actual mask functions:
  call create_mask_fsi(time,Insect,beams)

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
  real (kind=pr) :: delta, GradientERF

  select case (iSmoothing)
  case ("erf")
      !-------------------------------------------------
      ! error function as non-oscilatory shape
      !-------------------------------------------------
      ! h - delta (gradient thickness)
      ! t - thickness (radius)
      GradientERF = dabs( ( dexp(-(2.d0*1.d0)**2)  - 1.d0 )/dsqrt(pi) )
      delta = h*GradientERF
      f = 0.5d0*( erf( (t-x)/delta ) + erf( (x+t)/delta )  )
  case ("cos")
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
  case default
      !-------------------------------------------------
      write(*,*) "Smoothing parameter not rightly set", iSmoothing
      call abort()
  end select
  
end subroutine smoothstep