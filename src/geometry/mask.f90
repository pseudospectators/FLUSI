! Wrapper for different (possibly time-dependend) mask functions
subroutine create_mask(time, beam)
  use mpi
  use vars
  use solid_model
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr) :: eps_inv
  type(solid), intent(in) :: beam

  ! Attention: mask is reset here (and not in subroutines)
  mask = 0.d0

  ! Actual mask functions:
  select case(method)
  case("fsi")
    call create_mask_fsi(time, beam)
  case("mhd")
    call create_mask_mhd()
  end select

  ! Attention: division by eps is done here, not in subroutines.
  eps_inv = 1.d0/eps
  mask = mask*eps_inv  

end subroutine create_mask


! Wrapper to set imposed velocity
subroutine update_us(ub)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  ! do not create any mask when not using penalization
  if (iPenalization==1) then  
    select case(method)
    case("fsi")
        call update_us_fsi(ub)
    case("mhd")
        call update_us_mhd()
    case default    
        if(mpirank == 0) then
          write (*,*) "Error: unkown method in update_us; stopping."
          call suicide()
        endif
    end select  
  endif
end subroutine update_us






! This subroutine returns the value f of a smooth step function
! The sharp step function would be 1 if x<=t and 0 if x>t
! h is the semi-size of the smoothing area, so
! f is 1 if x<=t-h
! f is 0 if x>t+h
! f is variable (smooth) in between
subroutine smoothstep(f,x,t,h)
  use fsi_vars
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
      call suicide()
  end select
  
end subroutine smoothstep


! Set the penalization velocity for the given fields (f1,f2,f3) to the
! steady-state of the Taylor-Couette case.
subroutine taylor_couette_u_us(f1,f2,f3)
  use mpi
  use vars
  implicit none
  
  real(kind=pr),intent(inout)::f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz
  
  if(mpirank == 0) then
     if(r1 >= r2) then
        write (*,*) "r1 >= r2 is not allowed in Taylor-Coette flow; stopping."
        stop
     endif
  endif

  f1=0.d0
  f2=0.d0
  f3=0.d0
  
  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x + y*y)
        
        if(r <= R1) then
           do iz=ra(3),rb(3)
              ! Velocity field:
              f1(ix,iy,iz)=-omega1*y
              f2(ix,iy,iz)=omega1*x
              f3(ix,iy,iz)=0.d0
           enddo
        endif
        if(r >= R2) then
           do iz=ra(3),rb(3)
              ! NB: We assume that the outer wall is not moving.
              f1(ix,iy,iz)=0.d0
              f2(ix,iy,iz)=0.d0
              f3(ix,iy,iz)=0.d0
           enddo
        endif

     enddo
  enddo


end subroutine taylor_couette_u_us


