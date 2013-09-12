! Wrapper for different (possibly time-dependend) mask functions
subroutine create_mask(time)
  use mpi_header
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr) :: eps_inv

  ! do not create any mask when not using penalization
  if (iPenalization==1) then  
    ! Attention: mask is reset here (and not in subroutines)
    mask = 0.d0

    ! Actual mask functions:
    select case(method)
    case("fsi")
        call Create_Mask_fsi(time)
    case("mhd")
        call Create_Mask_mhd()
    case default    
        if(mpirank == 0) then
          write (*,*) "Error: unkown method in create_mask; stopping."
          stop
        endif
    end select

    ! Attention: division by eps is done here, not in subroutines.
    eps_inv=1.d0/eps
    mask=mask*eps_inv  
  endif
end subroutine create_mask


! Wrapper to set imposed velocity
subroutine update_us(ub)
  use mpi_header
  use vars
  implicit none

  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  ! do not create any mask when not using penalization
  if (iPenalization==1) then  
    select case(method)
    case("fsi")
        call update_us_fsi(ub)
    case("mhd")
        call update_us_mhd(ub)
    case default    
        if(mpirank == 0) then
          write (*,*) "Error: unkown method in update_us; stopping."
          stop
        endif
    end select  
  endif
end subroutine update_us


! Spherical obstacle
subroutine draw_sphere 
  use mpi_header
  use fsi_vars
  implicit none

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, tmp, R, N_smooth

  N_smooth = 2.d0

  do ix=ra(1),rb(1)
     do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
           x=dble(ix)*dx
           y=dble(iy)*dy
           z=dble(iz)*dz
           R = dsqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
           if ( R <= 0.5d0*length+2.d0*N_smooth*max(dx,dy,dz) ) then
              call SmoothStep (tmp, R, 0.5d0*length , N_smooth*max(dx,dy,dz))
              mask (ix, iy, iz) = tmp
           endif

           if ((ix==0).and.(iy==0).and.(iz==0)) then
              mask (ix,iy,iz)=1.d0
           endif
           if ((ix==nx-1).and.(iy==ny-1).and.(iz==nz-1)) then
              mask (ix,iy,iz)=1.d0
           endif
        enddo
     enddo
  enddo
end subroutine draw_sphere


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
  real (kind=pr) :: delta, GradientERF!, SmoothStep
!   real (kind=pr) :: a,b,c,d
  !--polynomial coefficients:
  !   a =  1.0 / (4.0*(h**3))
  !   b = -3.0*t / (4.0*(h**3))
  !   c =  3.0*(t+h)*(t-h)/(4.0*(h**3))
  !   d =  ((t+h)**2)*(2.0*h-t)/(4.0*(h**3))

  !  if (x<=t-h) then
  !    f = 1.0
  !  elseif (((t-h)<x).and.(x<(t+h))) then
  ! !     f = a*(x**3) + b*(x**2) + c*x + d
  ! !     f = 1.0 - (x-t+h)/(2.0*h)
  !    f = 0.5*(1.+cos((x-t+h)*pi/(2.0*h)) )
  !  else
  !    f = 0.0
  !  endif

  !-----------------------------------
  ! version 14 - error function as non-oscilatory shape
  !-----------------------------------
  ! h - delta (gradient thickness)
  ! t - thickness (radius)
  GradientERF = abs( ( exp(-(2.0*1.0)**2)  - 1.0 )/sqrt(pi) )
  delta = h*GradientERF
  f = 0.5*( erf( (t-x)/delta ) + erf( (x+t)/delta )  )
end subroutine smoothstep


! Set the penalization velocity for the given fields (f1,f2,g3) to the
! steady-state of the Taylor-Couette case.
subroutine taylor_couette_u_us(f1,f2,f3)
use mpi_header
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


