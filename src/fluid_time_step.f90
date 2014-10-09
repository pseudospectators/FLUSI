!-------------------------------------------------------------------------------
! Warpper calling different individual time steppers
!-------------------------------------------------------------------------------
subroutine FluidTimestep(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use p3dfft_wrapper
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  
  real(kind=pr)::t1
  t1=MPI_wtime()

  !-----------------------------------------------------------------------------
  ! Call fluid advancement subroutines.
  !-----------------------------------------------------------------------------
  select case(iTimeMethodFluid)
  case("RK2")
      call RK2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("RK4")
      call RK4(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("AB2")
      if(time%it == 0) then
        call EE1(time,u,nlk,work,mask,mask_color,us,Insect,beams)
      else
        call AB2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
      endif
!   case ("AB2_rigid_solid")
!       call AB2_rigid_solid(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,workc,&
!            expvis,press,0,Insect)      
!   case("Euler")
!       call euler(time,it,dt0,dt1,u,uk,nlk,vort,work,workc,expvis,press)
!   case("FSI_AB2_iteration")
!       call FSI_AB2_iteration(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work, &
!            workc,expvis,press,beams)
!   case("FSI_AB2_staggered")
!       call FSI_AB2_staggered(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work, &
!            workc,expvis,press,beams)
  case("FSI_RK2_semiimplicit")
      call FSI_RK2_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case default
      if (root) write(*,*) "Error! iTimeMethodFluid unknown. Abort."
      call abort()
  end select
  
  
  !-----------------------------------------------------------------------------
  ! compute unsteady corrections in every time step.
  !-----------------------------------------------------------------------------
  if (unst_corrections==1) then
    call cal_unst_corrections ( time, mask, mask_color, us, Insect )  
  endif

  time_fluid=time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep





! !-------------------------------------------------------------------------------
! ! FSI scheme based on AB2/EE1 for the fluid, iterates coupling conditions.
! ! adapted from the 2D codes (V12), based on the PhD thesis of von Scheven
! !-------------------------------------------------------------------------------
! subroutine FSI_AB2_iteration(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,&
!            workc,expvis,press,beams)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   use solid_model
!   use insect_module
!   implicit none
! 
!   real(kind=pr),intent(inout) :: time,dt1,dt0
!   integer,intent (in) :: n0,n1,it
!   complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
!   complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
!   complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
!   real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
!   real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
!   real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
!   type(solid),dimension(1:nbeams),intent(inout) :: beams
!   
!   !-- iteration specific variables 
!   type(solid), dimension(1) :: beams_old
!   real(kind=pr),dimension(0:ns-1) :: deltap_new, deltap_old, bpress_old_iterating
!   real(kind=pr)::bruch, upsilon_new, upsilon_old, kappa2, ROC1,ROC2, norm
!   real(kind=pr)::omega_old, omega_new
!   type(diptera)::Insect_dummy
!   integer :: inter
!   logical :: iterate
!   
!   ! useful error messages
!   if (use_solid_model/="yes") then
!     write(*,*) "using FSI_AB2_iteration without solid model?"
!     call abort()
!   endif
!   
!   ! allocate extra space for velocity in Fourier space
!   if (.not.allocated(uk_old)) then
!     call alloccomplexnd(uk_old)    
!   endif
!   ! we need that to compute the pressure at preliminary times
!   if (.not.allocated(nlk_tmp)) then
!     allocate(nlk_tmp(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
!   endif
!   ! copy velocity at time level (n)
!   uk_old = uk
!   
!   ! initialize iteration variables
!   inter = 0
!   iterate = .true.
!   deltap_new = 0.d0
!   deltap_old = 0.d0
!   upsilon_new = 0.d0
!   upsilon_old = 0.d0
!   beams_old = beams
!   omega_old = 0.5d0
!   omega_new = 0.5d0
!   ! the passive scalar is advanced only once during the iteration process to keep the cost down.
!   if(use_passive_scalar==1) compute_scalar = .true.
!   
!   ! predictor for the pressure     
!   beams(1)%pressure_new = beams(1)%pressure_old
!   
!   ! begin main iteration loop
!   do while (iterate)     
!     !---------------------------------------------------------------------------
!     ! create mask
!     !---------------------------------------------------------------------------
!     call create_mask(time, Insect_dummy, beams)
!     
!     !---------------------------------------------------------------------------
!     ! advance fluid to from (n) to (n+1)
!     !---------------------------------------------------------------------------
!     if(use_passive_scalar==1) then 
!       ! advance the scalar only once when iterating.
!       if (inter==0) compute_scalar = .true.
!       if (inter/=0) compute_scalar = .false.
!     endif
!     
!     ! start from t(n) again, except for the scalar
!     uk(:,:,:,1:3) = uk_old(:,:,:,1:3)
!     
!     if(it == 0) then
!       call euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort, &
!            work,workc,expvis,press,inter)
!     else
!       call adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort, &
!            work,workc,expvis,press,inter)
!     endif
!     
!     !---------------------------------------------------------------------------
!     ! get forces at new time level
!     !---------------------------------------------------------------------------
!     bpress_old_iterating = beams(1)%pressure_new(0:ns-1) ! exit of the old iteration
!     call pressure_given_uk(time,u,uk,nlk_tmp,vort,work,workc,press)
!     call get_surface_pressure_jump (time, beams(1), press, timelevel="new")
!     
!     !---------------------------------------------------------------------------
!     ! relaxation
!     !---------------------------------------------------------------------------
!     ! whats the diff betw new interp press and last iteration's step?
!     deltap_new = bpress_old_iterating - beams(1)%pressure_new(0:ns-1)
!     if (inter==0) then
!       upsilon_new = 0.0d0
!       ! von scheven normalizes with the explicit scheme, which is what we do now
!       norm = sqrt(sum((beams(1)%pressure_new(0:ns-1)-beams(1)%pressure_old(0:ns-1))**2))    
!     else
!       bruch = (sum((deltap_old-deltap_new)*deltap_new)) &
!             / (sum((deltap_old-deltap_new)**2))
!       upsilon_new = upsilon_old + (upsilon_old-1.d0) * bruch
!     endif
!     kappa2 = 1.d0 - upsilon_new
!     ! new iteration pressure is old one plus star
!     beams(1)%pressure_new(0:ns-1) = (1.d0-kappa2)*bpress_old_iterating(0:ns-1) &
!         + kappa2*beams(1)%pressure_new(0:ns-1)
!                           
!     !---------------------------------------------------------------------------  
!     ! advance solid model from (n) to (n+1)
!     !---------------------------------------------------------------------------
!     beams_old(1)%pressure_new = beams(1)%pressure_new
!     beams = beams_old ! advance from timelevel n
!     call SolidSolverWrapper( time, dt1, beams )
!     
!     !---------------------------------------------------------------------------
!     ! convergence test
!     !---------------------------------------------------------------------------
!     ROC1 = dsqrt( sum((beams(1)%pressure_new(0:ns-1)-bpress_old_iterating)**2)) / dble(ns)
!     ROC2 = dsqrt( sum((beams(1)%pressure_new(0:ns-1)-bpress_old_iterating)**2)) / norm 
!     if (((ROC2<1.0e-3).or.(inter==100)).or.(it<2)) then
!       iterate = .false.
!     endif
!   
!     ! iterate
!     deltap_old = deltap_new
!     upsilon_old = upsilon_new
!     inter = inter + 1
!     omega_old=omega_new
!     
!     if (root) then
!       open (15, file='iterations_log.t',status='unknown',position='append')
!       write(15,'("t=",es12.4," dt=",es12.4," inter=",i3," ROC=",es15.8,&
!                 &" ROC2=",es15.8," p_end=",es15.8," kappa2=",es15.8)') &
!       time,dt1,inter,ROC1,ROC2,beams(1)%pressure_new(ns-1), kappa2
!       close (15)
!     endif
!   enddo
!   
!   ! dump iteration information to disk
!   if (root) then
!     open (15, file='iterations.t',status='unknown',position='append')
!     write(15,'(2(es15.8,1x),i3,2(es15.8,1x))') time, dt1, inter, ROC1, ROC2
!     close(15)
!   
!     ! mark end of time step
!     open (15, file='iterations_log.t',status='unknown',position='append')
!     write(15,'("---")')
!     close(15)
!   endif
!   
!   ! free work array
! !   deallocate (uk_old,nlk_tmp)
! end subroutine FSI_AB2_iteration



! !-------------------------------------------------------------------------------
! ! explicit FSi scheme, cheapest and simplest possible. Advances first the fluid
! ! then the solid, and computes the solid with the pressure from the old time 
! ! level (n), avoiding computing it at the new level, which saves about 50%
! ! with respect to FSI_AB2_semiimplicit. The latter may however be more accurat
! ! or stable.
! !-------------------------------------------------------------------------------
! subroutine FSI_AB2_staggered(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,&
!            workc,expvis,press,beams)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   use solid_model
!   use insect_module
!   implicit none
! 
!   real(kind=pr),intent(inout) :: time,dt1,dt0
!   integer,intent (in) :: n0,n1,it
!   complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
!   complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
!   complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
!   real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
!   real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
!   real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
!   type(solid),dimension(1:nbeams),intent(inout) :: beams
!   type(diptera)::Insect_dummy
!   
!   ! useful error messages
!   if (use_solid_model/="yes") then 
!     write(*,*) "using FSI_AB2_staggered without solid model?"
!     call abort()
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! create mask
!   !---------------------------------------------------------------------------
!   call create_mask(time,Insect_dummy, beams(1))
!   
!   !---------------------------------------------------------------------------
!   ! advance fluid to from (n) to (n+1)
!   !---------------------------------------------------------------------------
!   if(it == 0) then
!     call euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   else
!     call adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! get forces at old time level, since press is at t^n.
!   ! save to both beam%pressure_new and beam%pressure_old
!   !---------------------------------------------------------------------------
!   call get_surface_pressure_jump (time, beams(1), press)
!   
!   !---------------------------------------------------------------------------  
!   ! advance solid model from (n) to (n+1)
!   !---------------------------------------------------------------------------
!   call SolidSolverWrapper( time, dt1, beams )
!     
! end subroutine FSI_AB2_staggered


! !-------------------------------------------------------------------------------
! ! semi implicit explicit staggered scheme for FSI simulations, uses AB2 for the
! ! fluid (or euler on startup) and evaluates the pressure at both old and new 
! ! time level. since computing the pressure is almost as expensive as doing a full
! ! fluid time step, this scheme is twice as expensive as its explicit counterpart
! ! FSI_AB2_staggered.
! !-------------------------------------------------------------------------------
! subroutine FSI_AB2_semiimplicit(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,&
!            work,workc,expvis,press,beams)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   use solid_model
!   use insect_module
!   implicit none
! 
!   real(kind=pr),intent(inout) :: time,dt1,dt0
!   integer,intent (in) :: n0,n1,it
!   complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
!   complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
!   complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
!   real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
!   real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
!   real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
!   type(solid),dimension(1:nbeams),intent(inout) :: beams
!   type(diptera)::Insect_dummy
!   
!   ! useful error messages
!   if (use_solid_model/="yes") then 
!    write(*,*) "using FSI_AB2_staggered without solid model?"
!    call abort()
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! create mask
!   !---------------------------------------------------------------------------
!   call create_mask(time,Insect_dummy, beams(1))
!   
!   !---------------------------------------------------------------------------
!   ! advance fluid to from (n) to (n+1)
!   !---------------------------------------------------------------------------
!   if(it == 0) then
!     call euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   else
!     call adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! get forces at old/new time level
!   !---------------------------------------------------------------------------
!   ! TODO: do we need that? not for BDF I think
!   call get_surface_pressure_jump (time, beams(1), press, timelevel="old")
!   ! we can overwrite NLK(n1) since this is the one overwritten in the next call
!   call pressure_given_uk(time,u,uk,nlk(:,:,:,:,n1),vort,work,workc,press)
!   call get_surface_pressure_jump (time, beams(1), press, timelevel="new")
!   
!   !---------------------------------------------------------------------------  
!   ! advance solid model from (n) to (n+1)
!   !---------------------------------------------------------------------------
!   call SolidSolverWrapper( time, dt1, beams )
!     
! end subroutine FSI_AB2_semiimplicit


!-------------------------------------------------------------------------------
! semi implicit explicit staggered scheme for FSI simulations
! We first interpolate the pressure at the old timelevel t^n
! then advance the fluid to the next time level t^n+1 by using a suitable
! integrator. The interpolated pressure at the new time level is then the input
! to the solid solver which is advanced then
! during the fluid time step, the mask function is held constant
!-------------------------------------------------------------------------------
subroutine FSI_RK2_semiimplicit(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  ! useful error messages
  if (use_solid_model/="yes") then 
   write(*,*) "using FSI_RK2_semiimplicit without solid model?"
   call abort()
  endif
  
  !---------------------------------------------------------------------------
  ! get forces at old time level
  !---------------------------------------------------------------------------  
  ! since the BDF scheme evaluates the beam RHS only at the new time level (and
  ! not at the old one at all) the pressure at t^n is not required when using BDF2
  if ( TimeMethodSolid .ne. "BDF2") then
    call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="old")
  endif
  
  !---------------------------------------------------------------------------
  ! advance fluid to from (n) to (n+1)
  !---------------------------------------------------------------------------
  call RK2 (time,u,nlk,work,mask,mask_color,us,Insect,beams)
  
  !---------------------------------------------------------------------------
  ! get forces at new time level
  !---------------------------------------------------------------------------  
  call get_surface_pressure_jump (time%time, beams(1), u(:,:,:,4), timelevel="new")
  
  !---------------------------------------------------------------------------  
  ! advance solid model from (n) to (n+1)
  !---------------------------------------------------------------------------
  call SolidSolverWrapper ( time%time, time%dt_new, beams )
    
end subroutine FSI_RK2_semiimplicit


! !-------------------------------------------------------------------------------
! ! FIXME: add documentation: which arguments are used for what?
! !-------------------------------------------------------------------------------
! subroutine FSI_RK2(time,it,dt0,dt1,u,uk,nlk,vort,work,workc,expvis,press,beams)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   use solid_model
!   implicit none
! 
!   real(kind=pr),intent (inout) :: time,dt1,dt0
!   integer,intent (in) :: it
!   complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
!   complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
!   complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
!   real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
!   real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
!   real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
!   type(solid),dimension(1),intent(inout) :: beams
!   type(solid),dimension(1):: beams_old
!   integer :: i
!   ! Compute integrating factor, only done if necessary (i.e. time step
!   ! has changed)
!   if (dt1 .ne. dt0) call cal_vis(dt1,expvis)
!   
!   call create_mask(time, beams(1))
!   call cal_nlk(time,it,nlk(:,:,:,:,0),uk,u,vort,work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
!   call get_surface_pressure_jump (time, beams(1), work, timelevel="old")
!   call adjust_dt(dt1,u)
! 
!   beams_old = beams
!   
!   ! multiply the RHS with the viscosity
!   do i=1,3
!     nlk(:,:,:,i,0)=nlk(:,:,:,i,0)*expvis(:,:,:,1)
!   enddo
! 
!   !-- Do the actual euler step. note nlk is already multiplied by vis
!   do i=1,3
!     uk(:,:,:,i)=(uk(:,:,:,i)*expvis(:,:,:,1) + dt1*nlk(:,:,:,i,0))
!   enddo
! 
!   ! RHS using the euler velocity
!   call cal_nlk(time,it,nlk(:,:,:,:,1),uk,u,vort,work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))) 
!   ! the pressure is now at the intermediate time step
!   call get_surface_pressure_jump (time, beams(1), work, timelevel="new")
!   call SolidSolverWrapper( time, dt1, beams )
!   call create_mask(time, beams(1))
!   ! RHS with the modified beam
!   call cal_nlk(time,it,nlk(:,:,:,:,1),uk,u,vort,work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))) 
!   ! advance to new time level
!   do i=1,nd
!      uk(:,:,:,i)=uk(:,:,:,i) +0.5*dt1*(-nlk(:,:,:,i,0) + nlk(:,:,:,i,1) )
!   enddo
!   call pressure_given_uk(uk,work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
!   call get_surface_pressure_jump (time, beams(1), work, timelevel="new")
!   beams_old(1)%pressure_new=beams(1)%pressure_new
!   beams = beams_old
!   call SolidSolverWrapper( time, dt1, beams )
! end subroutine FSI_RK2


!-------------------------------------------------------------------------------
! Standard RK2 scheme for integrating the RHS defined in cal_nlk
! The scheme is also referred to as "midpoint rule"
! INPUT:
!       time: a derived datatype that contains current time, iteration number
!             time step an so on (definition see vars.f90)
!       u: The solution vector with "neq" components at time t
!       work: work array, currently unused
!       mask: the mask function (geometry of the problem)
!-------------------------------------------------------------------------------
subroutine RK2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  type(timetype) :: t
  t=time
    
  ! define the time step we want to use
  call adjust_dt(u,time%dt_new)
  
  ! now NLK2 is old velocity
  nlk(:,:,:,:,2) = u
  
  ! compute rhs at old time and make half a time step into the future
  call cal_nlk(time, u, nlk(:,:,:,:,1), work, mask, mask_color, us, Insect, beams)
  u = u + 0.5d0*time%dt_new * nlk(:,:,:,:,1)
  t%time = time%time + 0.5d0*time%dt_new
  
  
  ! now NLK1 is old velocity
  nlk(:,:,:,:,1) = nlk(:,:,:,:,2)
  
  ! using the solution at t^(n+1/2), we compute the RHs again
  call cal_nlk(t, u, nlk(:,:,:,:,2), work, mask, mask_color, us, Insect, beams)

  ! final step, advancing the midpoint derivative by a whole time step
  u = nlk(:,:,:,:,1) + time%dt_new * nlk(:,:,:,:,2)
end subroutine RK2

!-------------------------------------------------------------------------------

subroutine RK4(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  type(timetype) :: t
  t = time
  
  call adjust_dt(u,time%dt_new)
  
  ! NLK 5th register holds old velocity
  nlk(:,:,:,:,5) = u
  
  !-- Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,u,nlk(:,:,:,:,1),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + 0.5d0*time%dt_new*nlk(:,:,:,:,1)
  t%time = time%time + 0.5d0*time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,2),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + 0.5d0*time%dt_new*nlk(:,:,:,:,2)
  t%time = time%time + 0.5d0*time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,3),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + time%dt_new * nlk(:,:,:,:,3)
  t%time = time%time + time%dt_new
  call cal_nlk(t,u,nlk(:,:,:,:,4),work,mask,mask_color,us,Insect,beams)
  
  u = nlk(:,:,:,:,5) + time%dt_new/6.d0*(nlk(:,:,:,:,1)+2.d0*nlk(:,:,:,:,2)&
      +2.d0*nlk(:,:,:,:,3)+nlk(:,:,:,:,4))
  
end subroutine RK4

!-------------------------------------------------------------------------------

subroutine EE1(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 

  call adjust_dt(u,time%dt_new)
  
  !-- Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,u,nlk(:,:,:,:,1),work,mask,mask_color,us,Insect,beams)
  
  !-- advance in time
  u = u + time%dt_new * nlk(:,:,:,:,1)
  
end subroutine EE1

!-------------------------------------------------------------------------------

subroutine AB2(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  real(kind=pr)::b10,b11,dt1,dt0
  integer::i
  
  if (time%n0==time%n1 .or. time%n0<1 .or. time%n0>2 .or. time%n1<1 .or. time%n1>2) then
    if(mpirank==0) write(*,*) time%n0, time%n1
    if(mpirank==0) write(*,*) "In AB2, the registers n0,n1 are not well defned"
    call abort()
  endif
  
  call adjust_dt(u,time%dt_new)

  !-- Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,u,nlk(:,:,:,:,time%n0),work,mask,mask_color,us,Insect,beams)

  dt1 = time%dt_new
  dt0 = time%dt_old
  
  b10 = dt1/dt0*(0.5d0*dt1 + dt0)
  b11 = -0.5d0*dt1*dt1/dt0

  !-- Advance in time
  u = u + b10*nlk(:,:,:,:,time%n0) + b11*nlk(:,:,:,:,time%n1)
  
end subroutine AB2

!-------------------------------------------------------------------------------

! ! time stepper for rigid soldi FSI, for example insect takeoff or falling sphere
! subroutine AB2_rigid_solid(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,workc,&
!                            expvis,press,iter,Insect)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   use insect_module
!   implicit none
! 
!   real(kind=pr),intent(inout)::time,dt1,dt0
!   integer,intent(in)::n0,n1,it,iter
!   complex(kind=pr),intent(inout) ::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
!   complex(kind=pr),intent(inout)::&
!        nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
!   ! the workc array is not always allocated, ensure allocation before using
!   complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)        
!   real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
!   real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
!   real(kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
!   ! pressure array. this is with ghost points for interpolation
!   real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
!   real(kind=pr)::b10,b11,t1
!   integer::i
!   type(diptera),intent(inout)::Insect
! 
!   if ((SolidDyn%idynamics/=1).and.(mpirank==0)) then
!     write(*,*) "ERROR AB2_rigid_solid and flag SolidDyn%idynamics/=1"
!     write(*,*) "it makes no sense to do that, change iFluidTimeMethod=AB2"
!     call abort()
!   endif
!   
!   if ((method/="fsi").and.(mpirank==0)) then
!     write(*,*) "AB2_rigid_solid is an FSI method and not suitable for MHD"
!     call abort()
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! advance fluid to from (n) to (n+1)
!   !---------------------------------------------------------------------------
!   if(it == 0) then
!     call euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   else
!     call adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort, &
!          work,workc,expvis,press,0)
!   endif
!   
!   !---------------------------------------------------------------------------
!   ! compute hydrodynamic forces for rigid solid solver
!   !---------------------------------------------------------------------------
!   t1 = MPI_wtime()
!   if (unst_corrections ==1) then
!     call cal_unst_corrections ( time, dt0, Insect )
!   endif
!   call cal_drag ( time, u, Insect ) ! note u is OLD time level
!   time_drag = time_drag + MPI_wtime() - t1
!   
!   !---------------------------------------------------------------------------
!   ! solve Newton's second law
!   !---------------------------------------------------------------------------
!   ! note dt0 is OLD time step t(n)-t(n-1)
!   ! advance in time ODEs that describe rigid solids
!   call rigid_solid_time_step(time,dt0,dt1,it,Insect)
!   
! 
! end subroutine AB2_rigid_solid


!-------------------------------------------------------------------------------
! Set the time step based on the CFL condition and penalization
! stability contidion. The following limitations exist:
! 1 - CFL condition (for umax and the speed of sound, c_0)
! 2 - fixed time step dt_fixed, ignoring all other constraints, if set in params
! 3 - penalization restriction dt<eps
! 4 - maximum time step dt_max, if set in params
! 5 - dt is smaller than tsave and tintegral
!-------------------------------------------------------------------------------
subroutine adjust_dt(u,dt1)
  use vars
  use mpi
  use basic_operators
  implicit none

  real(kind=pr), intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer::mpicode
  real(kind=pr), intent(out)::dt1
  real(kind=pr)::umax

  if (dt_fixed>0.0) then
     !-- fix the time step no matter what. the result may be unstable.
     dt1=dt_fixed
  else
     !-- FSI runs just need to respect CFL for velocity
     umax = fieldmaxabs(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
     
     !-- Adjust time step at 0th process
     if(mpirank == 0) then
        if(is_nan(umax)) then
           write(*,*) "Evolved field contains a NAN: aborting run."
           call abort()
        endif
     
        !-- Impose the CFL condition.
        if (umax >= 1.0d-8) then
           dt1=min(dx,dy,dz)*cfl/umax
        else
           !-- umax is very very small
           dt1=1.0d-3
        endif

        !-- Impose penalty stability condition: dt cannot be larger than eps
        if (iPenalization > 0) dt1=min(0.99*eps,dt1) 
        
        ! Don't jump past save-points: if the time-step is larger than
        ! the time interval between outputs, decrease the time-step.
        if(tsave > 0.d0 .and. dt1 > tsave) then
           dt1=tsave
        endif
        if(tintegral > 0.d0 .and. dt1 > tintegral) then
           dt1=tintegral
        endif

        ! CFL condition for speed of sound
        dt1 = min( dt1, min(dx,dy,dz)*cfl/c_0 )
        
        !-- impose max dt, if specified
        if (dt_max>0.d0) dt1=min(dt1,dt_max)
     endif

     ! Broadcast time step to all processes
     call MPI_BCAST(dt1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
  endif
  
end subroutine adjust_dt