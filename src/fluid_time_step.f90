! Wrapper for different time marching methods
! FIXME: add documentation: which arguments are used for what?  What
! are their dimensions?
!-------------------------------------------------------------------------------
! INPUT
! time, dt0, dt1: time, dt0=t^n - t^n-1, dt1= t^n+1-t^n
! u             work array (input value is overwritten)
! uk            velocity in F-space at time level n
! vort          work array (input value is overwritten)
! expvis        integrating factor(s). if dt did not change, the input value is 
!               used to advance in time. otherwise, the exponential term is 
!               recomputed with the new time step
! work          work array 
! nlk           right hand side work array. holds the RHS at time level (n-1)
!               which we need for AB2. The second part of it is the work array 
!               used to store the new RHS at (t)
!-------------------------------------------------------------------------------
! OUTPUT
! u             velocity in phys space at time level (n) (th OLD level)
! uk            velocity in fourier space at the new time level (n+1)
! nlk           contains the RHS at time (n) and (n-1). in the next time step,
!               one can overwrite (n-1) and (n) becomes (n-1)
! work          pressure in phys space, if the flag use_solid_model==yes, and 
!               the penalized vorticity when using sponge==yes. Otherwise, not
!               used.
! expvis        the integrating factor(s) which are updated if the dt changes
! vort          the vorticity at time level (n) in phys space
!-------------------------------------------------------------------------------
subroutine FluidTimestep(time,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis,it)
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,n1,it
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real (kind=pr) :: t1

  t1=MPI_wtime()  
  ! Note that in the new version, dealiasing is done in cal_vis.

  ! Call fluid advancement subroutines.
  select case(iTimeMethodFluid)
  case("RK2")
     call RungeKutta2(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  case("AB2")
!      if(it == 0) then
     if(time == 0.d0) then
        call euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort,work,expvis)
     else
        call adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis)
     end if
  case("Euler")
     call euler(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  case default
     if (mpirank == 0) write(*,*) "Error! iTimeMethodFluid unknown. Abort."
     stop
  end select

  ! Force zero mode for mean flow
  if(method == "fsi") call set_mean_flow(uk,time)

  ! Set the divergence of the magnetic field to zero to avoid drift.
  if(method == "mhd") call div_field_nul(uk(:,:,:,4),uk(:,:,:,5),uk(:,:,:,6))

  time_fluid=time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep


! FIXME: add documentation: which arguments are used for what?
subroutine rungekutta2(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: it
  complex (kind=pr),intent (inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent (inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  integer :: i,j,l

  ! Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,it,nlk(:,:,:,:,0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! multiply the RHS with the viscosity
  do j=1,nf
     do i=1,3
        l=i+3*(j-1)
        nlk(:,:,:,l,0)=nlk(:,:,:,l,0)*expvis(:,:,:,j)
     enddo
  enddo

  ! Compute integrating factor, only done if necessary (i.e. time step
  ! has changed)
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  !-- Do the actual euler step. note nlk is already multiplied by vis
  do j=1,nf
     do i=1,3
        l=i+3*(j-1)
        uk(:,:,:,l)=(uk(:,:,:,l)*expvis(:,:,:,j) + dt1*nlk(:,:,:,l,0))
     enddo
  enddo

  ! RHS using the euler velocity
  call cal_nlk(time,it,nlk(:,:,:,:,1),uk,u,vort,work ) 
  call adjust_dt(dt1,u)

  ! do the actual time step. note the minus sign!!
  ! in the original formulation, it reads 
  ! u^n+1=u^n + dt/2*( N(u^n)*vis + N(u_euler) )
  ! but we don't want to save u_euler seperately, we want to overwrite
  ! u^n with it!  so the formulation reads
  ! u^n+1=u_euler - dt*N(u^n)*vis + dt/2*( N(u^n)*vis + N(u_euler) )
  !-- which yields simply
  !-- u^n+1=u_euler + dt/2*( -N(u^n)*vis + N(u_euler) )
  do i=1,nd
     uk(:,:,:,i)=uk(:,:,:,i) +0.5*dt1*(-nlk(:,:,:,i,0) + nlk(:,:,:,i,1) )
  enddo
end subroutine rungekutta2


! This is standard Euler-explicit time marching. It does not serve as
! startup scheme for AB2.
! FIXME: add documentation: which arguments are used for what?
subroutine euler(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: it
  complex (kind=pr),intent(inout):: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout):: &
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr),intent(inout) :: work (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  integer :: i,j,l

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,1),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  ! Multiply be integrating factor (always!)
  do j=1,nf
     do i=1,3
        l=i+3*(j-1)
        uk(:,:,:,l)=(uk(:,:,:,l) + dt1*nlk(:,:,:,l,1))*expvis(:,:,:,j)
     enddo
  enddo
end subroutine euler


! Note this is not an optimized Euler. It only does things we need for AB2.
! FIXME: add documentation: which arguments are used for what?
subroutine euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort,work,expvis)
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,it
  complex (kind=pr),intent(inout) ::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  integer :: i,j,l

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,n0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  ! Multiply be integrating factor (always!)
  do j=1,nf
     do i=1,3
        l=i+3*(j-1)
        uk(:,:,:,l)=(uk(:,:,:,l) + dt1*nlk (:,:,:,l,n0))*expvis(:,:,:,j)
        if (it==0) then
        nlk(:,:,:,l,n0)=nlk (:,:,:,l,n0)*expvis(:,:,:,j)
        endif
     enddo
  enddo

  if (mpirank ==0) write(*,'(A)') "*** info: did startup euler............"
end subroutine euler_startup


! FIXME: add documentation: which arguments are used for what?
subroutine adamsbashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,n1,it
  complex (kind=pr),intent(inout) ::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real (kind=pr) :: b10,b11
  integer :: i,j,a

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,n0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Calculate velocity at new time step 
  ! (2nd order Adams-Bashforth with exact integration of diffusion term)
  b10=dt1/dt0*(0.5*dt1 + dt0)
  b11=-0.5*dt1*dt1/dt0

  ! compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  ! Multiply be integrating factor (always!) 
  do j=1,nf
  do i=1,3
    a=i+3*(j-1)
    uk(:,:,:,a)=(uk(:,:,:,a)+b10*nlk(:,:,:,a,n0)+b11*nlk(:,:,:,a,n1))*expvis(:,:,:,j)
    if (it==0) then
    nlk(:,:,:,a,n0)=nlk(:,:,:,a,n0)*expvis(:,:,:,j)
    endif
  enddo
  enddo
end subroutine adamsbashforth


!-------------------------------------------------------------------------------
! Set the time step based on the CFL condition and penalization
! stability contidion. If dt_fixed is set, 
subroutine adjust_dt(dt1,u)
  use vars
  use mpi
  implicit none

  real (kind=pr), intent(in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  integer :: mpicode
  real (kind=pr), intent (out) :: dt1
  real(kind=pr) :: umax

!  if (dt_fixed>0.0) then
!     dt1=dt_fixed
!  else

     ! Determine the maximum velocity/magnetic field value, divided by
     ! grid size to determine the CFL condition.
     ! call maxabs1(umax,u) ! Max in each direction
     call maxabs(umax,u) ! Max magnitude

     !--Adjust time step at 0th process
     if(mpirank == 0) then
        if(.NOT.(umax.eq.umax)) then
           write(*,*) "Evolved field contains a NAN: aborting run."
           stop
        endif
     
     ! Impose the CFL condition.
        if (umax >= 1.0d-8) then
           dt1=cfl/umax
        else
           dt1=1.0d-2
        endif

        ! Round the time-step to one digit to reduce calls to cal_vis
        call truncate(dt1,dt1) 

        ! Impose penalty stability condition: dt cannot be less than 1/eps/
        if (iPenalization > 0) dt1=min(0.99*eps,dt1) 
        ! time step is smaller than eps 
        
        ! Don't jump past save-points: if the time-step is larger than
        ! the time interval between outputs, decrease the time-step.
!         if(tsave > 0.d0 .and. dt1 > tsave) then
!            dt1=tsave
!         endif
!         if(tintegral > 0.d0 .and. dt1 > tintegral) then
!            dt1=tintegral
!         endif
     endif
if (dt_fixed>0.d0) dt1=min(dt1,dt_fixed)
     ! Broadcast time step to all processes
     call MPI_BCAST(dt1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
!  endif
  
end subroutine adjust_dt


! FIXME: add documentation
subroutine truncate(a,b)
  ! rounds time step (from 1.246262e-2 to 1.2e-2)
  use vars
  implicit none

  real(kind=pr) :: a,b
  character (len=7) :: str

  write (str,'(es7.1)') a
  read (str,*) b
end subroutine truncate


! Force zero mode for mean flow
subroutine set_mean_flow(uk,time)
  use mpi
  use fsi_vars
  implicit none
  
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent (inout) :: time

  if(iMeanFlow == 1) then
     ! Force zero mode for mean flow
     ! TODO: this might not always select the proper mode; it could be
     ! better to determine if 0 is between ca(i) and cb(i) for i=1,2,3
     if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
        uk(0,0,0,1)=Uxmean
        uk(0,0,0,2)=Uymean
        uk(0,0,0,3)=Uzmean
     endif
  endif
end subroutine set_mean_flow


! Set umax to be the maximum of the magnitude of the velocity divided
! by the grid spacing.
subroutine maxabs(umax,ub)
  use vars
  use mpi
  implicit none

  real (kind=pr), intent(in) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(out) :: umax
  real(kind=pr),dimension(nf) :: uloc
  integer i,j
  integer :: mpicode
  
  ! Find the max velocity for each field.
  do i=1,nf
     j=(i-1)*3
     uloc(i)=maxval(&
          ub(:,:,:,1+j)*ub(:,:,:,1+j)/(dx*dx)+&
          ub(:,:,:,2+j)*ub(:,:,:,2+j)/(dy*dy)+&
          ub(:,:,:,3+j)*ub(:,:,:,3+j)/(dz*dz)&
          )
     uloc(i)=dsqrt(uloc(i))
  enddo

  ! Make uloc(1) the max local velocity for all fields.
  do i=2,nf
     if(uloc(1) < uloc(i)) uloc(1)=uloc(i) 
  enddo

  ! Set umax to be the maximum for all fields over all procs.
  call MPI_REDUCE(uloc(1),umax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpicode)
end subroutine maxabs

! Set umax to be the max velocity in each direction divided by the
! grid spacing.
subroutine maxabs1(umax,ub)
  use vars
  use mpi
  implicit none

  real (kind=pr), intent(in) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(out) :: umax
  real(kind=pr),dimension(nd) :: u_loc,u_loc_red
  integer :: i,j  
  integer :: mpicode
  
  do j=0,nf-1
     u_loc(1+3*j)=maxval(abs(ub(:,:,:,1+3*j)))/dx
     u_loc(2+3*j)=maxval(abs(ub(:,:,:,2+3*j)))/dy
     u_loc(3+3*j)=maxval(abs(ub(:,:,:,3+3*j)))/dz
  enddo

  call MPI_REDUCE(u_loc,u_loc_red,nd,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
       MPI_COMM_WORLD,mpicode)
  
  umax=0.d0
  do i=1,nd
     umax=max(umax,u_loc(i))
  enddo
end subroutine maxabs1

