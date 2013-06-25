! Wrapper for different time marching methods
! FIXME: add documentation: which arguments are used for what?  What
! are their dimensions?
subroutine FluidTimestep(time,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis,it)
  use mpi_header
  use vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,n1,it
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real (kind=pr) :: t1

  t1=MPI_wtime()  
  ! Note that in the new version, dealiasing is done in cal_vis.

  ! Call fluid advancement subroutines.
  select case(iTimeMethodFluid)
  case("RK2")
     call RungeKutta2(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  case("AB2")
     if(it == 0) then
        call Euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort,work,expvis)
     else
        call AdamsBashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis)
     end if
  case("Euler")
     call Euler(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  case default
     if (mpirank == 0) write(*,*) "Error! iTimeMethodFluid unknown. Abort."
  end select

  ! Force zero mode for mean flow
  if(iMeanFlow == 1) call set_mean_flow(uk,time)

  time_fluid=time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep


! FIXME: add documentation: which arguments are used for what?
subroutine RungeKutta2(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  use mpi_header
  use fsi_vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: it
  complex (kind=pr),intent (inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent (inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: i

  ! Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  call cal_nlk(time,it,nlk(:,:,:,:,0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! multiply the RHS with the viscosity
  do i=1,3
     nlk(:,:,:,i,0)=nlk(:,:,:,i,0)*expvis
  enddo

  ! Compute integrating factor, only done if necessary (i.e. time step
  ! has changed)
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  !-- Do the actual euler step. note nlk is already multiplied by vis
  do i=1,3
     uk(:,:,:,i)=(uk(:,:,:,i)*expvis + dt1*nlk(:,:,:,i,0))
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
  do i=1,3
     uk(:,:,:,i)=uk(:,:,:,i) +0.5*dt1*(-nlk(:,:,:,i,0) + nlk(:,:,:,i,1) )
  enddo
end subroutine RungeKutta2


! This is standard Euler-explicit time marching. It does not serve as
! startup scheme for AB2.
! FIXME: add documentation: which arguments are used for what?
subroutine Euler(time,it,dt0,dt1,u,uk,nlk,vort,work,expvis)
  use mpi_header
  use fsi_vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: it
  complex (kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout):: &
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent(inout) :: work (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: i

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,1),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis (dt1,expvis)
  endif

  ! Multiply be integrating factor (always!)
  do i=1,3
     uk(:,:,:,i)=(uk(:,:,:,i) + dt1*nlk(:,:,:,i,1))*expvis
  enddo
end subroutine Euler


! Note this is not an optimized Euler. It only does things we need for AB2.
! FIXME: add documentation: which arguments are used for what?
subroutine Euler_startup(time,it,dt0,dt1,n0,u,uk,nlk,vort,work,expvis)
  use mpi_header
  use fsi_vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,it
  complex (kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout):: &
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: i

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,n0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis (dt1,expvis)
  endif

  ! Multiply be integrating factor (always!)
  do i=1,3
     uk(:,:,:,i)=(uk(:,:,:,i) + dt1*nlk (:,:,:,i,n0))*expvis
     nlk (:,:,:,i,n0)=nlk (:,:,:,i,n0)*expvis
  enddo

  if (mpirank ==0) write(*,'(A)') "*** info: did startup euler............"
end subroutine Euler_startup


! FIXME: add documentation: which arguments are used for what?
subroutine AdamsBashforth(time,it,dt0,dt1,n0,n1,u,uk,nlk,vort,work,expvis)
  use mpi_header
  use fsi_vars
  implicit none

  real (kind=pr),intent (inout) :: time,dt1,dt0
  integer,intent (in) :: n0,n1,it
  complex (kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent(inout) :: expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real (kind=pr) :: b10,b11
  integer :: i

  ! Calculate fourier coeffs of nonlinear rhs and forcing
  call cal_nlk(time,it,nlk(:,:,:,:,n0),uk,u,vort,work)
  call adjust_dt(dt1,u)

  ! Calculate velocity at new time step 
  ! (2nd order Adams-Bashforth with exact integration of diffusion term)
  b10=dt1/dt0*(0.5*dt1 + dt0)
  b11=-0.5*dt1**2/dt0

  ! compute integrating factor, if necesssary
  if (dt1 .ne. dt0) then
     call cal_vis(dt1,expvis)
  endif

  ! Multiply be integrating factor (always!)
  do i=1,3
     uk(:,:,:,i)=(uk(:,:,:,i) +b10*nlk(:,:,:,i,n0) +b11*nlk(:,:,:,i,n1))*expvis
     nlk(:,:,:,i,n0)=nlk(:,:,:,i,n0)*expvis
  enddo
end subroutine AdamsBashforth


! Set the time step based on the CFL condition and penalization
! stability contidion.
subroutine adjust_dt(dt1,u)
  use vars
  use mpi_header
  implicit none

  real (kind=pr), intent(in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  integer :: mpicode
  real (kind=pr), intent (out) :: dt1
  real(kind=pr),dimension(nf) :: uloc
  real(kind=pr) :: umax
  integer i,j

  if (dt_fixed>0.0) then
     dt1=dt_fixed
  else
     
     ! Find the max velocity for each field.
     do i=1,nf
        j=(i-1)*3
        uloc(i)=maxval(&
             u(:,:,:,1+j)*u(:,:,:,1+j)+&
             u(:,:,:,2+j)*u(:,:,:,2+j)+&
             u(:,:,:,3+j)*u(:,:,:,3+j)&
             )
        uloc(i)=dsqrt(uloc(i))
     enddo

     ! Make uloc(1) the max local velocity for all fields.
     do i=2,nf
        if(uloc(1) < uloc(i)) uloc(1)=uloc(i) 
     enddo

     ! Set umax to be the maximum for all fields over all procs.
     call MPI_REDUCE(uloc(1),umax,1,mpireal,MPI_MAX,0,MPI_COMM_WORLD,mpicode)
     
     !--Adjust time step at 0th process
     if ( mpirank == 0 ) then

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
     endif

     ! Broadcast time step to all processes
     call MPI_BCAST(dt1,1,mpireal,0,MPI_COMM_WORLD,mpicode)
  endif
end subroutine adjust_dt


! FIXME: add documentation
subroutine truncate(a,b)
  ! rounds time step (from 1.246262e-2 to 1.2e-2)
  use fsi_vars
  implicit none

  real(kind=pr) :: a,b
  character (len=7) :: str

  write (str,'(es7.1)') a
  read (str,*) b
end subroutine truncate


! Force zero mode for mean flow
subroutine set_mean_flow(uk,time)
  use mpi_header
  use vars
  implicit none
  
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  real (kind=pr),intent (inout) :: time

  ! Force zero mode for mean flow
  if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
     uk(0,0,0,1)=Ux
     uk(0,0,0,2)=Uy
     uk(0,0,0,3)=Uz
  endif
end subroutine set_mean_flow
