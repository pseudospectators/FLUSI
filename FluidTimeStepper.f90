subroutine FluidTimestep ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, it, GlobIntegrals )
  !-----------------------------------------------------------------------------------
  ! WRAPPER for different time marching methods
  !-----------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1,it
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  type(Integrals), intent (out) :: GlobIntegrals
  real (kind=pr) :: t1
  
  t1 = MPI_wtime()  
  ! note that in the new version, dealiasing is done in cal_vis  
  
  !-----------------------------------------------------------------------------------
  ! call fluid advancement subroutines
  !-----------------------------------------------------------------------------------
  if (iTimeMethodFluid == "RK2") then  
  
    call RungeKutta2    ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )    
    
  elseif ((iTimeMethodFluid == "AB2").and.(it==0)) then
  
    call Euler_startup  ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )    
    
  elseif (iTimeMethodFluid == "AB2") then  
  
    call AdamsBashforth ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )  
    
  elseif (iTimeMethodFluid == "Euler") then  
  
    call Euler          ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )      
    
  else
    if (mpirank == 0) then
      write(*,*) "Error! iTimeMethodFluid unknown. Suicide!"
      stop
    endif
  endif 
    
  !-----------------------------------------------------------------------------------
  ! Force zero mode for mean flow
  !-----------------------------------------------------------------------------------
  if ( (iMeanFlow == 1).and.((ca(1) == 0) .and. (ca(2) == 0) .and. (ca(3) == 0)) ) then
    uk(0, 0, 0, 1) = Ux + Ax * time
    uk(0, 0, 0, 2) = Uy + Ay * time
    uk(0, 0, 0, 3) = Uz + Az * time
  endif
  
 
  time_fluid = time_fluid + MPI_wtime() - t1
end subroutine FluidTimestep





subroutine RungeKutta2 ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1,it !unused; only for compataility
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  type(Integrals), intent (out) :: GlobIntegrals

  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  !-----------------------------------------------------------------------------------
  call cal_nlk (time,it, dt1, nlk(:,:,:,:,0), uk, u, vort, work, GlobIntegrals )
  ! multiply the RHS with the viscosity
  nlk (:,:,:,1,0) = nlk (:,:,:,1,0) * workvis  
  nlk (:,:,:,2,0) = nlk (:,:,:,2,0) * workvis  
  nlk (:,:,:,3,0) = nlk (:,:,:,3,0) * workvis  
  
  
  !-----------------------------------------------------------------------------------
  ! compute integrating factor, only done if necessary (i.e. time step has changed)
  !-----------------------------------------------------------------------------------
  if (dt1 .ne. dt0) then
    call cal_vis (dt1, workvis)
  endif

  !-----------------------------------------------------------------------------------
  !-- do the actual euler step. note nlk is already multiplied by vis
  !----------------------------------------------------------------------------------- 
  uk(:,:,:,1) = (uk(:,:,:,1)*workvis + dt1 * nlk (:,:,:,1,0))
  uk(:,:,:,2) = (uk(:,:,:,2)*workvis + dt1 * nlk (:,:,:,2,0))
  uk(:,:,:,3) = (uk(:,:,:,3)*workvis + dt1 * nlk (:,:,:,3,0))
  
  !-----------------------------------------------------------------------------------
  !-- RHS using the euler velocity
  !----------------------------------------------------------------------------------- 
  call cal_nlk (time,it, dt1, nlk(:,:,:,:,1), uk, u, vort, work ) 
  
  !-----------------------------------------------------------------------------------
  !-- do the actual time step. note the minus sign!!
  !-- in the original formulation, it reads u^n+1 = u^n + dt/2 * ( N(u^n)*vis + N(u_euler) )
  !-- but we don't want to save u_euler seperately, we want to overwrite u^n with it!
  !-- so the formulation reads u^n+1 = u_euler - dt*N(u^n)*vis + dt/2 * ( N(u^n)*vis + N(u_euler) )
  !-- which yields simply
  !-- u^n+1 = u_euler + dt/2 * ( -N(u^n)*vis + N(u_euler) )
  !-----------------------------------------------------------------------------------   
  uk(:,:,:,1) = uk(:,:,:,1) + 0.5*dt1*( -nlk(:,:,:,1,0) + nlk(:,:,:,1,1) )
  uk(:,:,:,2) = uk(:,:,:,2) + 0.5*dt1*( -nlk(:,:,:,2,0) + nlk(:,:,:,2,1) )
  uk(:,:,:,3) = uk(:,:,:,3) + 0.5*dt1*( -nlk(:,:,:,3,0) + nlk(:,:,:,3,1) )   

end subroutine




subroutine Euler ( time, it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )
  !-----------------------------------------------------------------------------------
  !-- This is standard Euler-explicit time marching. does not serve as startup scheme for AB2!!
  !-----------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1,it !unused; only for compataility
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: t1,t2,t3,t4
  type(Integrals), intent (out) :: GlobIntegrals
  
  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk ( time,it, dt1, nlk(:,:,:,:,1), uk, u, vort, work, GlobIntegrals )

  !-----------------------------------------------------------------------------------
  ! compute integrating factor, if necesssary
  !-----------------------------------------------------------------------------------
  if (dt1 .ne. dt0) then
    call cal_vis (dt1, workvis)
  endif
  
  !-----------------------------------------------------------------------------------
  ! Multiply be integrating factor (always!)
  !-----------------------------------------------------------------------------------
  uk(:,:,:,1) = (uk(:,:,:,1) + dt1 * nlk (:,:,:,1,1)) * workvis
  uk(:,:,:,2) = (uk(:,:,:,2) + dt1 * nlk (:,:,:,2,1)) * workvis
  uk(:,:,:,3) = (uk(:,:,:,3) + dt1 * nlk (:,:,:,3,1)) * workvis

end subroutine





subroutine Euler_startup ( time,it, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )
  !-----------------------------------------------------------------------------------
  !-- note this is not an optimized Euler. it does things we need only for AB2
  !-----------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1,it
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: t1,t2,t3,t4
  type(Integrals), intent (out) :: GlobIntegrals
  
  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk ( time,it, dt1, nlk(:,:,:,:,n0), uk, u, vort, work, GlobIntegrals )

  !-----------------------------------------------------------------------------------
  ! compute integrating factor, if necesssary
  !-----------------------------------------------------------------------------------
  if (dt1 .ne. dt0) then
    call cal_vis (dt1, workvis)
  endif
  
  !-----------------------------------------------------------------------------------
  ! Multiply be integrating factor (always!)
  !-----------------------------------------------------------------------------------
  uk(:,:,:,1) = (uk(:,:,:,1) + dt1 * nlk (:,:,:,1,n0)) * workvis
  uk(:,:,:,2) = (uk(:,:,:,2) + dt1 * nlk (:,:,:,2,n0)) * workvis
  uk(:,:,:,3) = (uk(:,:,:,3) + dt1 * nlk (:,:,:,3,n0)) * workvis
  nlk (:,:,:,1,n0) = nlk (:,:,:,1,n0) * workvis
  nlk (:,:,:,2,n0) = nlk (:,:,:,2,n0) * workvis
  nlk (:,:,:,3,n0) = nlk (:,:,:,3,n0) * workvis
  
  
  if (mpirank ==0) write(*,'(A)') "*** info: did startup euler............"
end subroutine



subroutine AdamsBashforth ( time, it, dt0, dt1, n0,n1, u, uk, nlk, vort, work, workvis, GlobIntegrals )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1,it
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: b10,b11
  type(Integrals), intent (out) :: GlobIntegrals

  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk ( time,it, dt1, nlk(:,:,:,:,n0), uk, u, vort, work, GlobIntegrals )
    
  !-----------------------------------------------------------------------------------
  !--Calculate velocity at new time step 
  !--(2nd order Adams-Bashforth with exact integration of diffusion term)
  !-----------------------------------------------------------------------------------
  b10 = dt1/dt0 * (0.5*dt1 + dt0)
  b11 = -0.5 * dt1**2 / dt0

  !-----------------------------------------------------------------------------------
  ! compute integrating factor, if necesssary
  !-----------------------------------------------------------------------------------
  if (dt1 .ne. dt0) then
    call cal_vis(dt1, workvis)
  endif
    
  !-----------------------------------------------------------------------------------
  ! Multiply be integrating factor (always!)
  !-----------------------------------------------------------------------------------
  uk(:,:,:,1) = (uk(:,:,:,1) + b10*nlk(:,:,:,1,n0) + b11*nlk(:,:,:,1,n1)) * workvis
  uk(:,:,:,2) = (uk(:,:,:,2) + b10*nlk(:,:,:,2,n0) + b11*nlk(:,:,:,2,n1)) * workvis
  uk(:,:,:,3) = (uk(:,:,:,3) + b10*nlk(:,:,:,3,n0) + b11*nlk(:,:,:,3,n1)) * workvis
  nlk (:,:,:,1,n0) = nlk (:,:,:,1,n0) * workvis
  nlk (:,:,:,2,n0) = nlk (:,:,:,2,n0) * workvis
  nlk (:,:,:,3,n0) = nlk (:,:,:,3,n0) * workvis

end subroutine