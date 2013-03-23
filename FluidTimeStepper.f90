subroutine FluidTimestep( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis, it )
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
  real (kind=pr) :: t1
  t1 = MPI_wtime()
  
  if (iTimeMethodFluid == 2) then  
    call RungeKutta2 ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )    
  elseif ((iTimeMethodFluid == 1).and.(it==0)) then  
    call Euler_startup ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )
    if (mpirank ==0) then
    write(*,'(A)') "*** info: did startup euler............"
    endif
  elseif (iTimeMethodFluid == 1) then  
    call AdamsBashforth ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )  
  elseif (iTimeMethodFluid == 3) then  
    call Euler ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )      
  endif 
  
  
  !-----------------------------------------------------------------------------------
  ! Force zero mode for mean flow
  !-----------------------------------------------------------------------------------
  ! Force zero mode for mean flow
  if ( (iMeanFlow == 1).and.((ca(1) == 0) .and. (ca(2) == 0) .and. (ca(3) == 0)) ) then
    uk(0, 0, 0, 1) = Ux + Ax * time
    uk(0, 0, 0, 2) = Uy + Ay * time
    uk(0, 0, 0, 3) = Uz + Az * time
  endif
  
  !-----------------------------------
  ! dealiasing
  !-----------------------------------
  if (iDealias==1) then
    call dealiase_velocity (uk)
  endif
  
  time_fluid = time_fluid + MPI_wtime() - t1
end subroutine





subroutine RungeKutta2 ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1 !unused; only for compataility
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: t1,t2,t3,t4
   

  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing (for the euler step)
  !-----------------------------------------------------------------------------------
  call cal_nlk (dt1, nlk(:,:,:,:,0), uk, u, vort, work )
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
  call cal_nlk (dt1, nlk(:,:,:,:,1), uk, u, vort, work ) 
  
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




subroutine Euler ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )
  !-----------------------------------------------------------------------------------
  !-- This is standard Euler-explicit time marching. does not serve as startup scheme for AB2!!
  !-----------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1 !unused; only for compataility
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: t1,t2,t3,t4
  
  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk (dt1, nlk(:,:,:,:,1), uk, u, vort, work )

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





subroutine Euler_startup ( time, dt0, dt1, n0, n1, u, uk, nlk, vort, work, workvis )
  !-----------------------------------------------------------------------------------
  !-- note this is not an optimized Euler. it does things we need only for AB2
  !-----------------------------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: t1,t2,t3,t4
  
  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk (dt1, nlk(:,:,:,:,n0), uk, u, vort, work )

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


end subroutine



subroutine AdamsBashforth ( time, dt0, dt1, n0,n1, u, uk, nlk, vort, work, workvis )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: time, dt1,dt0
  integer, intent (in) :: n0, n1
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (inout):: nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort, u
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: workvis
  real (kind=pr) :: b10,b11

  !-----------------------------------------------------------------------------------
  !--Calculate fourier coeffs of nonlinear rhs and forcing
  !-----------------------------------------------------------------------------------
  call cal_nlk ( dt1, nlk(:,:,:,:,n0), uk, u, vort, work )
    
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



subroutine dealiase_velocity (uk)
  ! Matthieu 21/07/2011 !  Elliptical dealiasing (truncation 2/3 rule)
  ! with STRIDE-1 ordering (X,Y,Z)_physical->(Z,X,Y)_Fourier
  ! Thomas 22/03/2013, modified to use no work array
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  integer :: ix, iy, iz
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: uk
  real(kind=pr) :: kx2, ky2, kz2, kx_trunc, ky_trunc, kz_trunc


  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)

  do iz = ca(1), cb(1)
     kz2=dble(modulo(iz+nz/2,nz)-nz/2)/kz_trunc
     kz2=kz2*kz2
     do iy = ca(3), cb(3)
        ky2=dble(modulo(iy+ny/2,ny)-ny/2)/ky_trunc
        ky2=ky2*ky2
        do ix = ca(2), cb(2)
           kx2=dble(ix)/kx_trunc
           kx2=kx2*kx2
           if (kx2 + ky2 + kz2  .ge. 1.d0) then
              uk(iz, ix, iy, 1:3) = dcmplx(0.d0,0.d0)
           endif
        enddo
     enddo
  enddo

end subroutine dealiase_velocity
