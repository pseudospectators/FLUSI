module SolidSolver
use share_vars
use motion
implicit none

! here we could define some variables

 contains

!-------------------------------------------------------------------------------------
!   SOLID SOLVER WRAPPER
!-------------------------------------------------------------------------------------
subroutine SolidSolverWrapper ( time, dt, beams )
  use share_vars
  implicit none
  real (kind=pr), intent (in) ::                     dt, time
  type(solid), dimension(1:iBeam), intent (inout) ::    beams
  integer :: i
  
  if (time>T_release) then ! it is not nessesaire to solve the solid equation when the beam is still held fixed  
     !-------------------------------------------
     ! the beams are released, call IBES solvers
     !-------------------------------------------
     ! all implicit solvers are in one subroutine
     do i = 1, iBeam
     call IBES_solver (time, dt, beams(i))     
     enddo
  else 
     !-------------------------------------------
     ! the beams are not yet released, but its leading edges may move
     !-------------------------------------------
     do i = 1, iBeam
     call integrate_position (time+dt, beams(i))
     enddo
  endif
  
  !-------------------------------------------
  ! compute energies and stuff
  !-------------------------------------------
  do i = 1, iBeam
  call SolidEnergies( beams(i) )
  enddo

end subroutine SolidSolverWrapper






!-------------------------------------------------------------------------------------
!   SOLID SOLVER INITIALIZATION
!-------------------------------------------------------------------------------------
subroutine InitializeSolidSolver( beams )
  use share_vars
  implicit none
  integer :: i
  type(solid), dimension(1:iBeam), intent (inout) :: beams
  
  ! marks all beams to be in the very first time step  
  ! the solver then uses CN2 instead of BDF2, because the old old time level
  ! t_n-1 is not available
  do i=1,iBeam
    beams(i)%StartupStep = .true.
  enddo
  
end subroutine InitializeSolidSolver





!-------------------------------------------------------------------------------------
!   energies and stuff for beams
!-------------------------------------------------------------------------------------
subroutine SolidEnergies( beam )
  use share_vars
  implicit none
  type(solid), intent (inout) :: beam
  real (kind=pr), dimension (0:ns-1) :: theta_s
  
  ! for beam elastic energy, we need theta_s over the beam 
  call Differentiate1D ( beam%theta, theta_s, ns, ds, 1)  
  
  beam%E_kinetic = mue*0.5d0*ds*sum( beam%vx**2 + beam%vy**2 )
  beam%E_pot     = mue*grav*ds *sum( beam%y-beam%y0 )
  beam%E_elastic = eta*0.5d0*ds*sum( theta_s**2 )
  
  beam%Inertial_Force(1) = mue*ds* sum( beam%ax ) 
  beam%Inertial_Force(2) = mue*ds* sum( beam%ay ) 
  
end subroutine SolidEnergies



!-------------------------------------------------------------------------------------
!   SOLID SOLVER ROUTINES
!-------------------------------------------------------------------------------------
subroutine IBES_solver ( time, dt, beam_solid )! note this is actuall only ONE beam!!
  use share_vars
  use mkl_pardiso
  use MKL_PARDISO_PRIVATE
  use mkl95_lapack  
  use mkl95_precision  
  use CompressMatrixCSR
  use omp_lib
  use PerformanceMeasurement
  implicit none
  real (kind=pr), intent (in) :: 				dt, time
  type (solid), intent(inout) ::                                beam_solid
  real (kind=pr), dimension (0:ns-1) ::				T
  real (kind=pr), dimension (0:ns-1) ::				old_rhs, theta_old,T_dummy
  real (kind=pr), dimension (0:ns-1) ::	    	                pressure_old, pressure_new, tau_beam_old, tau_beam_new
  real (kind=pr), dimension (0:ns-1, 1:6) :: 	                beam
  real (kind=pr), dimension (-1:ns+1) :: 			theta_guess
  real (kind=pr), dimension (-1:ns-1) :: 			T_guess
  real (kind=pr), dimension (0:ns-1, 1:6) :: 			beam_old, beam_guess
  real (kind=pr), dimension (0:ns-1, 1:6) :: 	                beam_oldold
  real (kind=pr), dimension (1:2*ns+4) :: 			x_guess, x, x_delta, F
  real (kind=pr), dimension (1:2*ns+4,1:2*ns+4)  :: 		J,J2, J2_norm
  real (kind=pr), dimension (1:ns+3) :: 			theta_act
  real (kind=pr), dimension (1:ns+1) :: 			T_act
  real (kind=pr) :: 						alpha, alpha_t, alpha_tt, err, A1, A2, K1, K2, T_s0, theta_ss0, theta_s0, C2,C1,C3,C4
  real (kind=pr) :: 						err_rel, R
  real (kind=pr) ::					        dt_old
  real (kind=pr), parameter ::					error_stop = 1.0e-7
  integer :: 							n,N_nonzero,k,iter
  integer, save ::						iCalls=10
  logical :: 							ActuallyBDF2=.false., iterate=.true.
  real (kind=pr), dimension(1:6) :: 				LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  

  !******************************************************* 
  ! NOTE: 2013, this is extended to take more than one beam into account.
  ! however, its too hard to reprogram everything, so now
  ! we force it to be compatible
  !*******************************************************
  beam(:,1) = beam_solid%x
  beam(:,2) = beam_solid%y
  beam(:,3) = beam_solid%vx
  beam(:,4) = beam_solid%vy
  beam(:,5) = beam_solid%theta
  beam(:,6) = beam_solid%theta_dot
  pressure_old = beam_solid%pressure_old
  pressure_new = beam_solid%pressure_new
  tau_beam_old = beam_solid%tau_old
  tau_beam_new = beam_solid%tau_new
  dt_old       = beam_solid%dt_old
  beam_oldold  = beam_solid%beam_oldold
  !*******************************************************
  
  call mouvement ( time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
    
  beam_old = beam
  beam_guess = beam  

  !--------------------------------------------------------------------
  ! 	Startup time step: use CN2 as first step of BDF2 (pay attention to call InitializeSolidSolver() before the run)
  !-------------------------------------------------------------------- 
  if (beam_solid%StartupStep) then	! this is the first time step 
    beam_solid%StartupStep = .false.	! we're about to do the first step
    if (TimeMethodSolid==BDF2) then	! if we deal with BDF2
      ActuallyBDF2 = .true.		! Remember to switch back to BDF2 (at the end of the step)
      TimeMethodSolid = CrankNicholson	! use CrankNicholson for the first step
    endif
  endif
  

  !--------------------------------------------------------------------
  ! 	Initial guess for the beam at the new time level. 
  ! 	An euler-explicit step is made, ghostpoints are added 
  !--------------------------------------------------------------------
  
  call EE1 (time, dt, beam_guess, pressure_old, T, tau_beam_old, beam_solid)
  ! now we compute the ghostpoints Theta(-1); Theta(ns); Theta(ns+1); T(-1)
  theta_guess(0:ns-1) = beam_guess(0:ns-1,5) !because of this, temp(1) = 0 (second point added for boundarys)
  !-- leading edge boundary conditions (constants)
  K2 = pressure_new(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))
  C2 = 2.0*ds*K2 - T(0)*theta_guess(1) + (eta/ds**2)*(10.0*theta_guess(0)-12.0*theta_guess(1)+6.0*theta_guess(2)-theta_guess(3) )
  !-- theta_extended(0) is the first virtual node. 
  theta_guess(-1) = C2 / ( (3.0*eta/ds**2)-T(0) )
  !-- solve last two boundary conditions: theta_s(ns-1) = theta_ss(ns-1) = 0
  A1 = beam_guess(ns-3,5)/12. - 2.*beam_guess(ns-2,5)/3.
  A2 =-beam_guess(ns-3,5)/12. + 4.*beam_guess(ns-2,5)/3. - 5.*beam_guess(ns-1,5)/2.
  !-- from these equations we can calculate the values of theta(ns) and theta(ns+1) which are added points
  theta_guess(ns) =  3. * (A1-A2)/2.
  theta_guess(ns+1) = 12. * (2.*A1-A2)
  

  ! ------ set virtual node for Tension T
  ! constant of the BC, see masters thesis
  K1 = mue*(LeadingEdge(5)*cos(alpha)+LeadingEdge(6)*sin(alpha)+grav*sin(alpha))  - tau_beam_new(0) !technically, the guess is at the new time level
  ! value of the inhomogenous Neumann BC at the leading edge
  theta_ss0 = ( 2.0*theta_guess(0) - 5.0*theta_guess(1) +4.0*theta_guess(2) - 1.0*theta_guess(3) )/ds**2
  theta_s0  = (-1.5*theta_guess(0) + 2.0*theta_guess(1) -0.5*theta_guess(2)  )/ds  
  T_s0 = -eta * theta_ss0*theta_s0 + K1 
  T_guess(-1) = T(1) - 2.0*ds*T_s0 !T_s = 0 means T(-1) = T(1)  
  T_guess(0:ns-1) = T(0:ns-1)  
  
  
  ! ------- rearrange the guess in the vector of unknowns x
  x_guess(1:ns+3)      = theta_guess(-1:ns+1)
  x_guess(ns+4:2*ns+4) = T_guess(-1:ns-1)  
  !now we have a complete initial guess for the variables

  !--------------------------------------------------------------------
  ! 	Calculate RHS vector @ t_n. This is required for the CN2 method only.
  !--------------------------------------------------------------------
  if (TimeMethodSolid == CrankNicholson) then
    theta_old = beam_old(:,5)
    old_rhs = beam_old(:,6)           
    call RHS_beameqn(time, theta_old, old_rhs, pressure_old, T_dummy, tau_beam_old, beam_solid)
  endif  
  !--------------------------------------------------------------------
  ! 	Newton iteration
  !--------------------------------------------------------------------
  err     = 1.0
  err_rel = 1.0
  
  x = x_guess
    
  iter = 0
  
  ! MODIFICATION 31.10.2012: the iteration now uses both relative and absolute error. If x is large, we can have trouble
  ! reaching a very small increment, and iterate forever. If x is small, then the relative criterion tries to go much below
  ! machine precision. So we use both, either with the same precision.
  iterate = .true.
  do while (iterate==.true.)
    theta_act = x(1:ns+3)
    T_act     = x(ns+4:2*ns+4)
    !-------------------------------------------------------------------------
    !	Calculate RHS vector
    !-------------------------------------------------------------------------
    call F_nonlinear( time, dt, dt_old, F, beam_old(:,5), beam_old(:,6), theta_act, T_act, pressure_new, old_rhs, beam_oldold(:,5), beam_oldold(:,6), tau_beam_new, beam_solid)
    F = -1.0*F !newton raphson is J*dx = -F

    !-------------------------------------------------------------------------
    !	Create Jacobi Matrix and Compress it to the sparse format.
    !-------------------------------------------------------------------------    
    call Jacobi(time, dt, dt_old, J, N_nonzero , T_act, theta_act, beam_old(:,5), beam_old(:,6), pressure_new, beam_oldold(:,5), beam_oldold(:,6), tau_beam_new, beam_solid)

    !-------------------------------------------------------------------------
    !	self-test. occasionally, check if Jacobian is okay
    !-------------------------------------------------------------------------     
    if (mod(iCalls,4607)==0) then
	call Jacobi_num(time, dt,dt_old, J2, k, T_act, theta_act, beam_old(:,5), beam_old(:,6), pressure_new, beam_oldold(:,5), beam_oldold(:,6), old_rhs, tau_beam_new, beam_solid)
	J2_norm = J2
	where (abs(J2_norm)<1.0e-7) J2_norm=1.0
	J2_norm = (J-J2)/J2_norm
	where (abs(J2_norm)<1.0e-5) J2_norm=0.0 ! delete small values
	if (maxval(J2_norm)>1.0e-2) then	
	    open (14, file = 'IBES_JACOBIAN', status = 'unknown', access = 'append') ! Append output data file
	    write (14,*) "---analytic---"
	    do n=1,ns*2+4
	    write (14,'(1x,3096(es8.1,1x))') J(:,n)
	    enddo
	    write (14,*) "---numeric---"
	    do n=1,ns*2+4
	    write (14,'(1x,3096(es8.1,1x))') J2(:,n)
	    enddo
	    write (14,*) "---difference---"
	    do n=1,ns*2+4
	    write (14,'(1x,3096(es8.1,1x))') J2_norm(:,n)
	    enddo   
	    call DisplayError
	    write(*,'(A)') "!!! A weird error in IBES occured. Now and then, I compute the exact Jacobian with finite differences."
	    write(*,'(A)') "    This time, the analytical Jacobian and the numerical one differ more than 1% from each other"
	    write(*,'(A)') "    This means something went wrong in IBES, and you should check it."
	    write(*,'(A)') "    Write an email to thomas.engels@mail.tu-berlin.de"
	    write (*,*) time, iter
	    close (14)
	    stop
	endif
	iCalls = 0
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SOLVE LINEAR SYSTEM
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call Solve_LGS ( J, F, x_delta, N_nonzero, "DIRECT")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    iter    = iter + 1
    x       = x + x_delta    
    err     = sqrt(sum(x_delta**2))
    err_rel = abs(sqrt(sum(x_delta**2)) / sqrt(sum(x**2)))

    !-----------------------------------------------------------
    ! CONVERGENCE CRITERION
    !-----------------------------------------------------------
    if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then ! absolute error criterion
      iterate = .false.
    endif
    
    !-----------------------------------------------------------
    ! EMERGENCY BRAKE
    !-----------------------------------------------------------
    if (iter>49) then
      write(*,*) "!!! ERROR: IBES performed like 500 iterations. this is not normal."
      stop
    endif
    
  enddo

  
  !--------------------------------------------------------------------
  !		cut ghostpoints 
  !--------------------------------------------------------------------
  beam(:,5)=x(2:ns+1)
  T = x(ns+5:2*ns+4)  !nessesairy if one does not use the EE1 predictor to guess
                      !attention! T(-1) is not part of T (Ghostpoint)
  !---------------------------------------------------------------------
  ! 		compute angular velocity (theta_dot)
  !---------------------------------------------------------------------  
  ! theta_dot was removed from the NL system and is now computed  
  if (TimeMethodSolid == EulerImplicit) then  
    C1=1.0 ! dt factor
    C2=0.0 ! factor for RHS
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == CrankNicholson) then  
    C1=2.0 ! dt factor
    C2=1.0 ! rhs old factor
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == BDF2) then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R) 	! dt factor
    C2 = 0.0			! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R) 	! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R) 	! factor before the THETA_DOT_N-1 term 
  endif
  
  beam(:,6) = (C1/dt) * ( beam(:,5) - C3*beam_old(:,5) - C4*beam_oldold(:,5) ) - C2*beam_old(:,6)
  
  !*******************************************************
  ! NOTE: 2013, this is extended to take more than one beam into account.
  ! however, its too hard to reprogram everything, so now
  ! we force it to be compatible
  !*******************************************************
  beam_solid%x = beam(:,1)
  beam_solid%y = beam(:,2)
  beam_solid%vx = beam(:,3)
  beam_solid%vy = beam(:,4)
  beam_solid%theta = beam(:,5)
  beam_solid%theta_dot = beam(:,6)
  beam_solid%pressure_old = pressure_old 
  beam_solid%pressure_new = pressure_new 
  beam_solid%tau_old = tau_beam_old
  beam_solid%tau_new = tau_beam_new
  !*******************************************************
  
  !---------------------------------------------------------------------
  ! 		get deflection line (integrate_position)
  !---------------------------------------------------------------------  
  ! this beam is the new one, so its time+dt ( to get mouvement at the right instant)
  call integrate_position ( time+dt, beam_solid )
  
  !---------------------------------------------------------------------
  !  accelerations
  !---------------------------------------------------------------------  
  ! we have the old velocity (t_n) and the new one 
  beam_solid%ax = (beam_solid%vx - beam_old(:,3)) / dt
  beam_solid%ay = (beam_solid%vy - beam_old(:,4)) / dt
  
  
  !---------------------------------------------------------------------
  !     	emergency brake
  !---------------------------------------------------------------------
  if (maxval(abs(beam(:,6)-beam_old(:,6) ))>100.0 ) then
    write (*,'(A)') "??????????????????????????????????????????????????????????????????????????????????????????????????????"
    call DisplayError()
    write (*,'(A)') " "
    write (*,'(A)') "!!! IBES-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.0"
    write (*,'(A)') "    That indicates a problem, probably artificial added mass instability."
    write (*,'(A)') "    Even though it won't make sense to continue, I make a backup, quit, and let you decide"   
    write (*,'("time=",es11.4)') time
    write (*,'(A)') "??????????????????????????????????????????????????????????????????????????????????????????????????????" 
    continue_timestep = .false.
  endif
  if (iter>1000) then
    call DisplayError()
    write (*,'(A)') "!!! IBES solver: It took more than 1000 iterations for a single time step. This indicates convergence problems"
    write (*,'(A)') "in the Jacobi-iteration. Reduce the time step and try again"
    continue_timestep = .false.
  endif
    
  
  !---------------------------------------------------------------------
  ! 		save number of iterations
  !---------------------------------------------------------------------    
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'IBES_iter', status = 'unknown', access = 'append') ! Append output data file
  write (14, '(es11.4,1x,i3)') time, iter
  close (14)
  
  !---------------------------------------------------------------------
  ! 		iterate ( skipped for CN2 and EI1 )
  !---------------------------------------------------------------------    
  
  if (ActuallyBDF2 == .true.) then 	! Remember to switch back to BDF2 
    TimeMethodSolid = BDF2 
    ActuallyBDF2 = .false.
  endif  

  iCalls = iCalls + 1			! count calls (to perform rare self-tests)

  ! The old beam at the current step is the oldold (t_n-1) at the next step
  beam_solid%beam_oldold = beam_old 
  ! for BDF2 with variable dt, we need the old time step
  ! this should be okay also when restarting, as the first CN2 step will provide us with the dt_old  
  beam_solid%dt_old = dt 	
  

end subroutine IBES_solver






subroutine Solve_LGS ( J, F, x, N_nonzero, solver)
  !--------------------------------------------
  ! solves the linear system J*x = F
  !--------------------------------------------
  use share_vars
  use CompressMatrixCSR
  use mkl_pardiso
  use MKL_PARDISO_PRIVATE
  use mkl95_lapack  
  use mkl95_precision
  use omp_lib
  implicit none
  real (kind=pr), dimension (1:2*ns+4,1:2*ns+4), intent(in) :: J
  real (kind=pr), dimension (1:2*ns+4), intent(out) :: x
  real (kind=pr), dimension (1:2*ns+4), intent(out) :: F
  character(len=*), intent(in) :: solver
  integer, intent(in) :: N_nonzero
  
  real (kind=pr), allocatable, dimension(:) :: values
  integer, allocatable, dimension(:) :: columns, rows
  integer, dimension(2*ns+4) :: perm
  integer, dimension(64) :: iparams
  integer :: error, ipiv(1:2*ns+4), i
  type (MKL_PARDISO_HANDLE) :: PT(64) 
  
  real (kind=pr), dimension (1:2*ns+4,1:2*ns+4) :: J2
  real (kind=pr), dimension (1:2*ns+4) :: F2
  
  x = 0.0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (solver == "PARDISO") then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    !-----------------------------------------------------------------
    !--Compress the Jacobian into CSR (Compressed sparse row) format
    !-----------------------------------------------------------------
    allocate( values(1:N_nonzero), columns(1:N_nonzero) )
    allocate( rows(1:2*ns+5) ) !rows for the CSR format has a fixed dimension
    call CompressMatrix(J, values, columns, rows)  
    
    !-----------------------------------------------------------------
    ! initialize PARDISO
    !-----------------------------------------------------------------    
    !pardiso internal adresses
    pt%dummy = 0 
    perm = 0  
    iparams = 0   !use standard parameters for pardiso
    iparams(1)=1  !don't use std parameters, specified below
    iparams(2)=3  !0 minum degree 2 nestes METIS 3 OPENMP
    iparams(3)=omp_get_num_threads()  !mkl_get_max_threads() ! number of threads
    iparams(4)=0
    iparams(5)=0
    iparams(6)=0
    iparams(8)=0
    iparams(10)=13
    iparams(11)=1
    iparams(13)=1
    iparams(18)=0
    iparams(19)=0
    iparams(21)=1
    iparams(27)=1
    iparams(28)=0
    iparams(35)=0
    iparams(60)=0 ! in core o/ outofcore(2)
    
    
    !-----------------------------------------------------------------
    ! SOLVE using PARDISO
    !-----------------------------------------------------------------  
    call pardiso( pt, 1, 1, 11, 13, 2*ns+4, values, rows, columns, perm, 1, iparams, 0, F, x, error)
!   call pardiso(pt, maxfct, mnum, type, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)

    if (error .ne. 0) then
      write(*,*) "!!! Crutial: PARDISO error.", error
      stop 
    endif
    
    
    !-----------------------------------------------------------------
    ! CLEANUP
    !-----------------------------------------------------------------      
    !release pardiso internal memory
    call pardiso( pt, 1, 1, 11, -1, 2*ns+4, values, rows, columns, perm, 1, iparams, 0, F, x, error)    
    deallocate ( rows,values,columns )
    
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  elseif (solver == "DIRECT") then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    J2 = transpose(J)
    F2 = F
    call dgetrf( 2*ns+4, 2*ns+4, J2, 2*ns+4, ipiv, error )
    if (error .ne. 0) then
      write(*,*) "!!! Crutial: dgetrf error.", error
      stop
    endif
    
    
    call dgetrs( 'N', 2*ns+4, 1, J2, 2*ns+4, ipiv, F2, 2*ns+4, error )
    if (error .ne. 0) then 
      write(*,*) "!!! Crutial: dgetrs error.", error
      stop
    endif
    
    x = F2
    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  else
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
      write (*,*) "!!! linear solver: method unkown"
      stop 
  endif
  
  !!!!!!!!!!!!!!!!!!!!
  ! check for NaN's
  !!!!!!!!!!!!!!!!!!!!
  do i = 1, 2*ns + 4
  if (isnan(x(i))) then
    write (*,*) "Solver: "//solver
    write (*,*) "OhOh. Something went wrong in the Linear solver for the solid"
    write (*,*) "May I suggest to kill yourself as you'll never find this mistake?"
    stop
  endif
  enddo
  
end subroutine

!##################################################################################################################################################################

subroutine F_nonlinear (time, dt, dt_old,  F, theta_old, theta_dot_old, theta, T, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam, beam_solid  )
  use share_vars
  implicit none  
  ! returns the RHS of the nonlinear eqn set F(x) for a given beam (the iterating one) and the (fixed) previous one
  real (kind=pr), dimension (0:ns-1), intent (in) :: theta_old, theta_dot_old , old_rhs, theta_oldold, theta_dot_oldold
  real (kind=pr), dimension (-1:ns+1), intent (in) :: theta
  type (solid) :: beam_solid
  real (kind=pr), dimension (-1:ns-1), intent (in) :: T
  real (kind=pr), dimension (0:ns-1), intent (in) :: p, tau_beam
  real (kind=pr), dimension (0:ns-1) :: theta_dot_new, tau_s
  real (kind=pr), dimension (1:(2*ns+4)), intent (out) :: F  
  real (kind=pr), intent (in) :: time, dt, dt_old
  real (kind=pr) :: K1,K2, C1,C2, C3,C4,R
  real (kind=pr) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  integer :: i
  
  call Differentiate1D (tau_beam, tau_s, ns, ds, 1)
  
  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  !---------new extended boundarys (constants)
  K1 = mue*(LeadingEdge(5)*cos(alpha)+LeadingEdge(6)*sin(alpha)+grav*sin(alpha)) - tau_beam(0)
  K2 = p(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))

  F = 0. !to check if an index is not set. is not the case.
  
  if (TimeMethodSolid == EulerImplicit) then  
    C1=1.0 ! dt factor
    C2=0.0 ! factor for RHS
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == CrankNicholson) then  
    C1=2.0 ! dt factor
    C2=1.0 ! rhs old factor
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == BDF2) then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R) 	! dt factor
    C2 = 0.0			! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R) 	! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R) 	! factor before the THETA_DOT_N-1 term 
  endif
  
  theta_dot_new = (C1/dt) * ( theta(0:ns-1) - C3*theta_old(0:ns-1) - C4*theta_oldold ) - C2*theta_dot_old
  
  call Check_Vector_NAN ( theta_old, "theta_old" )
  call Check_Vector_NAN ( theta_dot_old, "theta_dot_old" )
  call Check_Vector_NAN ( theta, "theta" )
  call Check_Vector_NAN ( T, "T" )
  call Check_Vector_NAN ( p, "p" )
  call Check_Vector_NAN ( old_rhs, "old_rhs" )
  call Check_Vector_NAN ( theta_oldold, "theta_oldold" )
  call Check_Vector_NAN ( theta_dot_oldold, "theta_dot_oldold" )
  call Check_Vector_NAN ( tau_beam, "tau_beam" )
  call Check_Vector_NAN ( theta_dot_new, "theta_dot_new" )
  call Check_Vector_NAN ( tau_s, "tau_s" )
  call Check_Vector_NAN ( theta_oldold, "theta_oldold" )
  
  
  
  
  ! first we set the 8 special eqns at the beginning 
  F(1) = theta(0) 
  F(2) = (T(1)-T(-1))/(2.0*ds)  +  (eta/(2.0*ds**3))*(theta(-1)-2.0*theta(0)+theta(1))*(theta(1)-theta(-1))  -  K1
  F(3) = T(0)*(theta(1)-theta(-1))/(2.0*ds)   -   (eta/(2.0*ds**3))*(-3.0*theta(-1)+10.0*theta(0)-12.0*theta(1)+6.0*theta(2)-theta(3))  -  K2
  F(4) = ( (1.0/12.0)*theta(ns-3) - (2.0/3.0)*theta(ns-2) + (2.0/3.0)*theta(ns) - (1.0/12.0)*theta(ns+1) )/ds
  F(5) = ( (-1.0/12.0)*theta(ns-3) + (4.0/3.0)*theta(ns-2) - 2.5*theta(ns-1) + (4.0/3.0)*theta(ns) + (-1.0/12.0)*theta(ns+1) )/ds**2
  F(6) = T(ns-1) 
  
  F(7) = ( T(-1)-2.0*T(0)+T(1) )/(ds**2) &
    - T(0)*((theta(1)-theta(-1))/(2.0*ds) )**2 & 
    + p(0)*((theta(1)-theta(-1))/(2.0*ds)) &
    + (eta/(ds**4))*(theta(1)-theta(-1))*(-2.5*theta(0)+9.0*theta(1)-12.0*theta(2)+7.0*theta(3)-1.5*theta(4)) & 
    + eta*((theta(1)-2.0*theta(0)+theta(-1))/(ds**2))**2 &
    + mue*(theta_dot_new(0) + alpha_t)**2  &
    - tau_s(0)
  
  F(8) = C3*theta_dot_old(ns-1) + C4*theta_dot_oldold(ns-1) - theta_dot_new(ns-1) &
    + (dt/(C1*mue)) * ( -(3.*p(ns-1)-4.*p(ns-2)+p(ns-3))/(2.*ds) & 
    - eta*(theta(ns+1)-4.*theta(ns)+6.*theta(ns-1)-4.*theta(ns-2)+theta(ns-3))/(ds**4) & 
    + (T(ns-1) + eta*( (theta(ns)-theta(ns-2))/(2.*ds) )**2)*(theta(ns)-2.*theta(ns-1)+theta(ns-2))/(ds**2) & 
    + (3.*T(ns-1)-4.*T(ns-2)+T(ns-3))*(theta(ns)-theta(ns-2))/(2.*ds**2) &
    - alpha_tt*(mue**2) &
    - sigma*theta_dot_new(ns-1) & 
    + tau_beam(ns-1)*(theta(ns)-theta(ns-2))/(2.0*ds) &
    + C2*mue*old_rhs(ns-1) ) 

  do i=1, ns-2 
    ! first block: eqn for theta (ns+2 eqn's) [F]
    F( 8+i ) = C3*theta_dot_old(i) + C4*theta_dot_oldold(i) - theta_dot_new(i) &
    + (dt/(C1*mue))*( -1.0*(p(i+1)-p(i-1))/(2.*ds) & 
    - eta*(theta(i+2)-4.*theta(i+1)+6.*theta(i)-4.*theta(i-1)+theta(i-2))/(ds**4) & 
    + (T(i) + eta*( (theta(i+1)-theta(i-1))/(2.*ds) )**2)*(theta(i+1)-2.*theta(i)+theta(i-1))/(ds**2) & 
    + (T(i+1) - T(i-1))*(theta(i+1)-theta(i-1))/(2.0*ds**2) &
    - alpha_tt*(mue**2) & 
    - sigma*theta_dot_new(i) &
    + tau_beam(i)*(theta(i+1)-theta(i-1))/(2.0*ds)&
    + C2*mue*old_rhs(i))  !mue*RHS because in RHS_beameqn term is divided by MUE!!!!
   
    ! 2nd block: eqn's for T (ns+2 eqn's) [G]
    F( 8+(ns-2)+i ) = (T(i-1)-2.0*T(i)+T(i+1))/(ds**2) - T(i)*( (theta(i+1)-theta(i-1))/(2.0*ds) )**2 & !alles ok
    + p(i)*(theta(i+1)-theta(i-1))/(2.0*ds) &
    + eta*(theta(i+1)-theta(i-1))*(theta(i+2)-2.0*theta(i+1)+2.0*theta(i-1)-theta(i-2))/(2.0*ds**4) &
    + eta*( (theta(i+1)-2.0*theta(i)+theta(i-1))/(ds**2) )**2 &
    + mue*( theta_dot_new(i) + alpha_t )**2 &
    - tau_s(i)
    
  enddo 

  !------------------------------------------
  ! check for NaN's
  !------------------------------------------
  call Check_Vector_NAN (F, "F_nonlinear")

end subroutine F_nonlinear

!##################################################################################################################################################################

subroutine Jacobi(time, dt, dt_old, J, N_nonzero , T, theta, theta_old, theta_dot_old, p, theta_oldold, theta_dot_oldold, tau_beam, beam_solid)
  use share_vars
  implicit none
  real (kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real (kind=pr), dimension(-1:ns-1), intent (in) :: T
  real (kind=pr), intent (in) :: time, dt, dt_old
  type(solid) :: beam_solid
  real (kind=pr), dimension(-1:ns+1), intent (in) :: theta
  real (kind=pr), dimension(0:ns-1), intent (in) :: p, theta_old, theta_dot_old, theta_oldold, theta_dot_oldold, tau_beam
  real (kind=pr), dimension(0:ns-1) ::theta_dot_new
  integer, intent (out) :: N_nonzero
  integer :: i, T0_index, k, l, m
  real (kind=pr) :: alpha, alpha_t, alpha_tt, C1,C2,C3,C4,D,R
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  
  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  !---------------------------------------
  ! ---indexing:
  ! theta natural index -1...ns+1
  ! dF/dtheta_i begins @ 1
  !    dTheta_-1   = 1
  !    dTheta_0    = 2
  !    dTheta_ns+1 = ns+4
  !    dTheta_ns   = ns+3
  !    dT_-1       = ns+5
  !    dT_0        = ns+6
  !    dT_1        = ns+7
  !    dT_ns-1     =2ns+4
  !    dT_ns-2     =2ns+3
  !---------------------------------------
  
  if (TimeMethodSolid == EulerImplicit) then  
    C1=1.0 								! dt factor
    C2=0.0 								! factor for RHS
    C3=1.0 								! factor before the THETA_DOT_N term
    C4=0.0 								! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == CrankNicholson) then  
    C1=2.0 								! dt factor
    C2=1.0 								! rhs old factor
    C3=1.0 								! factor before the THETA_DOT_N term
    C4=0.0 								! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == BDF2) then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R) 	! dt factor
    C2 = 0.0			! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R) 	! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R) 	! factor before the THETA_DOT_N-1 term 
  endif
  D  = C1/dt ! this is the only entry in the jacobian for theta_dot_new
  theta_dot_new = (C1/dt) * ( theta(0:ns-1) - C3*theta_old(0:ns-1) - C4*theta_oldold ) - C2*theta_dot_old
  
  
  call Check_Vector_NAN ( theta_old, "theta_old" )
  call Check_Vector_NAN ( theta_dot_old, "theta_dot_old" )
  call Check_Vector_NAN ( theta, "theta" )
  call Check_Vector_NAN ( T, "T" )
  call Check_Vector_NAN ( p, "p" )
  call Check_Vector_NAN ( theta_oldold, "theta_oldold" )
  call Check_Vector_NAN ( theta_dot_oldold, "theta_dot_oldold" )
  call Check_Vector_NAN ( tau_beam, "tau_beam" )
  call Check_Vector_NAN ( theta_dot_new, "theta_dot_new" )
  call Check_Vector_NAN ( theta_oldold, "theta_oldold" )
  
  
 
  T0_index = ns+4 							!T0 means here the first T = T(-1)
  J = 0.0								! initialize J
  !-- set first 6 eqns 
  J(2,1) = 1.0
  !----eqn 2 (special BC 1)
  J(1,2) = (eta/(2.0*ds**3))*(-2.0*theta(-1)+2.0*theta(0))		! dF2 / dTheta(-1)
  J(2,2) =-(eta/(ds**3))*(theta(1)-theta(-1))				! dF2 / dTheta(0)
  J(3,2) = (eta/(2.0*ds**3))*(-2.0*theta(0) +2.0*theta(1))		! dF2 / dTheta(1)
  J(T0_index,2)   = -0.5/ds						! dF2 / dT(-1)
  J(T0_index+2,2) =  0.5/ds						! dF2 / dT(+1)
  !----eqn 3 (special BC 2)
  J(1,3) = (-0.5*eta/(ds**3))*(-3.0) - T(0)/(2.0*ds)			! dF3 / dTheta(-1)
  J(2,3) = (-0.5*eta/(ds**3))*(10.0)					! dF3 / dTheta(0)
  J(3,3) = (-0.5*eta/(ds**3))*(-12.0)+ T(0)/(2.0*ds)			! dF3 / dTheta(+1)
  J(4,3) = (-0.5*eta/(ds**3))*(6.0)					! dF3 / dTheta(+2)
  J(5,3) = (-0.5*eta/(ds**3))*(-1.0)					! dF3 / dTheta(+3)
  J(T0_index+1,3) = (theta(1)-theta(-1))/(2.0*ds)			! dF3 / dT(0)
  
  
  ! the following to eqns have been changed; actually, you could multiply them by any factor, but its easier to read like this
  J(ns+3,4) = (-1.0/12.0)/ds
  J(ns+2,4) = ( 2.0/3.0) /ds
  J(ns,4)   = (-2.0/3.0) /ds
  J(ns-1,4) = ( 1.0/12.0)/ds

  J(ns+3,5) = (-1.0/12.0)/ds**2
  J(ns+2,5) = ( 4.0/3.0 )/ds**2 !
  J(ns+1,5) = (-5.0/2.0 )/ds**2 !
  J(ns  ,5) = ( 4.0/3.0 )/ds**2 !
  J(ns-1,5) = (-1.0/12.0)/ds**2 !

  J(2*ns+4,6) = 1.0 
!tout va bien
  !-------------7th eqn: seperated 1
  J(1,7) = T(0)*(theta(1) -theta(-1))/(2.0*ds**2) - p(0)/(2.0*ds) - (eta/(ds**4))*(-2.5*theta(0)+9.0*theta(1)-12.0*theta(2) +7.0*theta(3)-1.5*theta(4)) &
           +2.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4)
  J(2,7) = -2.5*eta*(theta(1)-theta(-1))/(ds**4) - 4.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4) &
           + 2.0*mue*(theta_dot_new(0)+alpha_t)*D
  J(3,7) = -T(0)*(theta(1) -theta(-1))/(2.0*ds**2) + p(0)/(2.0*ds) + (eta/(ds**4))*(-2.5*theta(0)+9.0*theta(1)-12.0*theta(2) +7.0*theta(3)-1.5*theta(4)) &
           +9.0*eta*(theta(1)-theta(-1))/(ds**4) + 2.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4)
  J(4,7) = -12.0*eta*(theta(1)-theta(-1))/(ds**4)
  J(5,7) =   7.0*eta*(theta(1)-theta(-1))/(ds**4)
  J(6,7) = - 1.5*eta*(theta(1)-theta(-1))/(ds**4)
  J(ns+4,7) = 1.0/(ds**2)
  J(ns+5,7) = (-2.0-((theta(1)-theta(-1))**2)/4.0)/(ds**2)
  J(ns+6,7) = 1.0/(ds**2)
  
  !--------------8th eqn
  J(2*ns+4,8) = (dt/(C1*mue*ds**2))* ((theta(ns)-2.0*theta(ns-1)+theta(ns-2)) + (3.0/2.0)*(theta(ns)-theta(ns-2)))
  J(2*ns+3,8) = -4.0*(dt/(C1*mue))*(theta(ns)-theta(ns-2))/(2.0*ds**2)
  J(2*ns+2,8) =  1.0*(dt/(C1*mue))*(theta(ns)-theta(ns-2))/(2.0*ds**2)
  J(ns+3,8)   = -dt*eta/((C1*mue)*ds**4)
  
  J(ns+2,8)   = (dt/(C1*mue*ds**2))   *   ( 4.0*eta/(ds**2) + (T(ns-1) + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2) &
	      + tau_beam(ns-1)*ds/(2.0) & 	! note we pulled 1/ds**2 out
              + eta*(theta(ns)-2.0*theta(ns-1)+theta(ns-2))*(theta(ns)-theta(ns-2))/(2.0*ds**2) &
              + (3.0*T(ns-1)-4.0*T(ns-2)+T(ns-3))/2.0         )
              
  J(ns  ,8)   = (dt/(C1*mue*ds**2))   *   ( 4.0*eta/(ds**2) + (T(ns-1) + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2) &
	      - tau_beam(ns-1)*ds/(2.0) &
              - eta*(theta(ns)-2.0*theta(ns-1)+theta(ns-2))*(theta(ns)-theta(ns-2))/(2.0*ds**2) &
              - (3.0*T(ns-1)-4.0*T(ns-2)+T(ns-3))/2.0         )

  J(ns+1,8)   = -D + (dt/(C1*mue))*( -6.0*eta/(ds**4) -2.0*(T(ns-1) + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2)/(ds**2) -D*sigma)              
              
  J(ns-1,8)   = -dt*eta/(C1*mue*ds**4)
  
  
!checked irreg points 16/01/11
  !-- all irregular points are filled. typing sucked. 
  do i=1, ns-2
    k = 8 + i !index for eqns F (row) -- that means first row is row 9
    l = i + 2 !index column (theta) -- starts @ theta_1 = index 3
    m = i + ns + 5 !index column (T) -- starts @ T(-1) = ns+5
    J( l ,k)  = -D + (dt/(C1*mue)) * (-6.0*eta/(ds**4) -(2.0/ds**2)*(T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 ) - D*sigma  )
    J( l-1,k) = (dt/(C1*mue))*( 4.0*eta/(ds**4) + (T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 )/(ds**2)  & 
		- tau_beam(i)/(2.0*ds) &
		- eta*(theta(i+1)-2.0*theta(i)+theta(i-1))*(theta(i+1)-theta(i-1))/(2.0*ds**4) - (T(i+1)-T(i-1))/(2.0*ds**2)  ) 
    J( l+1,k) = (dt/(C1*mue))*( 4.0*eta/(ds**4) + (T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 )/(ds**2)  & 
		+ tau_beam(i)/(2.0*ds) &
		+ eta*(theta(i+1)-2.0*theta(i)+theta(i-1))*(theta(i+1)-theta(i-1))/(2.0*ds**4) + (T(i+1)-T(i-1))/(2.0*ds**2)  ) 
    J( l-2,k) = -eta*dt/(C1*mue*ds**4) 
    J( l+2,k) = -eta*dt/(C1*mue*ds**4) 
    !block: dF/dT
    J( m,k)   =  (dt/(C1*mue))*(theta(i+1)-2.0*theta(i)+theta(i-1) )/(ds**2) 
    J( m-1,k) = -(dt/(C1*mue))*(theta(i+1)-theta(i-1) )/(2.0*ds**2) 
    J( m+1,k) =  (dt/(C1*mue))*(theta(i+1)-theta(i-1) )/(2.0*ds**2) 
  enddo
  do i=1, ns-2
    k = 8 + ns - 2 + i !index for eqns F (row) !starts @ ns+7
    l = i + 2 !index column (theta)
    m = i + ns + 5 !index column (T) -- starts @ T(-1) = ns+5
    !block dG/dTheta
    J( l ,k)  = -4.0*eta* (theta(i+1)-2.0*theta(i)+theta(i-1))/(ds**4) + 2.0*mue*(theta_dot_new(i)+alpha_t)*D
    
    J( l-1,k) = T(i)*(theta(i+1)-theta(i-1))/(2.0*ds**2) -p(i)/(2.0*ds) +eta*(theta(i+1)-theta(i-1))/(ds**4) & 
                 - eta* (theta(i+2)-2.0*theta(i+1)+2.0*theta(i-1)-theta(i-2))/(2.0*ds**4) +2.0*eta*(theta(i+1)-2.0*theta(i) +theta(i-1))/(ds**4) 
    J( l+1,k) = -T(i)*(theta(i+1)-theta(i-1))/(2.0*ds**2) +p(i)/(2.0*ds) -eta*(theta(i+1)-theta(i-1))/(ds**4) & 
                 + eta* (theta(i+2)-2.0*theta(i+1)+2.0*theta(i-1)-theta(i-2))/(2.0*ds**4) +2.0*eta*(theta(i+1)-2.0*theta(i) +theta(i-1))/(ds**4) 
    J( l-2,k) = -eta*(theta(i+1)-theta(i-1))/(2.0*ds**4) 
    J( l+2,k) =  eta*(theta(i+1)-theta(i-1))/(2.0*ds**4) 
    !block dG/dT
    J( m,k)   = -2.0/(ds**2) - ( (theta(i+1)-theta(i-1) )/(2.0*ds) )**2 
    J( m-1,k) = 1.0/(ds**2) 
    J( m+1,k) = 1.0/(ds**2) 
  enddo
  
  
  N_nonzero = 0
  do k = 1,2*ns+4
  do l = 1,2*ns+4
    !-------------------------------
    !-- count nonzero elements
    !-------------------------------
    if (J(k,l).ne.0.0) N_nonzero=N_nonzero+1
    !------------------------------------------
    ! check for NaN's
    !------------------------------------------
    if (isnan(J(k,l))) then
      write (*,*) "!!! !!! NaN in Jacobian..."
      stop
    endif  
  enddo
  enddo
  
  
  
end subroutine Jacobi




!##################################################################################################################################################################



! -------------------------------------------------------------------------------------------------------------------------------

subroutine GravityImpulse(time)
  use share_vars
  implicit none
  ! gives a little gravity impulse to pertubate the beam between T0 and T1 (sinusoidal)
  real (kind=pr), intent (in) ::  time
  real (kind=pr) :: T0,T1,a,b,c,d,k,t
  
  if (iImpulse==1) then
    T0=0.75
    T1=0.85
    grav = 1.00*sin ( pi*(time-T0)/(T1-T0) )
    
    if (time<T0) grav=0.0
    if (time>T1) grav=0.0
   elseif (iImpulse==2) then
    T0=0.75
    T1=0.85
    grav = 4.00*sin ( pi*(time-T0)/(T1-T0) )
    
    if (time<T0) grav=0.0
    if (time>T1) grav=0.0
  elseif (iImpulse==3) then
    T0=2.0
    T1=3.0
    if (time <= T0) then
      k = 0.0
    elseif ((time>=T0).and.(time<=T1)) then
      a = -20.0; b= 70.0; c=-84.0; d=35.0;
      t = (time-T0)/(T1-T0)
      k = a*t**7 + b*t**6 + c*t**5 + d*t**4
    elseif (time>T1) then
      k   = 1.0
    endif   
    grav = 0.7*k
  endif


end subroutine GravityImpulse

! -------------------------------------------------------------------------------------------------------------------------------




subroutine RHS_beameqn (time, theta, theta_dot, pressure_beam, T, tau_beam, beam_solid)
    !------------------------------------------------------------------------------
    ! Beam Equation right hand side at time_n
    ! Version 19.09.2012, completely debugged, gives exactly the same results as the matlab solver.
    ! INPUT
    ! time:		the time at which we compute the RHS. Note this is important: The leading edge motion 
    !			may be time dependent, (also in between a runge kutta step) so we have to call the subroutine
    !			at the right time
    ! pressure_beam	the pressure jump at time time
    ! OUTPUT:
    ! T			the tension in the beam at time time. we return it to form the complete initial guess for the implicit solvers
    !
    !------------------------------------------------------------------------------
    use mkl95_lapack
    use mkl95_precision
    use share_vars
    implicit none
    real (kind=pr), intent (in) :: 				time
    real (kind=pr) ::  						A1, A2, K2,C2
    real (kind=pr), dimension (0:ns-1), intent (in) :: 		pressure_beam, tau_beam
    real (kind=pr), dimension (0:ns-1), intent (out) :: 	T
    type(solid) :: beam_solid
    real (kind=pr), dimension (0:ns-1), intent (inout) :: 	theta, theta_dot
    real (kind=pr), dimension (0:ns+2) :: 			theta_extended, theta_extended_s, theta_extended_ss, theta_extended_sss, theta_extended_ssss
    real (kind=pr), dimension (0:ns-1) :: 			theta_s, theta_ss, theta_sss, theta_ssss, T_s, p_s
    real (kind=pr) :: 						alpha, alpha_t, alpha_tt
    real (kind=pr), dimension(1:6) :: 				LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  

    call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
    
    theta(0)     = 0.0 ! first boundary condition, angle here is predescribed by the motion protocol.
    theta_dot(0) = 0.0 ! modified: set to zero. was alpha, but is now angle in RELATIVE system (23.02.2011)
    
    !-- compute the tension in the beam
    call Tension (time, T, T_s, theta, theta_dot, pressure_beam, tau_beam, beam_solid)   
    
    ! -------------------------------------------------------------------------------------------------------
    ! Extend the beam with ghostpoints to fulfill the boundary conditions
    ! -------------------------------------------------------------------------------------------------------    
    !-- leading edge boundary conditions (constants)
    K2 = pressure_beam(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))
    C2 = 2.0*ds*K2 - T(0)*theta(1) + (eta/ds**2)*(10.0*theta(0)-12.0*theta(1)+6.0*theta(2)-theta(3) )

    !-- theta_extended(0) is the first virtual node. 
    theta_extended(0) 	 = C2 / ( (3.0*eta/ds**2)-T(0) )
    theta_extended(1:ns) = theta !because of this, theta_extended(1) = 0 (second point added for boundarys)

    !-- solve last two boundary conditions: theta_s(ns-1) = theta_ss(ns-1) = 0
    A1 = theta(ns-3)/12. - 2.*theta(ns-2)/3.
    A2 =-theta(ns-3)/12. + 4.*theta(ns-2)/3. - 5.*theta(ns-1)/2.
    !-- from these equations we can calculate the values of theta(ns) and theta(ns+1) which are added points
    theta_extended(ns+1) =  3. * (A1-A2)/2.
    theta_extended(ns+2) = 12. * (2.*A1-A2)

    !-- now we extended the grid by three points which are determined to fulfill the boundary conditions
    !-- dimension is now 0:ns+2 where points 1:ns correspond to the actual beam
    call Differentiate1D (theta_extended, theta_extended_s, ns+3, ds, 1)
    call Differentiate1D (theta_extended, theta_extended_ss, ns+3, ds, 2)
    call Differentiate1D (theta_extended, theta_extended_sss, ns+3, ds, 3)
    call Differentiate1D (theta_extended, theta_extended_ssss, ns+3, ds, 4)

    !-- cut the virtual nodes
    theta_s   =theta_extended_s(1:ns) 
    theta_ss  =theta_extended_ss(1:ns) 
    theta_sss =theta_extended_sss(1:ns) 
    theta_ssss=theta_extended_ssss(1:ns) 
    
    call Differentiate1D (pressure_beam, p_s, ns, ds, 1)

    ! ---------------------------------------------------------------------------------------------------------------------------------
    ! Last step: compute evolution equation for theta
    ! ---------------------------------------------------------------------------------------------------------------------------------
    theta     = theta_dot
    theta_dot = (-p_s - eta*theta_ssss + theta_ss*(T+eta*(theta_s**2)) + 2.*T_s*theta_s - mue*alpha_tt  - sigma*theta_dot + tau_beam*theta_s) / mue
    theta_dot(0) = 0.0  
 
end subroutine RHS_beameqn


! ----------------------------------------------------------------------------------------------------------------------------------------


subroutine Tension ( time, T, T_s, theta, theta_dot, pressure, tau_beam, beam_solid)
  ! note we use a TRIAG solver. the neumann condition is not in the matrix, it is set explicitly. therefore we return also T_s 
  use share_vars
  implicit none
  real (kind=pr), intent (in) ::				time
  real (kind=pr), dimension (0:ns-1) :: 			theta_s, theta_ss, theta_sss
  real (kind=pr), dimension (0:ns-1) :: 			tau_beam_s
  real (kind=pr), dimension (0:ns-1), intent (in) ::  		theta, theta_dot, pressure, tau_beam
  real (kind=pr), dimension (0:ns-1), intent (out) ::  		T, T_s
  type(solid) :: beam_solid
  real (kind=pr) ::						K1, T_s0
  real (kind=pr) :: 						alpha, alpha_t, alpha_tt
  real (kind=pr), dimension(1:6) :: 				LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  real (kind=pr), dimension (0:ns-1) :: 			diagonal_main
  real (kind=pr), dimension (0:ns-2) :: 			diagonal_up, diagonal_down
  real (kind=pr), dimension (0:ns-3) :: 			diagonal_up2
  real (kind=pr), dimension (0:ns-1) :: 			ipiv !used only for the MKL lib
  real (kind=pr), dimension (0:ns-1) ::				rhs
  integer :: 							i, info
  
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  
  call Differentiate1D (theta, theta_s, ns, ds, 1)
  call Differentiate1D (theta, theta_ss, ns, ds, 2)
  call Differentiate1D (theta, theta_sss, ns, ds, 3)
  call Differentiate1D (tau_beam, tau_beam_s, ns, ds, 1)
  
  
  ! constant of the BC, see masters thesis
  K1 = mue*(LeadingEdge(5)*cos(alpha)+LeadingEdge(6)*sin(alpha)+grav*sin(alpha)) - tau_beam(0)

  ! value of the inhomogenous Neumann BC at the leading edge
  T_s0 = -eta*theta_ss(0)*theta_s(0) + K1 

  ! the T-eqn with Neumann / Dirchlet conditions
  diagonal_main 	= -2.0/(ds**2) - theta_s**2
  diagonal_up   	=  1.0/(ds**2)
  diagonal_down 	=  1.0/(ds**2)  
  ! Dirichlet Condition
  diagonal_main(ns-1)	= 1.0 ! note this point is overwritten
  diagonal_down(ns-2) 	= 0.0 ! for the dirichlet condition   
  ! Neumann Condition
  diagonal_up(0) 	= 2.0/(ds**2)
  
  ! construct RHS for T-equation
  rhs(0:ns-1) 		= -(pressure*theta_s) -(2.*eta*theta_s*theta_sss)-(eta*theta_ss**2 ) -mue*(theta_dot + alpha_t)**2  - tau_beam_s
  ! Dirichlet
  rhs(ns-1) 		= 0.0
  ! Neumann
  rhs(0) 		= rhs(0) + T_s0 * 2.0/ds

  ! preconditioner - improve condition number of system by normalizing the diagonal
  do i=0,ns-1
    rhs(i) = rhs(i) / diagonal_main(i)
  enddo
  
  do i=0,ns-2
    diagonal_up(i)   = diagonal_up(i) / diagonal_main(i)
    diagonal_down(i) = diagonal_down(i) / diagonal_main(i+1)
  enddo  
  diagonal_main 	= 1.0


  ! solve system
  call dgttrf( ns, diagonal_down, diagonal_main, diagonal_up, diagonal_up2, ipiv, info )			! LU decomposition of matrix
  call dgttrs( 'N', ns, 1, diagonal_down, diagonal_main, diagonal_up, diagonal_up2, ipiv, rhs, ns, info )	! solve. now g is the solution  

  
  T = rhs  ! return tension
  
  call Differentiate1D ( T, T_s, ns, ds, 1)
  
  T_s (0) = T_s0 ! don't forget the Neumann Condition
 
  if (info.ne.0) then 
  write (*,*) '!!! MKL linear solver was experiencing trouble.', info
  stop
  endif
  
 
  
end subroutine Tension

! ----------------------------------------------------------------------------------------------------------------------------------------

subroutine Differentiate1D (f, f_derivative, N, dx, order)
  use share_vars
  implicit none
  integer, intent (in) 				:: order, n
  real, intent (in)    				:: dx
  real, dimension(0:N-1)			:: f, f_derivative
  real (kind=pr), dimension (0:N-1, 0:N-1) 	:: D1, D2, D3 ,D4
  
  call create_diff_matrices (D1, D2, D3, D4, N, dx)
  
  select case (order)
    case(1)
      f_derivative = matmul(D1,f)
    case(2)
      f_derivative = matmul(D2,f)
    case(3)
      f_derivative = matmul(D3,f)
    case(4)
      f_derivative = matmul(D4,f)
  end select
  
end subroutine Differentiate1D
  
! -------------------------------------------------------------------------------------------------------------------------------

subroutine create_diff_matrices (D1, D2, D3, D4, nn, dx)
  use share_vars
  implicit none
  integer :: i,N
  integer, intent (in) :: nn
  real(kind=pr), intent(in) :: dx
  real (kind=pr), dimension (0:nn-1, 0:nn-1), intent (out) :: D1, D2, D3 ,D4

  !------------------------------------
  ! create diff-matrices.
  ! first derivative: x'=D1 * x
  !------------------------------------

  D1=0.0
  D2=0.0
  D3=0.0
  D4=0.0

  do i=0,nn-1 !main diag
    D1(i,i) = 0. !main
    D2(i,i) =-2. !main
    D3(i,i) = 0. !main
    D4(i,i) = 6. !main   
  enddo
 
  do i=0,nn-2 !first upper/lower diag
    D1(i,i+1) = 1.  !up1
    D1(i+1,i) =-1.  !lo1
    D2(i,i+1) = 1.  !up1
    D2(i+1,i) = 1.  !lo1
    D3(i,i+1) =-2.  !up1
    D3(i+1,i) = 2.  !lo1
    D4(i,i+1) =-4.  !up1
    D4(i+1,i) =-4.  !lo1
  enddo

  do i=0,nn-3 !second upper/lower diag
    D3(i,i+2) = 1.   !up2
    D3(i+2,i) =-1.   !lo2
    D4(i,i+2) = 1.   !up2
    D4(i+2,i) = 1.   !lo2
  enddo
!   N=ns+2 !dimension of the matrices (we added 3 points for 3 boundary equations
    N = nn-1
!-----------------------------  
    D1(0,0)=  -3. 
    D1(0,1)=   4. 
    D1(0,2)=  -1. 
    D1(N,N-2)= 1. 
    D1(N,N-1)=-4. 
    D1(N,N)=   3. 
    D1=D1/(2.0* dx)     
!-----------------------------
    D2(0,0)=2. 
    D2(0,1)=-5. 
    D2(0,2)=4. 
    D2(0,3)=-1.     
    D2(N,N-3)=-1. 
    D2(N,N-2)=4. 
    D2(N,N-1)=-5. 
    D2(N,N)=2.     
    D2=D2/( dx**2 )     
!-----------------------------    
    D3(0,0)=-17.0/2. 
    D3(0,1)=71.0/2. 
    D3(0,2)=-59. 
    D3(0,3)=49. 
    D3(0,4)=-41.0/2.     
    D3(0,5)=7.0/2.   
    D3(1,0)=0. 
    D3(1,1)=-17.0/2. 
    D3(1,2)=71.0/2. 
    D3(1,3)=-59. 
    D3(1,4)=49. 
    D3(1,5)=-41.0/2.     
    D3(1,6)=7.0/2.  
    D3(N,N-5)=-7.0/2. 
    D3(N,N-4)=41.0/2. 
    D3(N,N-3)=-49. 
    D3(N,N-2)=59. 
    D3(N,N-1)=-71.0/2. 
    D3(N,N)=17.0/2.     
    D3(N-1,N-6)=-7.0/2. 
    D3(N-1,N-5)=41.0/2. 
    D3(N-1,N-4)=-49.0 
    D3(N-1,N-3)=59.0 
    D3(N-1,N-2)=-71.0/2. 
    D3(N-1,N-1)=17.0/2. 
    D3(N-1,N)=0.     
    D3=D3/(2.0* dx**3)     
!-----------------------------
    D4(0,0)=35.0/6. 
    D4(0,1)=-31. 
    D4(0,2)=137.0/2. 
    D4(0,3)=-242.0/3. 
    D4(0,4)=107.0/2.     
    D4(0,5)=-19.0 
    D4(0,6)=17.0/6.     
    D4(1,0)=0. 
    D4(1,1)=35.0/6. 
    D4(1,2)=-31. 
    D4(1,3)=137.0/2. 
    D4(1,4)=-242.0/3.0 
    D4(1,5)=107.0/2.     
    D4(1,6)=-19.0 
    D4(1,7)=17.0/6.0     
    D4(N,N-6)=17.0/6.0  
    D4(N,N-5)=-19.0 
    D4(N,N-4)=107.0/2.0 
    D4(N,N-3)=-242.0/3.0 
    D4(N,N-2)=137.0/2.0 
    D4(N,N-1)=-31.0 
    D4(N,N)=35.0/6.0 
    D4(N-1,N-7)=17.0/6. 
    D4(N-1,N-6)=-19. 
    D4(N-1,N-5)=+107.0/2. 
    D4(N-1,N-4)=-242.0/3. 
    D4(N-1,N-3)=+137.0/2. 
    D4(N-1,N-2)=-31.0 
    D4(N-1,N-1)=+35.0/6. 
    D4(N-1,N)=0.         
    D4=D4/( dx**4)     
 !-----------------------------
end subroutine create_diff_matrices

!##################################################################################################################################################################

subroutine RK4 (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid) 
use share_vars
  implicit none
  real (kind=pr), intent (in) :: dt_beam, time
  real (kind=pr), dimension (0:ns-1, 1:6), intent (inout) :: beam
  real (kind=pr), dimension (0:ns-1, 1:6) :: beam_old
  real (kind=pr), dimension (0:ns-1), intent (in) :: pressure_beam, tau_beam
  type (solid) :: beam_solid
  real (kind=pr), dimension (0:ns-1), intent (inout) :: T
  real (kind=pr), dimension (0:ns-1) :: theta, theta_dot, T1, T2, T3, T4
  real (kind=pr), dimension (0:ns-1) :: theta_1, theta_2, theta_3, theta_4
  real (kind=pr), dimension (0:ns-1) :: theta_dot_1, theta_dot_2, theta_dot_3, theta_dot_4
  beam_old=beam
  theta = beam(:,5)
  theta_dot = beam(:,6)

  theta_1 = 0.0
  theta_2 = 0.0
  theta_3 = 0.0
  theta_4 = 0.0
  theta_dot_1 = 0.0
  theta_dot_2 = 0.0
  theta_dot_3 = 0.0
  theta_dot_4 = 0.0
  
  T1 = T
  T2 = T
  T3 = T
  T4 = T
  
 
 !is just a runge-kutta 4th order
 !subroutine RHS_beameqn (time, theta , theta_dot, pressure_beam )
    theta_1 = theta
    theta_dot_1 = theta_dot
  call RHS_beameqn (time, theta_1 , theta_dot_1, pressure_beam, T1, tau_beam, beam_solid)
    theta_2     = theta + 0.5*dt_beam*theta_1
    theta_dot_2 = theta_dot + 0.5*dt_beam*theta_dot_1
  call RHS_beameqn (time+0.5*dt_beam, theta_2 , theta_dot_2, pressure_beam, T2, tau_beam, beam_solid)
    theta_3     = theta + 0.5*dt_beam*theta_2
    theta_dot_3 = theta_dot + 0.5*dt_beam*theta_dot_2
  call RHS_beameqn (time+0.5*dt_beam, theta_3 , theta_dot_3, pressure_beam, T3, tau_beam, beam_solid)
    theta_4     = theta + dt_beam*theta_3
    theta_dot_4 = theta_dot + dt_beam*theta_dot_3
  call RHS_beameqn (time+dt_beam, theta_4 , theta_dot_4, pressure_beam, T4, tau_beam, beam_solid)
    theta       = theta + dt_beam * (theta_1 + 2.*theta_2 + 2.*theta_3 + theta_4 )/6.0
    theta_dot   = theta_dot + dt_beam * (theta_dot_1 + 2.*theta_dot_2 + 2.*theta_dot_3 + theta_dot_4 )/6.0
    
  beam(:,5) = theta
  beam(:,6) = theta_dot
  
  T=T1

  !---------------------------------------------------------------------
  !     emergency brake
  !---------------------------------------------------------------------
  if (maxval(abs(beam(:,6)-beam_old(:,6) ))>100.0 ) then
    write (*,'(A)') "??????????????????????????????????????????????????????????????????????????????????????????????????????"
    call DisplayError()
    write (*,'(A)') " "
    write (*,'(A)') "!!! rk4-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.0"
    write (*,'(A)') "    That indicates a problem, probably artificial added mass instability."
    write (*,'(A)') "    Even though it won't make sense to continue, I make a backup, quit, and let you decide"   
    write (*,'("time=",es11.4, " dt=",es11.4)') time, dt_beam
    write (*,'(A)') "??????????????????????????????????????????????????????????????????????????????????????????????????????" 
    continue_timestep = .false.
  endif  
  
end subroutine RK4

! ---------------------------------------------------------------------------------------------------------------------------------------
  
subroutine EE1 (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid) 
use share_vars
  implicit none
  real (kind=pr), intent (in) :: dt_beam, time
  real (kind=pr), dimension (0:ns-1, 1:6), intent (inout) :: beam
  real (kind=pr), dimension (0:ns-1), intent(out) :: T
  real (kind=pr), dimension (0:ns-1), intent (in) :: pressure_beam, tau_beam
  type (solid) :: beam_solid
  real (kind=pr), dimension (0:ns-1) :: theta, theta_dot
  real (kind=pr), dimension (0:ns-1) :: theta_1
  real (kind=pr), dimension (0:ns-1) :: theta_dot_1

  theta     = beam(:,5)
  theta_dot = beam(:,6)

  theta_1     = theta
  theta_dot_1 = theta_dot
  
  call RHS_beameqn (time, theta_1 , theta_dot_1, pressure_beam, T, tau_beam, beam_solid)

  beam(:,5) = theta + dt_beam * theta_1
  beam(:,6) = theta_dot + dt_beam * theta_dot_1
end subroutine EE1


!##################################################################################################################################################################

subroutine beamfilter(time,theta)
  ! note this routine uses ONE BASED INDEXING: its easier because its developped in MATLAB
  use share_vars
  implicit none
  real (kind=pr), dimension (1:ns), intent(inout) :: theta
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension (1:ns) :: theta_old, theta_filtered
  real (kind=pr), dimension (:), allocatable :: theta_extended, tukey, theta_extended_k, spectrum
  integer :: i, n_ext, n_tukey, k, k1,k2
  real(kind=pr) :: filter
  

  
!   open (14, file = 'beam.in', status='old', action='read') ! Append output data file
!   do i=1,ns
!   read (14, *) theta(i)
!   enddo
!   close (14) 
  

  ! ---------------------------------------
  ! first, we extend the beam 
  ! ---------------------------------------
  n_ext = ns-12 ! n_ext is the number of points we added
  
  allocate ( tukey(1:ns+2*n_ext),theta_extended(1:ns+2*n_ext),theta_extended_k(1:ns+2*n_ext),spectrum(1:ns+2*n_ext) )
  
  ! the first part is antisymmetric: its the same values with a neg sign in reverse order
  theta_extended (1 : n_ext)               = -theta(n_ext+1:2:-1)
  ! embedd original beam
  theta_extended (n_ext+1 : ns+n_ext)      =  theta(1 : ns) 
  ! the last part is a mirror of the beam: note that due to BC, theta_s and theta_ss are zero
  theta_extended (ns+n_ext+1 : ns+2*n_ext) =  theta(ns-1 : (ns-1)-(n_ext-1) : -1)  ! skip point ns here, otherwise taken twice
  


  ! ---------------------------------------
  ! tukey window
  ! ---------------------------------------  
  n_tukey = 3*n_ext/4
  tukey = 1.0
  
  do i=1,n_tukey
    tukey(i) = 0.5*(1.0 - cos( pi*real(i-1)/real(n_tukey) ) )
    tukey(ns+2*n_ext-i+1) = tukey(i)
  enddo
  
  
  ! appley Tukey window to extended beam
  theta_extended = theta_extended * tukey
  
  ! go to fourier space
  call cofts (theta_extended, theta_extended_k, ns+2*n_ext, 1)
  
  spectrum=theta_extended_k
  
! write(*,'("theta",1600(es11.4,1x))') theta  
! write(*,'("theta_ext",1600(es11.4,1x))') theta_extended 
! write(*,'("tukey",1600(es11.4,1x))') tukey  
! write(*,'("theta_ex_k ",1600(es11.4,1x))') theta_extended_k
  
  ! apply filter
  k1=100
  k2=200
  do k = 1 , (ns+2*n_ext)/2  ! is k_max=756 nodes for ns=512

    if (k<k1) then
      filter = 1.0
    elseif ((k>=k1).and.(k<=k2)) then
      filter = 1.0 - 0.5*(1.0 - cos( pi*real(k-k1)/real(k2-k1) ) )
    else
      filter = 0.0
    endif
    
    theta_extended_k (2*k-1) = filter*theta_extended_k (2*k-1) !real part
    theta_extended_k (2*k)   = filter*theta_extended_k (2*k) !imag part
    
  end do

  
  ! inverse fourier transform
  call cofits (theta_extended_k, theta_extended, ns+2*n_ext, 1)
  
  ! cut added values
  theta_old = theta
  theta_filtered = theta_extended (n_ext+1 : ns+n_ext)
  
  !===================================================================================
!   blender = 1.0
! 
! !   do i=5,10
! !     blender(i) = real(i-5)/5.0
! !   enddo  
! !   blender( ns-10:ns-5 ) = blender( 10:5:-1 )
! 
!   n_blender = 20
!   n_BC = 4
!     blender(1:n_BC) = 0.0
!   blender(ns-n_BC:ns) = 0.0
!   do i=n_BC,n_blender+n_BC
!     blender(i) = 0.5*(1.0 - cos( pi*real(i-n_BC)/real(n_blender) ) )
!     blender(ns-i+1) = blender(i)
!   enddo
!   theta = blender*theta_filtered +(1.0-blender)*theta_old
  !===================================================================================
  
  theta = theta_filtered

!   write(*,'("BLENDER",1600(es11.4,1x))') blender
!   stop
  
! write(*,'("theta_k_filt",1600(es11.4,1x))') theta_extended_k  
! write(*,'("theta_neu",1600(es11.4,1x))') theta
! write(*,*) "---"
write(*,*) maxval(abs(theta-theta_old))

!   open (14, file = 'beam_spectrum', status = 'unknown', access = 'append') ! Append output data file
!   write (14, '(1600(es11.4,1x))') time, (sqrt(spectrum(2*k-1)**2+spectrum(2*k)**2), k = 1,(ns+2*n_ext)/2 )
!   close (14)

end subroutine

!##################################################################################################################################################################


subroutine Jacobi_num(time, dt,dt_old, J, N_nonzero , T, theta, theta_old, theta_dot_old, p, theta_oldold, theta_dot_oldold, old_rhs, tau_beam_new, beam_solid)
  use share_vars
  implicit none
  type (solid) :: beam_solid
  real (kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real (kind=pr), dimension(1:2*ns+4) :: F1,F2
  real (kind=pr), dimension(-1:ns-1), intent (in) :: T
  real (kind=pr), dimension(-1:ns-1) :: T1, T2
  real (kind=pr), intent (in) :: time, dt,dt_old
  real (kind=pr), dimension(-1:ns+1), intent (in) :: theta
  real (kind=pr), dimension(-1:ns+1) :: theta1, theta2
  real (kind=pr), dimension(0:ns-1), intent (in) :: p, theta_old, theta_dot_old, theta_oldold,theta_dot_oldold, old_rhs, tau_beam_new
  integer, intent (out) :: N_nonzero
  integer :: i, k, l

  !---------------------------------------
  !	NUMERIC COMPUTATION OF THE JACOBIAN
  !	-> I think the row/column ordering in J (in IBES) is ackward; pay attention.
  !---------------------------------------
 
  ! derive wrt to theta_new
  l=1
  do i=-1,ns+1
      theta1    = theta
      theta1(i) = theta1(i)-1.0e-6
      theta2    = theta
      theta2(i) = theta2(i)+1.0e-6
      
      call F_nonlinear (time, dt, dt_old, F1, theta_old, theta_dot_old, theta1, T, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid )
      call F_nonlinear (time, dt, dt_old, F2, theta_old, theta_dot_old, theta2, T, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid )
      
      J(l,:)= (F2-F1)/(2.0e-6)
      l=l+1

  enddo

  ! derive wrt to T_new
  do i=-1,ns-1
      T1=T
      T1(i) = T1(i)-1.0e-6
      T2=T
      T2(i) = T2(i)+1.0e-6
      
      call F_nonlinear (time, dt, dt_old, F1, theta_old, theta_dot_old, theta, T1, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid  )
      call F_nonlinear (time, dt, dt_old, F2, theta_old, theta_dot_old, theta, T2, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid  )
      
      J(l,:)= (F2-F1)/(2.0e-6)
      l=l+1
  enddo

    !-- count nonzero elements
  N_nonzero = 0
  do k = 1,2*ns+4
  do l = 1,2*ns+4
    if (J(k,l).ne.0.0) N_nonzero=N_nonzero+1
  enddo
  enddo
  
end subroutine Jacobi_num


subroutine Check_Vector_NAN(f, msg)
  use share_vars
  real(kind=pr), intent(in) :: f(:)
  character(len=*) :: msg
  integer :: a, i
  a = size(f)
  do i=1,a
    if (isnan(f(i))) then
     write (*,*) "??? SOLID SOLVER: Found NaN in vector at", i
     write (*,*) msg
     stop
     endif
  enddo
end subroutine

subroutine Check_Vector_NAN_try_correct(f, msg)
  use share_vars
  real(kind=pr), intent(inout) :: f(:)
  character(len=*) :: msg
  integer :: a, i
  a = size(f)
  do i=1,a
    if (isnan(f(i))) then
     write (*,*) "??? SOLID SOLVER: Found NaN in vector at", i
     write (*,*) msg
     write (*,*) "I will try to correct this and proceed, with fingers crossed!"
     f(i)=0.0
     endif
  enddo
end subroutine
 
end module



