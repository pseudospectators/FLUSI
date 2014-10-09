
module solid_model
  use vars ! for the precision statement
  use basic_operators 
  implicit none
  
  !----------------------------------------------
  ! module global variables
  !----------------------------------------------  
  integer,parameter :: nBeams = 1  
  integer,save :: ns
  ! see "type solid" about nsmax 06 Aug 2014
  integer,parameter :: nsmax = 200
  ! TODO: move these into the solid model datastructure
  real(kind=pr),save :: mue
  real(kind=pr),save :: eta
  real(kind=pr),save :: grav
  real(kind=pr),save :: sigma
  real(kind=pr),save :: t_beam, L_span, N_smooth
  real(kind=pr),save :: AngleBeam, frequ
  real(kind=pr),save :: ds
  real(kind=pr),save :: T_release, tau
  real(kind=pr),save :: R_cylinder
  real(kind=pr), parameter :: error_stop = 1.0e-6  
  character(len=strlen),save :: imposed_motion_leadingedge, TimeMethodSolid
  character(len=strlen),save :: has_cylinder
  character(len=strlen),save :: interp
  character(len=strlen),save :: infinite, plate_shape

  ! this is a hack to avoid allocating/deallacting these arrays in every time
  ! step. turing fails with memory issues if done otherwise. 7 Aug 2014
  real(kind=pr),save,dimension(:,:,:,:), allocatable :: surfaces
  real(kind=pr),save,dimension(:,:,:), allocatable :: p_surface, p_surface_local  
  real(kind=pr),save,dimension(:,:), allocatable :: active_points
  real(kind=pr),save,dimension(:),allocatable :: heights
  
  !----------------------------------------------
  ! Solid datatype
  !----------------------------------------------
  type Solid
    ! since 6 Aug 2014, these arrays are statically allocated (they lie thus on 
    ! the stack), since on turing we had problems with memory fragmentation. It
    ! is still possible to use ns<nsmax via the params file, but not ns>nsmax.
    real(kind=pr),dimension(0:nsmax) :: x,y,vx,vy
    real(kind=pr),dimension(0:nsmax) :: theta,theta_dot,ax,ay
    real(kind=pr),dimension(0:nsmax) :: pressure_old, pressure_new
    real(kind=pr),dimension(0:nsmax) :: tau_old, tau_new
    real(kind=pr),dimension(0:nsmax,1:6) :: beam_oldold
    real(kind=pr),dimension(1:2) :: Force, Force_unst, Force_press, Inertial_Force
    real(kind=pr) :: E_kinetic, E_pot, E_elastic
    real(kind=pr) :: x0, y0, AngleBeam, phase
    ! we need the previous time step for the BDF solver
    real(kind=pr) :: dt_old    
    ! these replace the save variables in the unst correction computation:
    real(kind=pr) :: drag_unst_new, drag_unst_old, lift_unst_new, lift_unst_old
    logical :: StartupStep, UnsteadyCorrectionsReady
  end type Solid
 
 
  
 contains
 
 
 !-----------------------------------------------------------------
 include "mouvement.f90"
 include "integrate_position.f90"
 include "init_beam.f90"
 include "save_beam.f90"
 include "BeamForces.f90"
 include "plate_geometry.f90"
 
!-------------------------------------------------------------------------------
!   solid solver main entry point
!-------------------------------------------------------------------------------
subroutine OnlySolidSimulation()
  use vars
  implicit none
  type(solid), dimension(1:nBeams) :: beams
  real (kind=pr) :: time
  integer :: it,nsave

  ! in case you forgot to set it
  if (itdrag==0) itdrag=10
  if (dt_fixed<=1.0d-10) dt_fixed=1.d-3
  
  write (*,*) "*** information: starting OnlySolidSimulation"
  time = 0.0
  it = 0  
  
  call show_solid_model_information
  
  !-- initialization
  call init_beams( beams )

  !--loop over time steps
  do while ((time<=tmax))
    !-- time stepping
    call SolidSolverWrapper( time, dt_fixed , beams )
        
    it   = it+1
    time = dble(it)*dt_fixed
    
    call SaveBeamData( time, beams )
  enddo

end subroutine
 
 
 
 
!-------------------------------------------------------------------------------
!   SOLID SOLVER WRAPPER
!-------------------------------------------------------------------------------
subroutine SolidSolverWrapper ( time, dt, beams )
  implicit none
  real(kind=pr), intent (in) ::  dt, time
  real(kind=pr) :: t0
  type(solid), dimension(1:nBeams), intent (inout) ::    beams
  integer :: i
  t0 = MPI_wtime()
  
  do i = 1, nBeams
      !------------------------------------------
      ! check if input values are okay
      !------------------------------------------
      if (Vector_isNAN(beams(i)%pressure_new).or.&
          Vector_isNAN(beams(i)%pressure_old).or.&
          Vector_isNAN(beams(i)%tau_new).or.&
          Vector_isNAN(beams(i)%tau_old) ) then
        if (root) write(*,*) "SolidSolver: input values contain NaNs"
        ! time to go..
        call abort()
      endif
      
      if (time>=T_release) then 
        !-------------------------------------------
        ! the beams are released, call IBES solvers
        !-------------------------------------------
        ! all implicit solvers are in one subroutine
        select case (TimeMethodSolid)
          case ("CN2","BDF2","EI1")
            call IBES_solver (time, dt, beams(i))    
          case ("RK4")
            call RK4_wrapper (time, dt, beams(i))
          case ("EE1")
            call EE1_wrapper (time, dt, beams(i))
          case default
            write(*,*) "SolidSolverWrapper::invalid value of TimeMethodSolid",&
                TimeMethodSolid
            call abort()
        end select
      else 
        !-------------------------------------------
        ! the beams are not yet released, but their leading edges may move
        !-------------------------------------------
        call integrate_position (time+dt, beams(i))
      endif
      
      !-------------------------------------------
      ! compute energies and stuff
      !-------------------------------------------
      call SolidEnergies( beams(i) )

      !-- check if everything seems okay, if not show beam and abort
      call show_beam_on_error( beams(i) )
  enddo
  time_solid = time_solid + MPI_wtime() - t0
end subroutine SolidSolverWrapper




!-------------------------------------------------------------------------------
!   SOLID SOLVER INITIALIZATION
!-------------------------------------------------------------------------------
subroutine InitializeSolidSolver( beams )
  implicit none
  integer :: i
  type(solid), dimension(1:nBeams), intent (inout) :: beams
  
  ! marks all beams to be in the very first time step  
  ! the solver then uses CN2 instead of BDF2, because the old old time level
  ! t_n-1 is not available
  do i=1,nBeams
    beams(i)%StartupStep = .true.
  enddo
  
  
end subroutine InitializeSolidSolver





!-------------------------------------------------------------------------------
!   energies and stuff for beams
!-------------------------------------------------------------------------------
subroutine SolidEnergies( beam )
  implicit none
  type(solid), intent (inout) :: beam
  real(kind=pr),dimension(0:ns-1) :: theta_s
  
  ! for beam elastic energy, we need theta_s over the beam 
  call Differentiate1D ( beam%theta(0:ns-1), theta_s, ns, ds, 1)  
  
  beam%E_kinetic = mue*0.5d0*ds*sum( beam%vx(0:ns-1)**2 + beam%vy(0:ns-1)**2 )
  beam%E_pot     = mue*grav*ds *sum( beam%y(0:ns-1)-beam%y0 )
  beam%E_elastic = eta*0.5d0*ds*sum( theta_s(0:ns-1)**2 )
  
  beam%Inertial_Force(1) = mue*ds* sum( beam%ax(0:ns-1) ) 
  beam%Inertial_Force(2) = mue*ds* sum( beam%ay(0:ns-1) ) 
  
end subroutine SolidEnergies



!-------------------------------------------------------------------------------
!   SOLID SOLVER ROUTINES
!-------------------------------------------------------------------------------
subroutine IBES_solver ( time, dt, beam_solid )! note this is actuall only ONE beam!!
  implicit none
  real(kind=pr), intent (in) :: dt, time
  type (solid), intent(inout) :: beam_solid
  real(kind=pr),dimension(0:ns-1) :: T
  real(kind=pr),dimension(0:ns-1) :: old_rhs, theta_old,T_dummy
  real(kind=pr),dimension(0:ns-1) :: pressure_old, pressure_new, tau_beam_old, tau_beam_new
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam
  real(kind=pr),dimension(-1:ns+1) :: theta_guess
  real(kind=pr),dimension(-1:ns-1) :: T_guess
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam_old, beam_guess
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam_oldold
  real(kind=pr),dimension(1:2*ns+4) :: x_guess, x, x_delta, F
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4) :: J,J2,J2_norm
  real(kind=pr),dimension(1:ns+3) :: theta_act
  real(kind=pr),dimension(1:ns+1) :: T_act
  real(kind=pr) :: alpha, alpha_t, alpha_tt, err, A1, A2, K1, K2, T_s0, theta_ss0, theta_s0, C2,C1,C3,C4
  real(kind=pr) :: err_rel, R
  real(kind=pr) :: dt_old  
  integer :: n,iter
  integer, save :: iCalls=10
  logical :: ActuallyBDF2=.false., iterate=.true.
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  
  T_dummy = 0.0
  
  !******************************************************* 
  ! NOTE: 2013, this is extended to take more than one beam into account.
  ! however, its too hard to reprogram everything, so now
  ! we force it to be compatible
  !*******************************************************
  beam(0:ns-1,1) = beam_solid%x(0:ns-1)
  beam(0:ns-1,2) = beam_solid%y(0:ns-1)
  beam(0:ns-1,3) = beam_solid%vx(0:ns-1)
  beam(0:ns-1,4) = beam_solid%vy(0:ns-1)
  beam(0:ns-1,5) = beam_solid%theta(0:ns-1)
  beam(0:ns-1,6) = beam_solid%theta_dot(0:ns-1)
  pressure_old(0:ns-1) = beam_solid%pressure_old(0:ns-1)
  pressure_new(0:ns-1) = beam_solid%pressure_new(0:ns-1)
  tau_beam_old(0:ns-1) = beam_solid%tau_old(0:ns-1)
  tau_beam_new(0:ns-1) = beam_solid%tau_new(0:ns-1)
  dt_old       = beam_solid%dt_old
  beam_oldold(0:ns-1,1:6)  = beam_solid%beam_oldold(0:ns-1,1:6)
  !*******************************************************
  
  call mouvement ( time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
    
  beam_old = beam
  beam_guess = beam  

  !--------------------------------------------------------------------
  !   Startup time step: use CN2 as first step of BDF2 (pay attention to call InitializeSolidSolver() before the run)
  !-------------------------------------------------------------------- 
  if (beam_solid%StartupStep) then  ! this is the first time step 
    beam_solid%StartupStep = .false.  ! we're about to do the first step
    if (TimeMethodSolid=="BDF2") then  ! if we deal with BDF2
      ActuallyBDF2 = .true.    ! Remember to switch back to BDF2 (at the end of the step)
      TimeMethodSolid = "CN2"  ! use "CN2" for the first step
    endif
  endif
  

  !--------------------------------------------------------------------
  !   Initial guess for the beam at the new time level. 
  !   An euler-explicit step is made, ghostpoints are added 
  !--------------------------------------------------------------------
  
  call EE1s (time, dt, beam_guess, pressure_old, T, tau_beam_old, beam_solid)
  ! now we compute the ghostpoints Theta(-1); Theta(ns); Theta(ns+1); T(-1)
  theta_guess(0:ns-1) = beam_guess(0:ns-1,5) !because of this, temp(1) = 0 (second point added for boundarys)
  !-- leading edge boundary conditions (constants)
  K2 = pressure_new(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))
  C2 = 2.0*ds*K2 - T(0)*theta_guess(1) + (eta/ds**2)*(10.d0*theta_guess(0)-12.0*theta_guess(1)+6.0*theta_guess(2)-theta_guess(3) )
  !-- theta_extended(0) is the first virtual node. 
  theta_guess(-1) = C2 / ( (3.0*eta/ds**2)-T(0) )
  !-- solve last two boundary conditions: theta_s(ns-1) = theta_ss(ns-1) = 0
  A1 = beam_guess(ns-3,5)/12. - 2.*beam_guess(ns-2,5)/3.
  A2 =-beam_guess(ns-3,5)/12. + 4.*beam_guess(ns-2,5)/3. - 5.*beam_guess(ns-1,5)/2.
  !-- from these equations we can calculate the values of theta(ns) and theta(ns+1) which are added points
  theta_guess(ns) =  3. * (A1-A2)/2.
  theta_guess(ns+1) = 12. * (2.*A1-A2)
  

  !-- set virtual node for Tension T
  ! constant of the BC, see masters thesis
  K1 = mue*(LeadingEdge(5)*cos(alpha)+LeadingEdge(6)*sin(alpha)+grav*sin(alpha))  - tau_beam_new(0) !technically, the guess is at the new time level
  ! value of the inhomogenous Neumann BC at the leading edge
  theta_ss0 = ( 2.0*theta_guess(0) - 5.0*theta_guess(1) +4.0*theta_guess(2) - 1.0*theta_guess(3) )/ds**2
  theta_s0  = (-1.5*theta_guess(0) + 2.0*theta_guess(1) -0.5*theta_guess(2)  )/ds  
  T_s0 = -eta * theta_ss0*theta_s0 + K1 
  T_guess(-1) = T(1) - 2.0*ds*T_s0 !T_s = 0 means T(-1) = T(1)  
  T_guess(0:ns-1) = T(0:ns-1)  
  
  
  !-- rearrange the guess in the vector of unknowns x
  x_guess(1:ns+3)      = theta_guess(-1:ns+1)
  x_guess(ns+4:2*ns+4) = T_guess(-1:ns-1)  
  !now we have a complete initial guess for the variables

  !--------------------------------------------------------------------
  !   Calculate RHS vector @ t_n. This is required for the CN2 method only.
  !--------------------------------------------------------------------
  if (TimeMethodSolid == "CN2") then
    theta_old = beam_old(:,5)
    old_rhs = beam_old(:,6)           
    call RHS_beameqn(time, theta_old, old_rhs, pressure_old, T_dummy, tau_beam_old, beam_solid)
  else
    theta_old = 0.0
    old_rhs=0.0
  endif  
  !--------------------------------------------------------------------
  !   Newton iteration
  !--------------------------------------------------------------------
  err     = 1.0
  err_rel = 1.0
  
  x = x_guess
    
  iter = 0
  
  ! MODIFICATION 31.10.2012: the iteration now uses both relative and absolute error. If x is large, we can have trouble
  ! reaching a very small increment, and iterate forever. If x is small, then the relative criterion tries to go much below
  ! machine precision. So we use both, either with the same precision.
  iterate = .true.
  do while (iterate)
    theta_act = x(1:ns+3)
    T_act     = x(ns+4:2*ns+4)
    !-------------------------------------------------------------------------
    !  Calculate RHS vector
    !-------------------------------------------------------------------------
    call F_nonlinear( time, dt, dt_old, F, beam_old(:,5), beam_old(:,6), theta_act, T_act, pressure_new,&
                      old_rhs, beam_oldold(:,5), beam_oldold(:,6), tau_beam_new, beam_solid)
    F = -1.0*F !newton raphson is J*dx = -F

    !-------------------------------------------------------------------------
    !  Create Jacobi Matrix and Compress it to the sparse format.
    !-------------------------------------------------------------------------    
    call Jacobi(time, dt, dt_old, J, T_act, theta_act, beam_old(:,5), beam_old(:,6), pressure_new,&
                beam_oldold(:,5), beam_oldold(:,6), tau_beam_new, beam_solid)

    !-------------------------------------------------------------------------
    !  self-test. occasionally, check if Jacobian is okay
    !-------------------------------------------------------------------------     
    if (mod(iCalls,4607)==0) then
      call Jacobi_num(time, dt,dt_old, J2, T_act, theta_act, beam_old(:,5), beam_old(:,6),&
                      pressure_new, beam_oldold(:,5), beam_oldold(:,6), old_rhs, tau_beam_new, beam_solid)
      J2_norm = J2
      where (abs(J2_norm)<1.0e-7) J2_norm=1.0
      J2_norm = (J-J2)/J2_norm
      where (abs(J2_norm)<1.0e-5) J2_norm=0.d0 ! delete small values
      if (maxval(J2_norm)>1.0e-2) then  
          open (14, file = 'IBES_JACOBIAN', status = 'replace') ! Append output data file
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
          write(*,'(A)') "IBES error, mismatch in Jacobian (numeric vs analytic)"
          write (*,*) time, iter
          close (14)
          call abort()
      endif
      iCalls = 0
    endif
    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! SOLVE LINEAR SYSTEM
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call Solve_LGS ( J, F, x_delta )
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    iter    = iter + 1
    x       = x + x_delta    
    err     = sqrt(sum(x_delta**2))
    err_rel = abs(sqrt(sum(x_delta**2)) / sqrt(sum(x**2)))

    !-----------------------------------------------------------
    ! CONVERGENCE CRITERION
    !-----------------------------------------------------------
    if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
      iterate = .false.
    endif
    
    !-----------------------------------------------------------
    ! EMERGENCY BRAKE
    !-----------------------------------------------------------
    if ((iter>1000).and.(root)) then      
      write(*,*) "!!! ERROR: IBES performed like 1000 iterations. this is not normal. time=", time
      call abort()
    endif
    
  enddo
  
  !--------------------------------------------------------------------
  !    cut ghostpoints 
  !--------------------------------------------------------------------
  beam(:,5)=x(2:ns+1)
  T = x(ns+5:2*ns+4)  !nessesairy if one does not use the EE1 predictor to guess
                      !attention! T(-1) is not part of T (Ghostpoint)
  !---------------------------------------------------------------------
  !     compute angular velocity (theta_dot)
  !---------------------------------------------------------------------  
  ! theta_dot was removed from the NL system and is now computed  
  if (TimeMethodSolid == "EI1") then  
    C1=1.0 ! dt factor
    C2=0.d0 ! factor for RHS
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.d0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "CN2") then  
    C1=2.0 ! dt factor
    C2=1.0 ! rhs old factor
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.d0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "BDF2") then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R)   ! dt factor
    C2 = 0.d0      ! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R)   ! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R)   ! factor before the THETA_DOT_N-1 term 
  endif
  
  beam(:,6) = (C1/dt) * ( beam(:,5) - C3*beam_old(:,5) - C4*beam_oldold(:,5) ) - C2*beam_old(:,6)
  
  !-----------------------------------------------------------------------------
  ! NOTE: 2013, this is extended to take more than one beam into account.
  ! however, its too hard to reprogram everything, so now we force it to be compatible
  !-----------------------------------------------------------------------------
  beam_solid%x(0:ns-1) = beam(:,1)
  beam_solid%y(0:ns-1) = beam(:,2)
  beam_solid%vx(0:ns-1) = beam(:,3)
  beam_solid%vy(0:ns-1) = beam(:,4)
  beam_solid%theta(0:ns-1) = beam(:,5)
  beam_solid%theta_dot(0:ns-1) = beam(:,6)
  beam_solid%pressure_old(0:ns-1) = pressure_old 
  beam_solid%pressure_new(0:ns-1) = pressure_new 
  beam_solid%tau_old(0:ns-1) = tau_beam_old
  beam_solid%tau_new(0:ns-1) = tau_beam_new
  
  !-----------------------------------------------------------------------------
  !     get deflection line (integrate_position)
  !-----------------------------------------------------------------------------  
  ! this beam is the new one, so its time+dt ( to get mouvement at the right instant)
  call integrate_position ( time+dt, beam_solid )
  
  !-----------------------------------------------------------------------------
  !  accelerations
  !-----------------------------------------------------------------------------  
  ! we have the old velocity (t_n) and the new one 
  beam_solid%ax(0:ns-1) = (beam_solid%vx(0:ns-1) - beam_old(:,3)) / dt
  beam_solid%ay(0:ns-1) = (beam_solid%vy(0:ns-1) - beam_old(:,4)) / dt
  
  !-----------------------------------------------------------------------------
  !       emergency brake
  !-----------------------------------------------------------------------------
  if ((maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0 ).and.(root)) then
    write (*,'(A)') "!!! IBES-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0"
    write (*,'(A)') "    That indicates a possible instability."
    write (*,'("time=",es11.4)') time
  endif
  
  !-----------------------------------------------------------------------------
  !     save number of iterations
  !-----------------------------------------------------------------------------    
  if (root) then
    open (14, file = 'IBES_iter.t', status = 'unknown', position='append')
    write (14, '(es11.4,1x,i3)') time, iter
    close (14)
  endif
  
  !-----------------------------------------------------------------------------
  !     iterate ( skipped for CN2 and EI1 )
  !-----------------------------------------------------------------------------    
  if (ActuallyBDF2 .eqv. .true.) then   ! Remember to switch back to BDF2 
    TimeMethodSolid = "BDF2" 
    ActuallyBDF2 = .false.
  endif  

  iCalls = iCalls + 1      ! count calls (to perform rare self-tests)

  ! The old beam at the current step is the oldold (t_n-1) at the next step
  beam_solid%beam_oldold(0:ns-1,1:6) = beam_old 
  ! for BDF2 with variable dt, we need the old time step
  ! this should be okay also when restarting, as the first CN2 step will provide us with the dt_old  
  beam_solid%dt_old = dt   

end subroutine IBES_solver






subroutine Solve_LGS ( J, F, x)
  !--------------------------------------------
  ! solves the linear system J*x = F
  !--------------------------------------------
  implicit none
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4), intent(in) :: J
  real(kind=pr),dimension(1:2*ns+4), intent(out) :: x
  real(kind=pr),dimension(1:2*ns+4), intent(in) :: F
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4) :: J2
  real(kind=pr) :: t0
  integer :: error, ipiv(1:2*ns+4)  
  t0 = MPI_wtime()
  
  J2 = transpose(J)
  call dgetrf( 2*ns+4, 2*ns+4, J2 , 2*ns+4, ipiv, error )
  if (error .ne. 0) then
    write(*,*) "!!! Crutial: dgetrf error.", error
    call abort()
  endif
  
  x = F
  call dgetrs( 'N', 2*ns+4, 1, J2, 2*ns+4, ipiv, x, 2*ns+4, error )
  if (error .ne. 0) then 
    write(*,*) "!!! Crutial: dgetrs error.", error
    call abort()
  endif
  
  time_LAPACK = time_LAPACK + MPI_wtime() - t0
end subroutine


!-------------------------------------------------------------------------------

subroutine F_nonlinear (time, dt, dt_old,  F, theta_old, theta_dot_old, theta, T, p,&
                        old_rhs, theta_oldold, theta_dot_oldold, tau_beam, beam_solid  )
  implicit none  
  ! returns the RHS of the nonlinear eqn set F(x) for a given beam (the iterating one) and the (fixed) previous one
  real(kind=pr),dimension(0:ns-1),intent (in) :: theta_old, theta_dot_old
  real(kind=pr),dimension(0:ns-1),intent (in) :: old_rhs, theta_oldold, theta_dot_oldold
  real(kind=pr),dimension(-1:ns+1), intent (in) :: theta
  type (solid) :: beam_solid
  real(kind=pr),dimension(-1:ns-1), intent (in) :: T
  real(kind=pr),dimension(0:ns-1), intent (in) :: p, tau_beam
  real(kind=pr),dimension(0:ns-1) :: theta_dot_new, tau_s
  real(kind=pr),dimension(1:(2*ns+4)), intent (out) :: F  
  real(kind=pr),intent (in) :: time, dt, dt_old
  real(kind=pr) :: K1,K2, C1,C2, C3,C4,R
  real(kind=pr) :: alpha, alpha_t, alpha_tt 
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  integer :: i
  
  call Differentiate1D (tau_beam, tau_s, ns, ds, 1)
  
  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
  !---------new extended boundarys (constants)
  K1 = mue*(LeadingEdge(5)*cos(alpha)+LeadingEdge(6)*sin(alpha)+grav*sin(alpha)) - tau_beam(0)
  K2 = p(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))

  F = 0. !to check if an index is not set. is not the case.
  
  if (TimeMethodSolid == "EI1") then  
    C1=1.0 ! dt factor
    C2=0.d0 ! factor for RHS
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.d0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "CN2") then  
    C1=2.0 ! dt factor
    C2=1.0 ! rhs old factor
    C3=1.0 ! factor before the THETA_DOT_N term
    C4=0.d0 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "BDF2") then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R)   ! dt factor
    C2 = 0.d0      ! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R)   ! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R)   ! factor before the THETA_DOT_N-1 term 
  endif
  
  theta_dot_new = (C1/dt)*(theta(0:ns-1)-C3*theta_old(0:ns-1)-C4*theta_oldold)-C2*theta_dot_old
  
  ! first we set the 8 special eqns at the beginning 
  F(1) = theta(0) 
  F(2) = (T(1)-T(-1))/(2.0*ds)  +  (eta/(2.0*ds**3))*(theta(-1)-&
         2.0*theta(0)+theta(1))*(theta(1)-theta(-1))  -  K1
  F(3) = T(0)*(theta(1)-theta(-1))/(2.0*ds)   -   &
         (eta/(2.0*ds**3))*(-3.0*theta(-1)+10.d0*theta(0)-12.0*theta(1)+6.0*theta(2)-theta(3))  -  K2
  F(4) = ( (1.0/12.0)*theta(ns-3) - (2.0/3.0)*theta(ns-2) + (2.0/3.0)*theta(ns)&
          - (1.0/12.0)*theta(ns+1) )/ds
  F(5) = ( (-1.0/12.0)*theta(ns-3) + (4.0/3.0)*theta(ns-2) - 2.5*theta(ns-1) + &
           (4.0/3.0)*theta(ns) + (-1.0/12.0)*theta(ns+1) )/ds**2
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
end subroutine F_nonlinear

!-------------------------------------------------------------------------------

subroutine Jacobi(time, dt, dt_old, J, T, theta, theta_old, theta_dot_old,&
                  p, theta_oldold, theta_dot_oldold, tau_beam, beam_solid)
  implicit none
  real(kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real(kind=pr), dimension(-1:ns-1), intent (in) :: T
  real(kind=pr), intent (in) :: time, dt, dt_old
  type(solid) :: beam_solid
  real(kind=pr), dimension(-1:ns+1), intent (in) :: theta
  real(kind=pr), dimension(0:ns-1), intent (in) :: p, theta_old, theta_dot_old, theta_oldold, theta_dot_oldold, tau_beam
  real(kind=pr), dimension(0:ns-1) ::theta_dot_new
  integer :: i, T0_index, k, l, m
  real(kind=pr) :: alpha, alpha_t, alpha_tt, C1,C2,C3,C4,D,R
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  
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
  
  if (TimeMethodSolid == "EI1") then  
    C1=1.0                 ! dt factor
    C2=0.d0                 ! factor for RHS
    C3=1.0                 ! factor before the THETA_DOT_N term
    C4=0.d0                 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "CN2") then  
    C1=2.0                 ! dt factor
    C2=1.0                 ! rhs old factor
    C3=1.0                 ! factor before the THETA_DOT_N term
    C4=0.d0                 ! factor before the THETA_DOT_N-1 term
  elseif (TimeMethodSolid == "BDF2") then  
    R  = dt / dt_old
    C1 = (1.+2.*R)/(1.+R)   ! dt factor
    C2 = 0.d0      ! rhs old factor
    C3 = ((1.+R)**2)/(1.+2.*R)   ! factor before the THETA_DOT_N term
    C4 = (-R**2 )/(1.+2.*R)   ! factor before the THETA_DOT_N-1 term 
  endif
  D  = C1/dt ! this is the only entry in the jacobian for theta_dot_new
  theta_dot_new = (C1/dt) * ( theta(0:ns-1) - C3*theta_old(0:ns-1) - C4*theta_oldold ) - C2*theta_dot_old
  
  
  T0_index = ns+4               !T0 means here the first T = T(-1)
  J = 0.d0                ! initialize J
  !-- set first 6 eqns 
  J(2,1) = 1.0
  !----eqn 2 (special BC 1)
  J(1,2) = (eta/(2.0*ds**3))*(-2.0*theta(-1)+2.0*theta(0))    ! dF2 / dTheta(-1)
  J(2,2) =-(eta/(ds**3))*(theta(1)-theta(-1))        ! dF2 / dTheta(0)
  J(3,2) = (eta/(2.0*ds**3))*(-2.0*theta(0) +2.0*theta(1))    ! dF2 / dTheta(1)
  J(T0_index,2)   = -0.5/ds            ! dF2 / dT(-1)
  J(T0_index+2,2) =  0.5/ds            ! dF2 / dT(+1)
  !----eqn 3 (special BC 2)
  J(1,3) = (-0.5*eta/(ds**3))*(-3.0) - T(0)/(2.0*ds)      ! dF3 / dTheta(-1)
  J(2,3) = (-0.5*eta/(ds**3))*(10.d0)          ! dF3 / dTheta(0)
  J(3,3) = (-0.5*eta/(ds**3))*(-12.0)+ T(0)/(2.0*ds)      ! dF3 / dTheta(+1)
  J(4,3) = (-0.5*eta/(ds**3))*(6.0)          ! dF3 / dTheta(+2)
  J(5,3) = (-0.5*eta/(ds**3))*(-1.0)          ! dF3 / dTheta(+3)
  J(T0_index+1,3) = (theta(1)-theta(-1))/(2.0*ds)      ! dF3 / dT(0)
  
  
  ! the following to eqns have been changed; actually, you could multiply them
  ! by any factor, but its easier to read like this
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

  !-------------7th eqn: seperated 1
  J(1,7) = T(0)*(theta(1) -theta(-1))/(2.0*ds**2) - p(0)/(2.0*ds) &
           -(eta/(ds**4))*(-2.5*theta(0)+9.0*theta(1)-12.0*theta(2) +7.0*theta(3)-1.5*theta(4)) &
           +2.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4)
  J(2,7) = -2.5*eta*(theta(1)-theta(-1))/(ds**4) - 4.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4) &
           + 2.0*mue*(theta_dot_new(0)+alpha_t)*D
  J(3,7) = -T(0)*(theta(1) -theta(-1))/(2.0*ds**2) + p(0)/(2.0*ds) &
           +(eta/(ds**4))*(-2.5*theta(0)+9.0*theta(1)-12.0*theta(2) +7.0*theta(3)-1.5*theta(4)) &
           +9.0*eta*(theta(1)-theta(-1))/(ds**4) + 2.0*eta*(theta(1)-2.0*theta(0)+theta(-1))/(ds**4)
  J(4,7) = -12.0*eta*(theta(1)-theta(-1))/(ds**4)
  J(5,7) =   7.0*eta*(theta(1)-theta(-1))/(ds**4)
  J(6,7) = - 1.5*eta*(theta(1)-theta(-1))/(ds**4)
  J(ns+4,7) = 1.0/(ds**2)
  J(ns+5,7) = (-2.0-((theta(1)-theta(-1))**2)/4.0)/(ds**2)
  J(ns+6,7) = 1.0/(ds**2)
  
  !--------------8th eqn
  J(2*ns+4,8) = (dt/(C1*mue*ds**2))* ((theta(ns)-2.0*theta(ns-1)+theta(ns-2)) &
              + (3.0/2.0)*(theta(ns)-theta(ns-2)))
  J(2*ns+3,8) = -4.0*(dt/(C1*mue))*(theta(ns)-theta(ns-2))/(2.0*ds**2)
  J(2*ns+2,8) =  1.0*(dt/(C1*mue))*(theta(ns)-theta(ns-2))/(2.0*ds**2)
  J(ns+3,8)   = -dt*eta/((C1*mue)*ds**4)
  
  J(ns+2,8)   = (dt/(C1*mue*ds**2))   *   ( 4.0*eta/(ds**2) + (T(ns-1) &
              + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2) &
              + tau_beam(ns-1)*ds/(2.0) &   ! note we pulled 1/ds**2 out
              + eta*(theta(ns)-2.0*theta(ns-1)+theta(ns-2))*(theta(ns)-theta(ns-2))/(2.0*ds**2) &
              + (3.0*T(ns-1)-4.0*T(ns-2)+T(ns-3))/2.0         )
              
  J(ns  ,8)   = (dt/(C1*mue*ds**2))   *   ( 4.0*eta/(ds**2) + (T(ns-1) &
              + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2) &
              - tau_beam(ns-1)*ds/(2.0) &
              - eta*(theta(ns)-2.0*theta(ns-1)+theta(ns-2))*(theta(ns)-theta(ns-2))/(2.0*ds**2) &
              - (3.0*T(ns-1)-4.0*T(ns-2)+T(ns-3))/2.0         )

  J(ns+1,8)   = -D + (dt/(C1*mue))*( -6.0*eta/(ds**4) -2.0*(T(ns-1) &
                + eta*( (theta(ns)-theta(ns-2))/(2.0*ds) )**2)/(ds**2) -D*sigma)              
              
  J(ns-1,8)   = -dt*eta/(C1*mue*ds**4)
  
  
!checked irreg points 16/01/11
  !-- all irregular points are filled. typing sucked. 
  do i=1, ns-2
    k = 8 + i !index for eqns F (row) -- that means first row is row 9
    l = i + 2 !index column (theta) -- starts @ theta_1 = index 3
    m = i + ns + 5 !index column (T) -- starts @ T(-1) = ns+5
    J( l ,k)  = -D + (dt/(C1*mue)) * (-6.0*eta/(ds**4) &
              -(2.0/ds**2)*(T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 ) - D*sigma  )
    J( l-1,k) = (dt/(C1*mue))*( 4.0*eta/(ds**4) + (T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 )/(ds**2)  & 
              - tau_beam(i)/(2.0*ds) &
              - eta*(theta(i+1)-2.0*theta(i)+theta(i-1))*(theta(i+1)-theta(i-1))/(2.0*ds**4) &
              - (T(i+1)-T(i-1))/(2.0*ds**2)  ) 
    J( l+1,k) = (dt/(C1*mue))*( 4.0*eta/(ds**4) + (T(i)+eta*((theta(i+1)-theta(i-1))/(2.0*ds))**2 )/(ds**2)  & 
              + tau_beam(i)/(2.0*ds) &
              + eta*(theta(i+1)-2.0*theta(i)+theta(i-1))*(theta(i+1)-theta(i-1))/(2.0*ds**4) &
              + (T(i+1)-T(i-1))/(2.0*ds**2)  ) 
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
                 - eta* (theta(i+2)-2.0*theta(i+1)+2.0*theta(i-1)-theta(i-2))/(2.0*ds**4) &
                 + 2.0*eta*(theta(i+1)-2.0*theta(i) +theta(i-1))/(ds**4) 
    J( l+1,k) = -T(i)*(theta(i+1)-theta(i-1))/(2.0*ds**2) +p(i)/(2.0*ds) -eta*(theta(i+1)-theta(i-1))/(ds**4) & 
                 + eta* (theta(i+2)-2.0*theta(i+1)+2.0*theta(i-1)-theta(i-2))/(2.0*ds**4) &
                 + 2.0*eta*(theta(i+1)-2.0*theta(i) +theta(i-1))/(ds**4) 
    J( l-2,k) = -eta*(theta(i+1)-theta(i-1))/(2.0*ds**4) 
    J( l+2,k) =  eta*(theta(i+1)-theta(i-1))/(2.0*ds**4) 
    !block dG/dT
    J( m,k)   = -2.0/(ds**2) - ( (theta(i+1)-theta(i-1) )/(2.0*ds) )**2 
    J( m-1,k) = 1.0/(ds**2) 
    J( m+1,k) = 1.0/(ds**2) 
  enddo
end subroutine Jacobi


!-------------------------------------------------------------------------------


subroutine GravityImpulse(time)
  implicit none
  ! gives a little gravity impulse to pertubate the beam between T0 and T1 (sinusoidal)
  real(kind=pr), intent (in) ::  time
  real(kind=pr) :: T0,T1,a,b,c,d,k,t
  
end subroutine GravityImpulse

! ------------------------------------------------------------------------------




subroutine RHS_beameqn (time, theta, theta_dot, pressure_beam, T, tau_beam, beam_solid)
    !---------------------------------------------------------------------------
    ! Beam Equation right hand side at time_n
    ! Version 19.09.2012, completely debugged, gives exactly the same results as the matlab solver.
    ! INPUT
    ! time:    the time at which we compute the RHS. Note this is important: The leading edge motion 
    !      may be time dependent, (also in between a runge kutta step) so we have to call the subroutine
    !      at the right time
    ! pressure_beam  the pressure jump at time time
    ! OUTPUT:
    ! T      the tension in the beam at time time. we return it to form the complete initial guess for the implicit solvers
    !
    !---------------------------------------------------------------------------
    implicit none
    real(kind=pr), intent (in) ::         time
    real(kind=pr) ::              A1, A2, K2,C2
    real(kind=pr),dimension(0:ns-1), intent (in) :: pressure_beam, tau_beam
    real(kind=pr),dimension(0:ns-1), intent (out) :: T
    type(solid) :: beam_solid
    real(kind=pr),dimension(0:ns-1), intent (inout) ::   theta, theta_dot
    real(kind=pr),dimension(0:ns+2) :: theta_extended, theta_extended_s, theta_extended_ss
    real(kind=pr),dimension(0:ns+2) :: theta_extended_sss, theta_extended_ssss
    real(kind=pr),dimension(0:ns-1) :: theta_s, theta_ss, theta_sss, theta_ssss, T_s, p_s
    real(kind=pr) :: alpha, alpha_t, alpha_tt
    real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  

    call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )
    
    theta(0)     = 0.d0 ! first boundary condition, angle here is predescribed by the motion protocol.
    theta_dot(0) = 0.d0 ! modified: set to zero. was alpha, but is now angle in RELATIVE system (23.02.2011)
    
    !-- compute the tension in the beam
    call Tension (time, T, T_s, theta, theta_dot, pressure_beam, tau_beam, beam_solid)   
    
    ! --------------------------------------------------------------------------
    ! Extend the beam with ghostpoints to fulfill the boundary conditions
    ! --------------------------------------------------------------------------
    !-- leading edge boundary conditions (constants)
    K2 = pressure_beam(0) + mue*(LeadingEdge(6)*cos(alpha)-LeadingEdge(5)*sin(alpha)+grav*cos(alpha))
    C2 = 2.0*ds*K2 - T(0)*theta(1) + (eta/ds**2)*(10.d0*theta(0)-12.0*theta(1)+6.0*theta(2)-theta(3) )

    !-- theta_extended(0) is the first virtual node. 
    theta_extended(0)    = C2 / ( (3.0*eta/ds**2)-T(0) )
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

    ! --------------------------------------------------------------------------
    ! Last step: compute evolution equation for theta
    ! --------------------------------------------------------------------------
    theta     = theta_dot
    theta_dot = (-p_s - eta*theta_ssss + theta_ss*(T+eta*(theta_s**2)) &
                 +2.*T_s*theta_s - mue*alpha_tt  &
                 -sigma*theta_dot + tau_beam*theta_s) / mue
    theta_dot(0) = 0.d0  
 
end subroutine RHS_beameqn


! ------------------------------------------------------------------------------


subroutine Tension ( time, T, T_s, theta, theta_dot, pressure, tau_beam, beam_solid)
  ! note we use a TRIAG solver. the neumann condition is not in the matrix, it is set explicitly. therefore we return also T_s 
  implicit none
  real(kind=pr), intent (in) :: time
  real(kind=pr),dimension(0:ns-1) :: theta_s, theta_ss, theta_sss
  real(kind=pr),dimension(0:ns-1) :: tau_beam_s
  real(kind=pr),dimension(0:ns-1), intent (in) ::  theta, theta_dot, pressure, tau_beam
  real(kind=pr),dimension(0:ns-1), intent (out) :: T, T_s
  type(solid) :: beam_solid
  real(kind=pr):: K1, T_s0
  real(kind=pr):: alpha, alpha_t, alpha_tt
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  real(kind=pr),dimension(0:ns-1) :: diagonal_main
  real(kind=pr),dimension(0:ns-2) :: diagonal_up, diagonal_down
  real(kind=pr),dimension(0:ns-3) :: diagonal_up2
  real(kind=pr),dimension(0:ns-1) :: ipiv !used only for the MKL lib
  real(kind=pr),dimension(0:ns-1) :: rhs
  integer :: i, info
  
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
  diagonal_main   = -2.0/(ds**2) - theta_s**2
  diagonal_up     =  1.0/(ds**2)
  diagonal_down   =  1.0/(ds**2)  
  ! Dirichlet Condition
  diagonal_main(ns-1)  = 1.0 ! note this point is overwritten
  diagonal_down(ns-2)   = 0.d0 ! for the dirichlet condition   
  ! Neumann Condition
  diagonal_up(0)   = 2.0/(ds**2)
  
  ! construct RHS for T-equation
  rhs(0:ns-1) = -(pressure*theta_s) -(2.*eta*theta_s*theta_sss)-(eta*theta_ss**2 ) &
              -mue*(theta_dot + alpha_t)**2  - tau_beam_s
  ! Dirichlet
  rhs(ns-1)     = 0.d0
  ! Neumann
  rhs(0)     = rhs(0) + T_s0 * 2.0/ds

  ! preconditioner - improve condition number of system by normalizing the diagonal
  do i=0,ns-1
    rhs(i) = rhs(i) / diagonal_main(i)
  enddo
  
  do i=0,ns-2
    diagonal_up(i)   = diagonal_up(i) / diagonal_main(i)
    diagonal_down(i) = diagonal_down(i) / diagonal_main(i+1)
  enddo  
  diagonal_main   = 1.0


  ! solve system LU decomposition of matrix
  call dgttrf( ns, diagonal_down, diagonal_main, diagonal_up, diagonal_up2, ipiv, info )
  call dgttrs( 'N', ns, 1, diagonal_down, diagonal_main, diagonal_up, diagonal_up2,&
               ipiv, rhs, ns, info )  ! solve. now g is the solution  

  
  T = rhs  ! return tension
  
  call Differentiate1D ( T, T_s, ns, ds, 1)
  
  T_s (0) = T_s0 ! don't forget the Neumann Condition
 
  if (info.ne.0) then 
  write (*,*) '!!! MKL linear solver was experiencing trouble.', info
  call abort()
  endif
end subroutine Tension

! ------------------------------------------------------------------------------

subroutine Differentiate1D (f, f_derivative, N, dx, order)
  implicit none
  integer, intent (in)         :: order, n
  real(kind=pr), intent (in)            :: dx
  real(kind=pr), dimension(0:N-1)      :: f, f_derivative
  real(kind=pr),dimension(0:N-1, 0:N-1)   :: D1, D2, D3 ,D4
  
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
    case default
      if (root) write(*,*) "Differentiate1D wrong choice"
      call abort()
  end select
  
end subroutine Differentiate1D
  
! ------------------------------------------------------------------------------

subroutine create_diff_matrices (D1, D2, D3, D4, nn, dx)
  implicit none
  integer :: i,N
  integer, intent (in) :: nn
  real(kind=pr), intent(in) :: dx
  real(kind=pr),dimension(0:nn-1, 0:nn-1), intent (out) :: D1, D2, D3 ,D4

  !------------------------------------
  ! create diff-matrices.
  ! first derivative: x'=D1 * x
  !------------------------------------

  D1=0.d0
  D2=0.d0
  D3=0.d0
  D4=0.d0

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


!-------------------------------------------------------------------------------
!   wrapper for new beam dataype for RK4
! This routine is called when really using RK4 for the beam (as opposed to use
! it just for predicting the next beam)
!-------------------------------------------------------------------------------
subroutine RK4_wrapper ( time, dt, beam_solid )! note this is actuall only ONE beam!!
  implicit none
  real(kind=pr), intent (in) :: dt, time
  type (solid), intent(inout) :: beam_solid
  real(kind=pr),dimension(0:ns-1) :: T
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam

  beam(0:ns-1,1) = beam_solid%x(0:ns-1)
  beam(0:ns-1,2) = beam_solid%y(0:ns-1)
  beam(0:ns-1,3) = beam_solid%vx(0:ns-1)
  beam(0:ns-1,4) = beam_solid%vy(0:ns-1)
  beam(0:ns-1,5) = beam_solid%theta(0:ns-1)
  beam(0:ns-1,6) = beam_solid%theta_dot(0:ns-1)
  
  call RK4s( time, dt, beam, beam_solid%pressure_old(0:ns-1), T, beam_solid%tau_old(0:ns-1), &
            beam_solid )
            
  beam_solid%x(0:ns-1)= beam(0:ns-1,1)
  beam_solid%y(0:ns-1) = beam(0:ns-1,2)
  beam_solid%vx(0:ns-1) = beam(0:ns-1,3)
  beam_solid%vy(0:ns-1) = beam(0:ns-1,4)
  beam_solid%theta(0:ns-1) = beam(0:ns-1,5)
  beam_solid%theta_dot(0:ns-1) = beam(0:ns-1,6)            
            
  call integrate_position( time, beam_solid )                
end subroutine RK4_wrapper




!-------------------------------------------------------------------------------
subroutine RK4s (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid) 
  implicit none
  real(kind=pr), intent (in) :: dt_beam, time
  real(kind=pr),dimension(0:ns-1, 1:6), intent (inout) :: beam
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam_old
  real(kind=pr),dimension(0:ns-1), intent (in) :: pressure_beam, tau_beam
  type (solid) :: beam_solid
  real(kind=pr),dimension(0:ns-1), intent (inout) :: T
  real(kind=pr),dimension(0:ns-1) :: theta, theta_dot, T1, T2, T3, T4
  real(kind=pr),dimension(0:ns-1) :: theta_1, theta_2, theta_3, theta_4
  real(kind=pr),dimension(0:ns-1) :: theta_dot_1, theta_dot_2, theta_dot_3, theta_dot_4
  beam_old=beam
  theta = beam(:,5)
  theta_dot = beam(:,6)

  theta_1 = 0.d0
  theta_2 = 0.d0
  theta_3 = 0.d0
  theta_4 = 0.d0
  theta_dot_1 = 0.d0
  theta_dot_2 = 0.d0
  theta_dot_3 = 0.d0
  theta_dot_4 = 0.d0
  
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
  if (maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0 ) then
    write (*,'(A)') "!!! rk4-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0"
    write (*,'(A)') "possible instability"
    write (*,'("time=",es11.4, " dt=",es11.4)') time, dt_beam
  endif  
  
end subroutine RK4s


!-------------------------------------------------------------------------------
!   wrapper for new beam dataype for RK4
! This routine is called when really using RK4 for the beam (as opposed to use
! it just for predicting the next beam)
!-------------------------------------------------------------------------------
subroutine EE1_wrapper ( time, dt, beam_solid )! note this is actuall only ONE beam!!
  implicit none
  real(kind=pr), intent (in) :: dt, time
  type (solid), intent(inout) :: beam_solid
  real(kind=pr),dimension(0:ns-1) :: T
  real(kind=pr),dimension(0:ns-1, 1:6) :: beam

  beam(:,1) = beam_solid%x(0:ns-1)
  beam(:,2) = beam_solid%y(0:ns-1)
  beam(:,3) = beam_solid%vx(0:ns-1)
  beam(:,4) = beam_solid%vy(0:ns-1)
  beam(:,5) = beam_solid%theta(0:ns-1)
  beam(:,6) = beam_solid%theta_dot(0:ns-1)
  
  call EE1s( time, dt, beam, beam_solid%pressure_old(0:ns-1), T, &
       beam_solid%tau_old(0:ns-1), beam_solid )

  beam_solid%x(0:ns-1) = beam(:,1)
  beam_solid%y(0:ns-1) = beam(:,2)
  beam_solid%vx(0:ns-1) = beam(:,3)
  beam_solid%vy(0:ns-1) = beam(:,4)
  beam_solid%theta(0:ns-1) = beam(:,5)
  beam_solid%theta_dot(0:ns-1) = beam(:,6)
  
  call integrate_position( time, beam_solid )         
end subroutine EE1_wrapper


!-------------------------------------------------------------------------------
  
subroutine EE1s (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid) 
  implicit none
  real(kind=pr), intent (in) :: dt_beam, time
  real(kind=pr),dimension(0:ns-1, 1:6), intent (inout) :: beam
  real(kind=pr),dimension(0:ns-1), intent(out) :: T
  real(kind=pr),dimension(0:ns-1), intent (in) :: pressure_beam, tau_beam
  type (solid) :: beam_solid
  real(kind=pr),dimension(0:ns-1) :: theta, theta_dot
  real(kind=pr),dimension(0:ns-1) :: theta_1
  real(kind=pr),dimension(0:ns-1) :: theta_dot_1

  theta     = beam(:,5)
  theta_dot = beam(:,6)

  theta_1     = theta
  theta_dot_1 = theta_dot
  
  call RHS_beameqn (time, theta_1 , theta_dot_1, pressure_beam, T, tau_beam, beam_solid)

  beam(:,5) = theta + dt_beam * theta_1
  beam(:,6) = theta_dot + dt_beam * theta_dot_1
 
end subroutine EE1s


!-------------------------------------------------------------------------------


subroutine Jacobi_num(time, dt,dt_old, J, T, theta, theta_old, &
       theta_dot_old, p, theta_oldold, theta_dot_oldold, old_rhs, tau_beam_new, beam_solid)
  implicit none
  type (solid) :: beam_solid
  real(kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real(kind=pr), dimension(1:2*ns+4) :: F1,F2
  real(kind=pr), dimension(-1:ns-1), intent (in) :: T
  real(kind=pr), dimension(-1:ns-1) :: T1, T2
  real(kind=pr), intent (in) :: time, dt,dt_old
  real(kind=pr), dimension(-1:ns+1), intent (in) :: theta
  real(kind=pr), dimension(-1:ns+1) :: theta1, theta2
  real(kind=pr), dimension(0:ns-1), intent (in) :: p, theta_old, theta_dot_old, &
  theta_oldold,theta_dot_oldold, old_rhs, tau_beam_new
  integer :: i, l

  !---------------------------------------
  !  NUMERIC COMPUTATION OF THE JACOBIAN
  !  -> I think the row/column ordering in J (in IBES) is ackward; pay attention.
  !---------------------------------------
 
  ! derive wrt to theta_new
  l=1
  do i=-1,ns+1
      theta1    = theta
      theta1(i) = theta1(i)-1.0e-6
      theta2    = theta
      theta2(i) = theta2(i)+1.0e-6
      
      call F_nonlinear (time, dt, dt_old, F1, theta_old, theta_dot_old, theta1, &
      T, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid )
      call F_nonlinear (time, dt, dt_old, F2, theta_old, theta_dot_old, theta2, &
      T, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid )
      
      J(l,:)= (F2-F1)/(2.0e-6)
      l=l+1

  enddo

  ! derive wrt to T_new
  do i=-1,ns-1
      T1=T
      T1(i) = T1(i)-1.0e-6
      T2=T
      T2(i) = T2(i)+1.0e-6
      
      call F_nonlinear (time, dt, dt_old, F1, theta_old, theta_dot_old, theta, &
      T1, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid  )
      call F_nonlinear (time, dt, dt_old, F2, theta_old, theta_dot_old, theta, &
      T2, p, old_rhs, theta_oldold, theta_dot_oldold, tau_beam_new, beam_solid  )
      
      J(l,:)= (F2-F1)/(2.0e-6)
      l=l+1
  enddo
end subroutine Jacobi_num


!-------------------------------------------------------------------------------


logical function Vector_isNAN(f)
  real(kind=pr), intent(in) :: f(:)
  integer :: a, i
  a = size(f)
  Vector_isNAN = .false.
  do i=1,a
    if (.not.(f(i).eq.f(i))) then
      Vector_isNAN = .true.
      return
    endif
  enddo
end function


!-------------------------------------------------------------------------------


subroutine show_solid_model_information
  implicit none
  write(*,*) "*****************************************"
  write(*,*) "solid solver debuging: showing all globals"
  write(*,*) "*****************************************"
  write(*,'("mpirank=",i5)') mpirank
  write(*,'("ns=",i3)') ns
  write(*,'("mue=",es12.4," eta=",es12.4," grav=",es12.4)') mue,eta,grav
  write(*,'("sigma=",es12.4," t_beam=",es12.4," ds=",es12.4)') sigma, t_beam, ds
  write(*,'("T_release=",es12.4," TimeMethodSolid=",A," tau=",es12.4)') &
  T_release, trim(TimeMethodSolid), tau
  write(*,*) "*******************end*******************"
end subroutine


!-------------------------------------------------------------------------------
 
subroutine show_beam( beam )
  use vars
  implicit none
  character(len=20):: str
  type(solid),intent(in) :: beam
  if (root) then
    write(str,'("(A,1x,",i3.3,"(f7.3,1x))")') ns
    write(*,str) "beam%x ", beam%x(0:ns-1)
    write(*,str) "beam%y ", beam%y(0:ns-1)
    write(*,str) "beam%vx", beam%vx(0:ns-1)
    write(*,str) "beam%vy", beam%vy(0:ns-1)
    write(*,*) "x0", beam%x0
    write(*,*) "y0", beam%y0
    write(*,*) "dt_old", beam%dt_old
    write(*,str) "beam%theta    ", beam%theta(0:ns-1)
    write(*,str) "beam%theta_dot", beam%theta_dot(0:ns-1)  
    write(*,str) "beam%pressure_new", beam%pressure_new(0:ns-1)
    write(*,str) "beam%pressure_old", beam%pressure_old(0:ns-1)
    write(*,str) "beam%tau_new", beam%tau_new(0:ns-1)
    write(*,str) "beam%tau_old", beam%tau_old(0:ns-1)
    write(*,str) "beam%beam_oldold(:,1)", beam%beam_oldold(0:ns-1,1)
    write(*,str) "beam%beam_oldold(:,2)", beam%beam_oldold(0:ns-1,2)
    write(*,str) "beam%beam_oldold(:,3)", beam%beam_oldold(0:ns-1,3)
    write(*,str) "beam%beam_oldold(:,4)", beam%beam_oldold(0:ns-1,4)
    write(*,str) "beam%beam_oldold(:,5)", beam%beam_oldold(0:ns-1,5)
    write(*,str) "beam%beam_oldold(:,6)", beam%beam_oldold(0:ns-1,6)
    write(*,*) "-->end beam<--"
  endif
end subroutine show_beam


!-------------------------------------------------------------------------------


subroutine show_beam_on_error( beam )
  use vars
  implicit none
  type(solid),intent(in) :: beam
  if ( Vector_isNAN(beam%x).or.&
       Vector_isNAN(beam%y).or.&
       Vector_isNAN(beam%vx).or.&
       Vector_isNAN(beam%vy).or.&
       Vector_isNAN(beam%theta).or.&
       Vector_isNAN(beam%theta_dot).or.&
       Vector_isNAN(beam%pressure_new).or.&
       Vector_isNAN(beam%pressure_old).or.&
       Vector_isNAN(beam%tau_new).or.&
       Vector_isNAN(beam%tau_old).or.&
       Vector_isNAN(beam%beam_oldold(:,1)).or.&
       Vector_isNAN(beam%beam_oldold(:,2)).or.&
       Vector_isNAN(beam%beam_oldold(:,3)).or.&
       Vector_isNAN(beam%beam_oldold(:,4)).or.&
       Vector_isNAN(beam%beam_oldold(:,5)).or.&
       Vector_isNAN(beam%beam_oldold(:,6)) &
     ) then
     
     call show_beam( beam )
     call abort()  
  endif
end subroutine show_beam_on_error


!-------------------------------------------------------------------------------


subroutine lapack_unit_test()
  !-----------------------------------------------------------------------------
  ! Minimal example for testing the Lapack library
  ! Solves the linear system
  ! I*x = b
  ! where I is the identity matrix, so the solution is trivially x=b
  !-----------------------------------------------------------------------------
  use vars
  implicit none
  integer, parameter:: n = 20
  real(kind=pr) :: a1(1:n,1:n),a2(1:n,1:n), b(1:n), x(1:n), err
  integer :: i,error
  integer :: ipiv(1:n)
  a2 = 0.d0
  a1 = 0.d0
  b  = 7.d0
  x  = 0.d0

  if (root) write(*,*) "--------------------------------"
  if (root) write(*,*) " Starting LAPACK unit test"
  
  ! create identity matrix
  do i=1,n
   a1(i,i) = 1.d0
  enddo

  !-----------------------------------------------------------------------------
  ! first step: factorization of the matrix a1 (which is overwritten on output
  ! so we make a copy of it)
  !-----------------------------------------------------------------------------  
  a2=a1 ! lapack overwrites a2 with the factorization
  call dgetrf ( n, n, a2, n, ipiv, error )
  
  !-----------------------------------------------------------------------------
  ! second step: backwards substitution
  !-----------------------------------------------------------------------------    
  x = b ! lapack overwrites x with the solution
  call dgetrs( 'N', n, 1, a2, n, ipiv, x, n, error ) 
  
  err = maxval( x - b )
  if ( err>1.d-12 ) then
    write (*,*) "LAPACK unit test failed"
    call abort()
  endif
  
  if (root) write(*,'(" Done. err=",es15.8)') err
  if (root) write(*,*) "--------------------------------"
end subroutine lapack_unit_test



!-------------------------------------------------------------------------
! dump runtime backup for the solid solver
!-------------------------------------------------------------------------
subroutine dump_solid_backup( time, beams, nbackup )
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  type(solid), dimension(1:nBeams), intent(in) :: beams
  integer, intent (in):: nbackup
  integer :: i
  character(len=24) :: filename
  
  !-- only root rank dumps backup
  if (root) then
    write(*,'(A)',advance='no') "Backuping solid solver..."
    write(filename,'("runtime_backup",i1,".fsi_bckp")') nbackup
    write(*,'(A)',advance='no') "file="//filename
    write(*,'(" time=",e11.4)',advance='no') time
    
    
    open(14,file=filename,status='replace',form='formatted')
    
    
    write(14,*) time
    write(14,*) ns, nBeams
    
    do i=1,nBeams      
      write(14,*) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy
      write(14,*) beams(i)%theta, beams(i)%theta_dot, beams(i)%ax, beams(i)%ay
      write(14,*) beams(i)%pressure_old, beams(i)%pressure_new
      write(14,*) beams(i)%tau_old, beams(i)%tau_new
      write(14,*) beams(i)%beam_oldold
      write(14,*) beams(i)%Force, beams(i)%Force_unst, beams(i)%Force_press
      write(14,*) beams(i)%Inertial_Force, beams(i)%x0, beams(i)%y0, beams(i)%AngleBeam
      write(14,*) beams(i)%phase, beams(i)%dt_old, beams(i)%drag_unst_new
      write(14,*) beams(i)%drag_unst_old, beams(i)%lift_unst_new, beams(i)%lift_unst_old
      write(14,*) beams(i)%StartupStep, beams(i)%UnsteadyCorrectionsReady
    enddo
    
    close(14)
    write(*,'(A)',advance='yes') "...DONE!"
  endif
end subroutine dump_solid_backup



!-------------------------------------------------------------------------
! read runtime backup for the solid solver
!-------------------------------------------------------------------------
subroutine read_solid_backup( beams, filename )
  use vars
  implicit none
  
  type(solid), dimension(1:nBeams), intent(inout) :: beams
  integer :: ns_file, nBeams_file, i
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time
  
  !-- all ranks read from file 
  if (root) write(*,'(A)',advance='no') "Reading in backup of solid solver: "//filename
  
  open(14,file=filename,status='old',form='formatted',action='read')
  read(14,*) time
  read(14,*) ns_file, nBeams_file

  if (root) write(*,'("(time=",e11.4,")")',advance='no') time
  
  if ((ns_file/=ns).or.(nBeams_file/=nBeams)) then
    write(*,*) "Cant retake solid backup: resolution ns or beam number doesnt match"
    call abort()
  endif
  
  do i =1, nBeams  
    read(14,*) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy
    read(14,*) beams(i)%theta, beams(i)%theta_dot, beams(i)%ax, beams(i)%ay
    read(14,*) beams(i)%pressure_old, beams(i)%pressure_new
    read(14,*) beams(i)%tau_old, beams(i)%tau_new
    read(14,*) beams(i)%beam_oldold
    read(14,*) beams(i)%Force, beams(i)%Force_unst, beams(i)%Force_press
    read(14,*) beams(i)%Inertial_Force, beams(i)%x0, beams(i)%y0, beams(i)%AngleBeam
    read(14,*) beams(i)%phase, beams(i)%dt_old, beams(i)%drag_unst_new
    read(14,*) beams(i)%drag_unst_old, beams(i)%lift_unst_new, beams(i)%lift_unst_old
    read(14,*) beams(i)%StartupStep, beams(i)%UnsteadyCorrectionsReady
  enddo
  
  close(14)

  if (root) write(*,'(A)',advance='yes') "...DONE!"
end subroutine read_solid_backup



!-------------------------------------------------------------------------
! read runtime backup for the solid solver (BINARY FILE)
!-------------------------------------------------------------------------
subroutine read_solid_backup_binary( beams, filename )
  use vars
  implicit none
  
  type(solid), dimension(1:nBeams), intent(inout) :: beams
  integer :: ns_file, nBeams_file, i
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time
  
  !-- all ranks read from file 
  if (root) write(*,'(A)',advance='no') "Reading in backup of solid solver: "//filename
  
  open(14,file=filename,status='old',form='unformatted',action='read')
  read(14) time
  read(14) ns_file, nBeams_file

  if (root) write(*,'("(time=",e11.4,")")',advance='no') time
  
  if ((ns_file/=ns).or.(nBeams_file/=nBeams)) then
    write(*,*) "Cant retake solid backup: resolution ns or beam number doesnt match"
    call abort()
  endif
  
  do i =1, nBeams  
    read(14) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy
    read(14) beams(i)%theta, beams(i)%theta_dot, beams(i)%ax, beams(i)%ay
    read(14) beams(i)%pressure_old, beams(i)%pressure_new
    read(14) beams(i)%tau_old, beams(i)%tau_new
    read(14) beams(i)%beam_oldold
    read(14) beams(i)%Force, beams(i)%Force_unst, beams(i)%Force_press
    read(14) beams(i)%Inertial_Force, beams(i)%x0, beams(i)%y0, beams(i)%AngleBeam
    read(14) beams(i)%phase, beams(i)%dt_old, beams(i)%drag_unst_new
    read(14) beams(i)%drag_unst_old, beams(i)%lift_unst_new, beams(i)%lift_unst_old
    read(14) beams(i)%StartupStep, beams(i)%UnsteadyCorrectionsReady
  enddo
  
  close(14)

  if (root) write(*,'(A)',advance='yes') "...DONE!"
end subroutine read_solid_backup_binary



subroutine convert_solid_bckp_ascii
  use vars
  implicit none
  type(solid), dimension(1:nBeams) :: beams
  
  inicond="nothing"
  call init_beams( beams )
  call read_solid_backup_binary( beams, "runtime_backup0.fsi_bckp" )
  call dump_solid_backup( 0.d0, beams, 0 )
  
end subroutine


  
end module solid_model


