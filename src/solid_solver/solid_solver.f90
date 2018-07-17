
module solid_model
  use vars ! for the precision statement
  ! we need the following line for presribed beams:
  use module_helpers
  use basic_operators
  implicit none

  !----------------------------------------------
  ! module global variables
  !----------------------------------------------
  integer,parameter :: nBeams = 1
  integer,save :: ns, debug_pressure
  ! see "type solid" about nsmax 06 Aug 2014
  integer,parameter :: nsmax = 300
  ! TODO: move these into the solid model datastructure
  real(kind=pr),save :: mue0
  real(kind=pr),save :: eta0
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

    ! coordinates, velocities and acceleration:
    real(kind=pr),dimension(0:nsmax) :: x,y,vx,vy,ax,ay
    ! note both theta and T include their ghosts [actually: theta(-1:ns+1), T(-1:ns-1)]
    real(kind=pr),dimension(-1:nsmax) :: theta, theta_dot, T
    ! external forces:
    real(kind=pr),dimension(0:nsmax) :: pressure_old, pressure_new
    real(kind=pr),dimension(0:nsmax) :: tau_old, tau_new
    ! material properties and their derivatives:
    real(kind=pr),dimension(0:nsmax) :: mu, mu_s, mu_star
    real(kind=pr),dimension(0:nsmax) :: zeta, zeta_s, zeta_ss, zeta_sss
    ! grid and width in rigid direction:
    real(kind=pr),dimension(0:nsmax) :: s, L_rigid
    ! values of theta and theta_dot at times t^(n) and t^(n-1)
    real(kind=pr),dimension(0:nsmax) :: theta_old, theta_oldold
    real(kind=pr),dimension(0:nsmax) :: theta_dot_old, theta_dot_oldold

    ! real(kind=pr),dimension(0:nsmax,1:6) :: beam_oldold
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
 include "aux.f90"
 include "prescribed_beam.f90"
 include "solid_solver_wrapper.f90"



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

  beam%E_kinetic = 0.5d0*ds*sum( beam%mu(0:ns-1) * (beam%vx(0:ns-1)**2 + beam%vy(0:ns-1)**2) )
  beam%E_pot     = grav*ds *sum( beam%mu(0:ns-1) * (beam%y(0:ns-1)-beam%y0) )
  beam%E_elastic = 0.5d0*ds*sum( beam%zeta(0:ns-1) * theta_s(0:ns-1)**2 )

  beam%Inertial_Force(1) = ds* sum( beam%mu(0:ns-1) * beam%ax(0:ns-1) )
  beam%Inertial_Force(2) = ds* sum( beam%mu(0:ns-1) * beam%ay(0:ns-1) )

end subroutine SolidEnergies


!-------------------------------------------------------------------------------
!   SOLID SOLVER ROUTINES
!-------------------------------------------------------------------------------
subroutine IBES_solver ( time, dt, beam )! note this is actuall only ONE beam!!
  implicit none
  real(kind=pr), intent (in) :: dt, time
  type (solid), intent(inout) :: beam
  real(kind=pr),dimension(0:ns-1) :: old_rhs
  real(kind=pr),dimension(1:2*ns+4) :: x, x_delta, F
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4) :: J
  real(kind=pr) :: err, C2,C1,C3,C4
  real(kind=pr) :: err_rel
  integer :: n,iter
  type(solid) :: beam_old
  logical :: ActuallyBDF2=.false., iterate=.true.

  beam_old = beam
  ! the current theta is at time t^n ( the time we call this routine )
  ! so afterwards it is going to be our old theta
  beam%theta_old(0:ns-1) = beam%theta(0:ns-1)
  beam%theta_dot_old(0:ns-1) = beam%theta_dot(0:ns-1)

  !-----------------------------------------------------------------------------
  ! Startup time step: use EE1 as first step of BDF2
  ! (pay attention to call InitializeSolidSolver() before the run)
  !-----------------------------------------------------------------------------
  if (beam%StartupStep) then  ! this is the first time step
    beam%StartupStep = .false.  ! we're about to do the first step
    if (TimeMethodSolid=="BDF2") then  ! if we deal with BDF2
      ActuallyBDF2 = .true.    ! Remember to switch back to BDF2 (at the end of the step)
      TimeMethodSolid = "EI1"  ! use "CN2" for the first step
    endif
  endif


  !-----------------------------------------------------------------------------
  ! Initial guess for the beam at the new time level.
  ! An euler-explicit step is made, ghostpoints are added
  !-----------------------------------------------------------------------------
  call initial_guess( time, dt, beam, x )


  !-----------------------------------------------------------------------------
  ! Newton iteration
  !-----------------------------------------------------------------------------
  ! this iterative process finds the solution "x" at the new time level
  ! note: during the iteration, only x is altered. the "beam" remains constant!
  !
  ! The iteration now uses both relative and absolute error. If x is large,
  ! we can have trouble reaching a very small increment, and iterate forever.
  ! If x is small, then the relative criterion tries to go much below
  ! machine precision. So we use both, either with the same precision.
  iterate = .true.
  err     = 1.0d0
  err_rel = 1.0d0
  iter    = 0

  do while (iterate)
    !  Calculate RHS vector
    call F_nonlinear( time, dt, F, x, beam )
    F = -F !newton raphson is J*dx = -F

    !  Create Jacobi Matrix
    call Jacobi( time, dt, J, x, beam )
    J = transpose(J)

    ! solve linear system
    call solve_linear_system ( J, F, x_delta )

    ! iterate
    iter    = iter + 1
    x       = x + x_delta
    err     = dsqrt(sum(x_delta**2))
    err_rel = abs(dsqrt(sum(x_delta**2)) / dsqrt(sum(x**2)))

    ! convergence test
    if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
      iterate = .false.
    endif

    ! emergency brake
    if (iter>1000) then
      call abort(333214,"!!! ERROR: IBES performed like 1000 iterations. this is not normal")
    endif
  enddo


  !--------------------------------------------------------------------
  ! update solid struct with new solution
  !--------------------------------------------------------------------
  ! we now have the solution vector x. note beam is completely untouched
  beam%theta(-1:ns+1) = x(1:ns+3) ! note we keep ghost nodes for fun
  beam%T(-1:ns-1)     = x(ns+4:2*ns+4)

  call time_marching_coefs(dt,beam%dt_old,C1,C2,C3,C4)
  beam%theta_dot(0:ns-1) = (C1/dt) * ( beam%theta(0:ns-1) &
                         - C3*beam%theta_old(0:ns-1) &
                         - C4*beam%theta_oldold(0:ns-1) ) &
                         - C2*beam%theta_dot_old(0:ns-1)

  !-----------------------------------------------------------------------------
  ! get deflection line (integrate_position)
  !-----------------------------------------------------------------------------
  ! this beam is the new one, so its time+dt ( to get mouvement at the right instant)
  call integrate_position ( time+dt, beam )

  !-----------------------------------------------------------------------------
  ! accelerations
  !-----------------------------------------------------------------------------
  ! we have the old velocity (t_n) and the new one
  beam%ax(0:ns-1) = (beam%vx(0:ns-1) - beam_old%vx(0:ns-1)) / dt
  beam%ay(0:ns-1) = (beam%vy(0:ns-1) - beam_old%vy(0:ns-1)) / dt

  !-----------------------------------------------------------------------------
  ! emergency brake
  !-----------------------------------------------------------------------------
  ! if ((maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0 ).and.(root)) then
  if ((maxval(abs(beam%theta_dot(0:ns-1)-beam%theta_dot_old(0:ns-1) ))>100.d0 ).and.(root)) then
    write (*,'(A)') "!!! IBES-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0"
    write (*,'(A)') "    That indicates a possible instability."
    write (*,'("time=",es11.4)') time
  endif

  !-----------------------------------------------------------------------------
  ! save number of iterations
  !-----------------------------------------------------------------------------
  if (root) then
    open (14, file = 'IBES_iter.t', status = 'unknown', position='append')
    write (14, '(es11.4,1x,i3)') time, iter
    close (14)
  endif

  !-----------------------------------------------------------------------------
  ! change marching scheme, if necessary
  !-----------------------------------------------------------------------------
  if (ActuallyBDF2 .eqv. .true.) then   ! Remember to switch back to BDF2
    TimeMethodSolid = "BDF2"
    ActuallyBDF2 = .false.
  endif

  ! we advanced the beam to time level (n+1). if we call this routine the next time
  ! our new solution will be "old" and what is now "old" will become "oldold"
  ! so iterate that now:
  beam%theta_oldold(0:ns-1) = beam%theta_old(0:ns-1)
  beam%theta_dot_oldold(0:ns-1) = beam%theta_dot_old(0:ns-1)

  ! for BDF2 with variable dt, we need the old time step
  ! this should be okay also when restarting, as the first CN2 step will provide us with the dt_old
  beam%dt_old = dt

end subroutine IBES_solver




!-------------------------------------------------------------------------------
! Non-linear funtion F(x), which is our root problem F(x)=0
! This new version is obtained using a lot of different tools and the equations
! should be exactly as they are written in my PHD thesis
!-------------------------------------------------------------------------------
subroutine F_nonlinear (time, dt, F, x, beam_solid)

  implicit none
  real(kind=pr),intent (in) :: time, dt
  type(solid), intent(in) :: beam_solid
  real(kind=pr),dimension(1:2*ns+4),intent(out) :: F
  real(kind=pr),dimension(1:2*ns+4),intent(in)  :: x

  real(kind=pr),dimension(0:ns-1)  :: theta_old, theta_dot_old
  real(kind=pr),dimension(0:ns-1)  :: old_rhs, theta_oldold, theta_dot_oldold
  real(kind=pr),dimension(-1:ns+1) :: theta
  real(kind=pr),dimension(-1:ns-1) :: T
  real(kind=pr),dimension(0:ns-1)  :: p, tau_beam, tau_s, p_s, K3, K4
  real(kind=pr),dimension(0:ns-1)  :: zeta, zeta_s, zeta_ss, zeta_sss, mu, mu_s, mu_star

  real(kind=pr) :: K1,K2,C1,C2,C3,C4,R
  real(kind=pr) :: alpha, alpha_t, alpha_tt
  real(kind=pr), dimension(1:6) :: LeadingEdge ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  integer :: i

  ! as we iterate to find theta^(n+1), also external loads are at (n+1)
  tau_beam = beam_solid%tau_new(0:ns-1)
  p        = beam_solid%pressure_new(0:ns-1)

  call time_marching_coefs(dt,beam_solid%dt_old,C1,C2,C3,C4)
  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )

  call Differentiate1D (tau_beam, tau_s, ns, ds, 1)
  call Differentiate1D (p, p_s, ns, ds, 1)

  ! shorten coding with local variables:

  ! these guys are what we iterate, so theta^(n+1) and T^(n+1)
  theta = x(1:ns+3)
  T     = x(ns+4:2*ns+4)

  ! old time levels:
  theta_old = beam_solid%theta(0:ns-1)
  theta_oldold = beam_solid%theta_oldold(0:ns-1)
  theta_dot_old = beam_solid%theta_dot(0:ns-1)
  theta_dot_oldold = beam_solid%theta_dot_oldold(0:ns-1)
  ! currently not implemented:
  old_rhs = 0.d0

  zeta     = beam_solid%zeta(0:ns-1)
  zeta_s   = beam_solid%zeta_s(0:ns-1)
  zeta_ss  = beam_solid%zeta_ss(0:ns-1)
  zeta_sss = beam_solid%zeta_sss(0:ns-1)
  mu       = beam_solid%mu(0:ns-1)
  mu_s     = beam_solid%mu_s(0:ns-1)
  mu_star  = beam_solid%mu_star(0:ns-1)

  !---------new extended boundarys (constants)
  K1 = mu(0)*( LeadingEdge(5)*dcos(alpha)+LeadingEdge(6)*dsin(alpha)+grav*dsin(alpha)) - tau_beam(0)
  K2 = mu(0)*(-LeadingEdge(5)*dsin(alpha)+LeadingEdge(6)*dcos(alpha)+grav*dcos(alpha)) + p(0)
  K3 = C3*theta_dot_old + C4*theta_dot_oldold + (dt/C1)*C2*old_rhs
  K4 = K3 + C2*theta_dot_old + (C1/dt)*(C3*theta_old+C4*theta_oldold)

  F(1) = theta(0)
  F(2) = (T(1) - T(-1))/(2*ds) - K1 + (zeta_s(0)*(theta(1) - theta(-1))**2)/(4*ds**2) &
       + (zeta(0)*(theta(1) - theta(-1))*(theta(1) - 2*theta(0) + theta(-1)))/(2*ds**3)
  F(3) = (T(0)*(theta(1) - theta(-1)))/(2*ds) - K2 - (zeta_ss(0)*(theta(1) &
       - theta(-1)))/(2*ds) - (2*zeta_s(0)*(theta(1) - 2*theta(0) + theta(-1)))/ds**2 &
       + (zeta(0)*(12*theta(1) - 10*theta(0) - 6*theta(2) + theta(3) + 3*theta(-1)))/(2*ds**3)
  F(4) = T(ns-1)
  F(5) = (zeta_s(ns-1)*(theta(ns) - theta(ns-2)))/(2*ds) - (zeta(ns-1)*((5*theta(ns-1))/2 &
       - (4*theta(ns))/3 - (4*theta(ns-2))/3 + theta(ns-3)/12 + theta(ns+1)/12))/ds**2
  F(6) = (zeta(ns-1)*((2*theta(ns))/3 - (2*theta(ns-2))/3 + theta(ns-3)/12 - theta(ns+1)/12))/ds

  i = 0
  F(7) = tau_s(i) + mu(i)*(C2*theta_dot_old(i) - alpha_t + (C1*(C3*theta_old(i) &
       - theta(i) + C4*theta_oldold(i)))/dt)**2 - mu_star(i)*(tau_beam(i) &
       + (zeta_s(i)*(theta(i-1)/2 - theta(i+1)/2)**2)/ds**2 - (zeta(i)*(theta(i-1)/2 &
       - theta(i+1)/2)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**3) &
       + (T(i-1) - 2*T(i) + T(i+1))/ds**2 + (zeta(i)*(theta(i-1) - 2*theta(i) &
       + theta(i+1))**2)/ds**4 + (mu_star(i)*(T(i-1)/2 - T(i+1)/2))/ds &
       - (p(i)*(theta(i-1)/2 - theta(i+1)/2))/ds - (T(i)*(theta(i-1)/2 &
       - theta(i+1)/2)**2)/ds**2 + (2*zeta_ss(i)*(theta(i-1)/2 - theta(i+1)/2)**2)/ds**2 &
       + (2*zeta(i)*(theta(i-1)/2 - theta(i+1)/2)*((5*theta(i))/2 - 9*theta(i+1) &
       + 12*theta(i+2) - 7*theta(i+3) + (3*theta(i+4))/2))/ds**4 - (5*zeta_s(i)*(theta(i-1)/2 &
       - theta(i+1)/2)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**3

  i = ns-1
  F(8) = K4(i) - (C1*theta(i))/dt - (dt*(alpha_tt*mu(i)**2 + p_s(i) - mu_star(i)*(p(i) &
      + (T(i)*(theta(i-1)/2 - theta(i+1)/2))/ds - (zeta_ss(i)*(theta(i-1)/2 - theta(i+1)/2))/ds &
      + (2*zeta_s(i)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**2 + (zeta(i)*(theta(i-1) &
      - theta(i-2)/2 - theta(i+1) + theta(i+2)/2))/ds**3) + (zeta(i)*(6*theta(i) &
      - 4*theta(i-1) + theta(i-2) - 4*theta(i+1) + theta(i+2)))/ds**4 &
      - (zeta_sss(i)*(theta(i-1)/2 - theta(i+1)/2))/ds - (T(i)*(theta(i-1) - 2*theta(i) &
      + theta(i+1)))/ds**2 + (tau_beam(i)*(theta(i-1)/2 - theta(i+1)/2))/ds &
      + (3*zeta_ss(i)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**2 &
      + (2*(theta(i-1)/2 - theta(i+1)/2)*((3*T(i))/2 + 2*T(i-1) - T(i-2)/2))/ds**2 &
      + (3*zeta_s(i)*(theta(i-1) - theta(i-2)/2 - theta(i+1) + theta(i+2)/2))/ds**3 &
      + (zeta_s(i)*(theta(i-1)/2 - theta(i+1)/2)**3)/ds**3 - (zeta(i)*(theta(i-1)/2 &
      - theta(i+1)/2)**2*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**4))/(C1*mu(i))

  do i=1, ns-2
    ! first block: eqn for theta (ns+2 eqn's) [F]
    F( 8+i ) = K4(i) - (C1*theta(i))/dt - (dt*(alpha_tt*mu(i)**2 + p_s(i) - mu_star(i)*(p(i) &
    + (T(i)*(theta(i-1)/2 - theta(i+1)/2))/ds - (zeta_ss(i)*(theta(i-1)/2 - theta(i+1)/2))/ds &
    + (2*zeta_s(i)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**2 + (zeta(i)*(theta(i-1) &
    - theta(i-2)/2 - theta(i+1) + theta(i+2)/2))/ds**3) + (zeta(i)*(6*theta(i) - 4*theta(i-1) &
    + theta(i-2) - 4*theta(i+1) + theta(i+2)))/ds**4 - (zeta_sss(i)*(theta(i-1)/2 - theta(i+1)/2))/ds &
    - (T(i)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**2 + (tau_beam(i)*(theta(i-1)/2 &
    - theta(i+1)/2))/ds + (3*zeta_ss(i)*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**2 &
    + (3*zeta_s(i)*(theta(i-1) - theta(i-2)/2 - theta(i+1) + theta(i+2)/2))/ds**3 &
    + (zeta_s(i)*(theta(i-1)/2 - theta(i+1)/2)**3)/ds**3 - (2*(theta(i-1)/2 &
    - theta(i+1)/2)*(T(i-1)/2 - T(i+1)/2))/ds**2 - (zeta(i)*(theta(i-1)/2 &
    - theta(i+1)/2)**2*(theta(i-1) - 2*theta(i) + theta(i+1)))/ds**4))/(C1*mu(i))

    ! 2nd block: eqn's for T (ns+2 eqn's) [G]
    F( 8+(ns-2)+i ) = tau_s(i) + mu(i)*(C2*theta_dot_old(i) - alpha_t + (C1*(C3*theta_old(i) - theta(i) &
    + C4*theta_oldold(i)))/dt)**2 - mu_star(i)*(tau_beam(i) + (zeta_s(i)*(theta(i-1)/2 &
    - theta(i+1)/2)**2)/ds**2 - (zeta(i)*(theta(i-1)/2 - theta(i+1)/2)*(theta(i-1) &
    - 2*theta(i) + theta(i+1)))/ds**3) + (T(i-1) - 2*T(i) + T(i+1))/ds**2 + (zeta(i)*(theta(i-1) &
    - 2*theta(i) + theta(i+1))**2)/ds**4 + (mu_star(i)*(T(i-1)/2 - T(i+1)/2))/ds - (p(i)*(theta(i-1)/2 &
    - theta(i+1)/2))/ds - (T(i)*(theta(i-1)/2 - theta(i+1)/2)**2)/ds**2 + (2*zeta_ss(i)*(theta(i-1)/2 &
    - theta(i+1)/2)**2)/ds**2 - (5*zeta_s(i)*(theta(i-1)/2 - theta(i+1)/2)*(theta(i-1) - 2*theta(i) &
    + theta(i+1)))/ds**3 - (2*zeta(i)*(theta(i-1)/2 - theta(i+1)/2)*(theta(i-1) - theta(i-2)/2 - theta(i+1) &
    + theta(i+2)/2))/ds**4
  enddo
end subroutine F_nonlinear


!-------------------------------------------------------------------------------


subroutine Jacobi(time, dt, J, x, beam_solid)
  implicit none

  real(kind=pr),intent(in) :: time, dt
  type(solid)  ,intent(in) :: beam_solid
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real(kind=pr),dimension(1:2*ns+4),intent(inout)  :: x

  real(kind=pr),dimension(0:ns-1)  :: theta_old, theta_dot_old
  real(kind=pr),dimension(0:ns-1)  :: old_rhs, theta_oldold, theta_dot_oldold
  real(kind=pr),dimension(-1:ns+1) :: theta
  real(kind=pr),dimension(-1:ns-1) :: T
  real(kind=pr),dimension(0:ns-1)  :: p, tau_beam, tau_s, p_s, K3, K4
  real(kind=pr),dimension(0:ns-1)  :: zeta, zeta_s, zeta_ss, zeta_sss, mu, mu_s, mu_star

  real(kind=pr) :: alpha, alpha_t, alpha_tt, C1,C2,C3,C4,R, K1,K2
  real(kind=pr), dimension(1:6) :: LeadingEdge
  integer :: i, T0_index, k, l, m, n
  real(kind=pr),dimension(1:2*ns+4,1:2*ns+4) :: J2,J2_norm
  integer, save :: iCalls=10

  J = 0.d0                ! initialize J

  ! as we iterate to find theta^(n+1), also external loads are at (n+1)
  tau_beam = beam_solid%tau_new(0:ns-1)
  p        = beam_solid%pressure_new(0:ns-1)

  call time_marching_coefs(dt,beam_solid%dt_old,C1,C2,C3,C4)
  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )

  call Differentiate1D (tau_beam, tau_s, ns, ds, 1)
  call Differentiate1D (p, p_s, ns, ds, 1)

  ! shorten coding with local variables:

  ! these guys are what we iterate, so theta^(n+1) and T^(n+1)
  theta = x(1:ns+3)
  T     = x(ns+4:2*ns+4)

  ! old time levels:
  theta_old = beam_solid%theta(0:ns-1)
  theta_oldold = beam_solid%theta_oldold(0:ns-1)
  theta_dot_old = beam_solid%theta_dot(0:ns-1)
  theta_dot_oldold = beam_solid%theta_dot_oldold(0:ns-1)
  ! currently not implemented:
  old_rhs = 0.d0

  zeta     = beam_solid%zeta(0:ns-1)
  zeta_s   = beam_solid%zeta_s(0:ns-1)
  zeta_ss  = beam_solid%zeta_ss(0:ns-1)
  zeta_sss = beam_solid%zeta_sss(0:ns-1)
  mu       = beam_solid%mu(0:ns-1)
  mu_s     = beam_solid%mu_s(0:ns-1)
  mu_star  = beam_solid%mu_star(0:ns-1)


  !---------new extended boundarys (constants)
  K1 = mu(0)*( LeadingEdge(5)*dcos(alpha)+LeadingEdge(6)*dsin(alpha)+grav*dsin(alpha)) - tau_beam(0)
  K2 = mu(0)*(-LeadingEdge(5)*dsin(alpha)+LeadingEdge(6)*dcos(alpha)+grav*dcos(alpha)) + p(0)
  K3 = C3*theta_dot_old + C4*theta_dot_oldold
  K4 = K3 + C2*theta_dot_old + (C1/dt)*(C3*theta_old+C4*theta_oldold)



  !---------------------------------------
  ! ---indexing:
  ! theta natural index -1...ns+1
  ! dF/dtheta_i begins @ 1
  !    dTheta_-1   = 1      ghost node LE
  !    dTheta_0    = 2      first regular
  !    dTheta_1    = 3
  !    dTheta_ns-3 = ns-1
  !    dTheta_ns-2 = ns
  !    dTheta_ns-1 = ns+1   last regular
  !    dTheta_ns   = ns+2   ghost1 TE
  !    dTheta_ns+1 = ns+3   ghost2 TE

  !    dT_-1       = ns+4
  !    dT_0        = ns+5
  !    dT_1        = ns+6
  !    dT_ns-2     =2ns+3
  !    dT_ns-1     =2ns+4
  !---------------------------------------

  !-- set first 6 eqns
  J(2,1) = 1.0d0

  !----eqn 2 (special BC 1)
  ! dF(2) / dtheta(-1)
  J(1,2) = (zeta(0)*(theta(0) - theta(-1)))/ds**3 - (zeta_s(0)*(theta(1) - theta(-1)))/(2*ds**2)
  ! dF(2) / dtheta(0)
  J(2,2) = -(zeta(0)*(theta(1) - theta(-1)))/ds**3
  ! dF(2) / dtheta(1)
  J(3,2) = (zeta_s(0)*(theta(1) - theta(-1)))/(2*ds**2) - (zeta(0)*(theta(0) - theta(1)))/ds**3
  ! dF(2) / dT(-1)
  J(ns+4,2) = -1/(2*ds)
  ! dF(2) / dT(1)
  J(ns+6,2) = 1/(2*ds)



  !----eqn 3 (special BC 2)
  ! dF(3) / dtheta(-1)
  J(1,3) = -(2*ds*zeta_s(0) - (3*zeta(0))/2 + ds**2*(T(0)/2 - zeta_ss(0)/2))/ds**3
  ! dF(3) / dtheta(0)
  J(2,3) = -(5*zeta(0) - 4*ds*zeta_s(0))/ds**3
  ! dF(3) / dtheta(1)
  J(3,3) = (6*zeta(0) - 2*ds*zeta_s(0) + ds**2*(T(0)/2 - zeta_ss(0)/2))/ds**3
  ! dF(3) / dtheta(2)
  J(4,3) = -(3*zeta(0))/ds**3
  ! dF(3) / dtheta(3)
  J(5,3) = zeta(0)/(2*ds**3)
  ! dF(3) / dT(0)
  J(ns+5,3) = (theta(1) - theta(-1))/(2*ds)

  !-----eqn 4
  ! dF4 / dT(ns-1)
  J(2*ns+4,4) = 1.d0

  !-----eqn 5
  ! dF(5) / dtheta(ns-3)
  J(ns-1,5) = -zeta(ns-1)/(12*ds**2)
  ! dF(5) / dtheta(ns-2)
  J(ns,5) = (8*zeta(ns-1) - 3*ds*zeta_s(ns-1))/(6*ds**2)
  ! dF(5) / dtheta(ns-1)
  J(ns+1,5) = -(5*zeta(ns-1))/(2*ds**2)
  ! dF(5) / dtheta(ns)
  J(ns+2,5) = (8*zeta(ns-1) + 3*ds*zeta_s(ns-1))/(6*ds**2)
  ! dF(5) / dtheta(ns+1)
  J(ns+3,5) = -zeta(ns-1)/(12*ds**2)


  ! dF(6) / dtheta(ns-3)
  J(ns-1,6) = zeta(ns-1)/(12*ds)
  ! dF(6) / dtheta(ns-2)
  J(ns  ,6) = -(2*zeta(ns-1))/(3*ds)
  ! dF(6) / dtheta(ns-1)
  J(ns+1,6) = 0
  ! dF(6) / dtheta(ns)
  J(ns+2,6) = (2*zeta(ns-1))/(3*ds)
  ! dF(6) / dtheta(ns+1)
  J(ns+3,6) = -zeta(ns-1)/(12*ds)

  i = 0
  ! dF(7) / dtheta(i-1)
  J(1,7) = ((5*zeta_s(i) - zeta(i)*mu_star(i))*(theta(i) - theta(i-1)))/ds**3 &
         - p(i)/(2*ds) - ((theta(i-1) - theta(i+1))*(T(i) - 2*zeta_ss(i) + zeta_s(i)*mu_star(i)))/(2*ds**2) &
         - (zeta(i)*(3*theta(i) - 4*theta(i-1) + 14*theta(i+1) - 24*theta(i+2) + 14*theta(i+3) &
         - 3*theta(i+4)))/(2*ds**4)
  ! dF(7) / dtheta(i)
  J(2,7) = ((5*zeta_s(i) - zeta(i)*mu_star(i))*(theta(i-1) - theta(i+1)))/ds**3 &
         - (zeta(i)*(3*theta(i-1) - 16*theta(i) + 13*theta(i+1)))/(2*ds**4) - (2*C1*mu(i)*(C1*C3*theta_old(i) &
         - alpha_t*dt - C1*theta(i) + C1*C4*theta_oldold(i) + C2*dt*theta_dot_old(i)))/dt**2
  ! dF(7) / dtheta(i+1)
  J(3,7) = p(i)/(2*ds) - ((5*zeta_s(i) - zeta(i)*mu_star(i))*(theta(i) - theta(i+1)))/ds**3 &
         + ((theta(i-1) - theta(i+1))*(T(i) - 2*zeta_ss(i) + zeta_s(i)*mu_star(i)))/(2*ds**2) &
         - (zeta(i)*(13*theta(i) + 14*theta(i-1) - 40*theta(i+1) + 24*theta(i+2) - 14*theta(i+3) &
         + 3*theta(i+4)))/(2*ds**4)
  ! dF(7) / dtheta(i+2)
  J(4,7) = (12*zeta(i)*(theta(i-1) - theta(i+1)))/ds**4
  ! dF(7) / dtheta(i+3)
  J(5,7) = -(7*zeta(i)*(theta(i-1) - theta(i+1)))/ds**4
  ! dF(7) / dtheta(i+4)
  J(6,7) = (3*zeta(i)*(theta(i-1) - theta(i+1)))/(2*ds**4)
  ! dF(7) / dT(i-1)
  J(ns+4,7) = (ds*mu_star(i) + 2)/(2*ds**2)
  ! dF(7) / dT(i)
  J(ns+5,7) = - (theta(i-1) - theta(i+1))**2/(4*ds**2) - 2/ds**2
  ! dF(7) / dT(i+1)
  J(ns+6,7) = -(ds*mu_star(i) - 2)/(2*ds**2)


  i = ns-1
  ! dF(8) / dtheta(i-2)
  J(ns-1,8) = -(dt*zeta(i) - (ds*dt*(3*zeta_s(i) - zeta(i)*mu_star(i)))/2)/(C1*ds**4*mu(i))
  ! dF(8) / dtheta(i-1)
  J(ns  ,8) = (dt*(zeta_sss(i) - tau_beam(i) + T(i)*mu_star(i) - zeta_ss(i)*mu_star(i)))/(2*C1*ds*mu(i)) &
            - (dt*(24*zeta_s(i) - 8*zeta(i)*mu_star(i) + 3*zeta_s(i)*theta(i-1)**2 + 3*zeta_s(i)*theta(i+1)**2 &
            - 6*zeta_s(i)*theta(i-1)*theta(i+1)))/(8*C1*ds**3*mu(i)) - (dt*(T(i) + 4*T(i-1) - T(i-2) &
            + 6*zeta_ss(i) - 4*zeta_s(i)*mu_star(i)))/(2*C1*ds**2*mu(i)) - (dt*zeta(i)*(4*theta(i)*theta(i-1) &
            - 4*theta(i)*theta(i+1) + 2*theta(i-1)*theta(i+1) - 3*theta(i-1)**2 + theta(i+1)**2 &
            - 16))/(4*C1*ds**4*mu(i))
  ! dF(8) / dtheta(i)
  J(ns+1,8) = - C1/dt - (dt*(6*zeta(i) + (zeta(i)*theta(i-1)**2)/2 + (zeta(i)*theta(i+1)**2)/2 &
            - zeta(i)*theta(i-1)*theta(i+1)) + ds**2*dt*(2*T(i) - 6*zeta_ss(i) &
            + 4*zeta_s(i)*mu_star(i)))/(C1*ds**4*mu(i))
  ! dF(8) / dtheta(i+1)
  J(ns+2,8) = ((dt*(32*zeta(i) - 2*zeta(i)*theta(i-1)**2 + 6*zeta(i)*theta(i+1)**2 &
            + 8*zeta(i)*theta(i)*theta(i-1) - 8*zeta(i)*theta(i)*theta(i+1) &
            - 4*zeta(i)*theta(i-1)*theta(i+1)))/8 + (ds*dt*(24*zeta_s(i) - 8*zeta(i)*mu_star(i) &
            + 3*zeta_s(i)*theta(i-1)**2 + 3*zeta_s(i)*theta(i+1)**2 - 6*zeta_s(i)*theta(i-1)*theta(i+1)))/8 &
            - (ds**3*dt*(4*zeta_sss(i) - 4*tau_beam(i) + 4*T(i)*mu_star(i) - 4*zeta_ss(i)*mu_star(i)))/8 &
            + (ds**2*dt*(20*T(i) + 16*T(i-1) - 4*T(i-2) - 24*zeta_ss(i) &
            + 16*zeta_s(i)*mu_star(i)))/8)/(C1*ds**4*mu(i))
  ! dF(8) / dtheta(i+2)
  J(ns+3,8) = -(dt*zeta(i) + (ds*dt*(3*zeta_s(i) - zeta(i)*mu_star(i)))/2)/(C1*ds**4*mu(i))

  ! dF(8) / dT(i-2)
  J(2*ns+2,8) = (dt*(theta(i-1) - theta(i+1)))/(2*C1*ds**2*mu(i))
  ! dF(8) / dT(i-1)
  J(2*ns+3,8) = -(2*dt*(theta(i-1) - theta(i+1)))/(C1*ds**2*mu(i))
  ! dF(8) / dT(i)
  J(2*ns+4,8) = -((dt*(4*theta(i) + theta(i-1) - 5*theta(i+1)))/2 - (ds*dt*(mu_star(i)*theta(i-1) &
              - mu_star(i)*theta(i+1)))/2)/(C1*ds**2*mu(i))

  ! loop over regular points for the evolution equation
  ! note index -1 is the ghost node for theta, then comes theta(0) which is dirichlet
  ! therefore, loop starts at i=1 (and no special treatment of the first point, enough ghosts)
  ! the modified point is (ns-1), which is eqn F8
  do i = 1, ns-2
    k = 8 + i !index for eqns F (row) -- that means first row is row 9
    l = i + 2 ! index of theta(i)
    m = i + ns + 5 !i index of T(i)

    ! for the first point, we derive equation 9 with respect to theta(1) and T(1)
    ! (and their neighboring points)

    ! dF( 8+i ) / dtheta(i-2)
    J(l-2,k) = -(dt*zeta(i) - (ds*dt*(3*zeta_s(i) - zeta(i)*mu_star(i)))/2)/(C1*ds**4*mu(i))
    ! dF( 8+i ) / dtheta(i-1)
    J(l-1,k) = (dt*(2*T(i) + T(i-1) - T(i+1) - 6*zeta_ss(i) + 4*zeta_s(i)*mu_star(i)))/(2*C1*ds**2*mu(i)) &
             - (dt*(24*zeta_s(i) - 8*zeta(i)*mu_star(i) + 3*zeta_s(i)*theta(i-1)**2 + 3*zeta_s(i)*theta(i+1)**2 &
             - 6*zeta_s(i)*theta(i-1)*theta(i+1)))/(8*C1*ds**3*mu(i)) + (dt*(zeta_sss(i) - tau_beam(i) &
             + T(i)*mu_star(i) - zeta_ss(i)*mu_star(i)))/(2*C1*ds*mu(i)) - (dt*zeta(i)*(4*theta(i)*theta(i-1) &
             - 4*theta(i)*theta(i+1) + 2*theta(i-1)*theta(i+1) - 3*theta(i-1)**2 + theta(i+1)**2 &
             - 16))/(4*C1*ds**4*mu(i))
    ! dF( 8+i ) / dtheta(i)
    J(l,k) = - C1/dt - (dt*(6*zeta(i) + (zeta(i)*theta(i-1)**2)/2 + (zeta(i)*theta(i+1)**2)/2 &
           - zeta(i)*theta(i-1)*theta(i+1)) + ds**2*dt*(2*T(i) - 6*zeta_ss(i) &
           + 4*zeta_s(i)*mu_star(i)))/(C1*ds**4*mu(i))
    ! dF( 8+i ) / dtheta(i+1)
    J(l+1,k) = (dt*(2*T(i) - T(i-1) + T(i+1) - 6*zeta_ss(i) + 4*zeta_s(i)*mu_star(i)))/(2*C1*ds**2*mu(i)) &
             + (dt*(24*zeta_s(i) - 8*zeta(i)*mu_star(i) + 3*zeta_s(i)*theta(i-1)**2 + 3*zeta_s(i)*theta(i+1)**2 &
             - 6*zeta_s(i)*theta(i-1)*theta(i+1)))/(8*C1*ds**3*mu(i)) - (dt*(zeta_sss(i) - tau_beam(i) &
             + T(i)*mu_star(i) - zeta_ss(i)*mu_star(i)))/(2*C1*ds*mu(i)) + (dt*zeta(i)*(4*theta(i)*theta(i-1) &
             - 4*theta(i)*theta(i+1) - 2*theta(i-1)*theta(i+1) - theta(i-1)**2 + 3*theta(i+1)**2 &
             + 16))/(4*C1*ds**4*mu(i))
    ! dF( 8+i ) / dtheta(i+2)
    J(l+2,k) = -(dt*zeta(i) + (ds*dt*(3*zeta_s(i) - zeta(i)*mu_star(i)))/2)/(C1*ds**4*mu(i))

    ! dF( 8+i ) / dT(i-1)
    J(m-1,k) = (dt*(theta(i-1) - theta(i+1)))/(2*C1*ds**2*mu(i))
    ! dF( 8+i ) / dT(i)
    J(m,k) = (dt*((theta(i-1) - 2*theta(i) + theta(i+1))/ds**2 + (mu_star(i)*(theta(i-1)/2 &
           - theta(i+1)/2))/ds))/(C1*mu(i))
    ! dF( 8+i ) / dT(i+1)
    J(m+1,k) = -(dt*(theta(i-1) - theta(i+1)))/(2*C1*ds**2*mu(i))
  enddo


  ! loop over regular points for tension equation
  ! point(0) is the modified stencil (eqn F7)
  do i = 1, ns-2
    k = 8 + ns - 2 + i !index for eqns F (row) !starts @ ns+7
    l = i + 2 !index column (theta)
    m = i + ns + 5 !index column (T) -- starts @ T(-1) = ns+5

    ! dF( 8+(ns-2)+i ) / dtheta(i-2)
    J(l-2,k) = (zeta(i)*(theta(i-1) - theta(i+1)))/(2*ds**4)
    ! dF( 8+(ns-2)+i ) / dtheta(i-1)
    J(l-1,k) = ((5*zeta_s(i) - zeta(i)*mu_star(i))*(theta(i) - theta(i-1)))/ds**3 - p(i)/(2*ds) &
            - ((theta(i-1) - theta(i+1))*(T(i) - 2*zeta_ss(i) + zeta_s(i)*mu_star(i)))/(2*ds**2) &
            - (zeta(i)*(8*theta(i) - theta(i-2) - 8*theta(i+1) + theta(i+2)))/(2*ds**4)
    ! dF( 8+(ns-2)+i ) / dtheta(i)
    J(l,k) = (10*zeta_s(i)*(theta(i-1)/2 - theta(i+1)/2))/ds**3 - (zeta(i)*(4*theta(i-1) - 8*theta(i) &
           + 4*theta(i+1)))/ds**4 - (2*C1*mu(i)*(C2*theta_dot_old(i) - alpha_t + (C1*(C3*theta_old(i) &
           - theta(i) + C4*theta_oldold(i)))/dt))/dt - (2*zeta(i)*mu_star(i)*(theta(i-1)/2 - theta(i+1)/2))/ds**3
    ! dF( 8+(ns-2)+i ) / dtheta(i+1)
    J(l+1,k) = p(i)/(2*ds) - ((5*zeta_s(i) - zeta(i)*mu_star(i))*(theta(i) - theta(i+1)))/ds**3 &
             + ((theta(i-1) - theta(i+1))*(T(i) - 2*zeta_ss(i) + zeta_s(i)*mu_star(i)))/(2*ds**2) &
             - (zeta(i)*(8*theta(i) - 8*theta(i-1) + theta(i-2) - theta(i+2)))/(2*ds**4)
    ! dF( 8+(ns-2)+i ) / dtheta(i+2)
    J(l+2,k) = -(zeta(i)*(theta(i-1) - theta(i+1)))/(2*ds**4)

    ! dF( 8+(ns-2)+i ) / dT(i-1)
    J(m-1,k) = (ds*mu_star(i) + 2)/(2*ds**2)
    ! dF( 8+(ns-2)+i ) / dT(i)
    J(m,k) = - (theta(i-1) - theta(i+1))**2/(4*ds**2) - 2/ds**2
    ! dF( 8+(ns-2)+i ) / dT(i+1)
    J(m+1,k) = -(ds*mu_star(i) - 2)/(2*ds**2)

  enddo


  !-------------------------------------------------------------------------
  !  self-test. occasionally, check if Jacobian is okay
  !-------------------------------------------------------------------------
  if (mod(iCalls,1967)==0) then
    call Jacobi_num(time, dt, J2, x, beam_solid)
    J2_norm = J2
    where (abs(J2_norm)<1.0d-7) J2_norm=1.0d0
    J2_norm = (J-J2)/J2_norm
    where (abs(J2_norm)<1.0d-5) J2_norm=0.d0 ! delete small values
    if (maxval(J2_norm)>1.0d-2) then
        open (14, file = 'IBES_JACOBIAN', status = 'replace')
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
        write (*,*) "time=", time
        close (14)
        call abort(1231, "IBES error, mismatch in Jacobian (numeric vs analytic)")
    endif
    iCalls = 0
  endif

  iCalls = iCalls + 1      ! count calls (to perform rare self-tests)
end subroutine Jacobi


!-------------------------------------------------------------------------------

! NOT YET UPDATED TO NON-CONST COEFFICIENTS
subroutine RHS_beameqn (time, theta, theta_dot, pressure_beam, T, tau_beam, beam_solid)

    implicit none
    real(kind=pr), intent (in) ::         time
    real(kind=pr) ::              A1, A2, K2,C2
    real(kind=pr),dimension(0:ns-1), intent (in) :: pressure_beam, tau_beam
    real(kind=pr),dimension(0:ns-1), intent (inout) :: T
    type(solid) :: beam_solid
    real(kind=pr),dimension(0:ns-1), intent (inout) ::   theta, theta_dot
    real(kind=pr),dimension(0:ns+2) :: theta_extended, theta_extended_s, theta_extended_ss
    real(kind=pr),dimension(0:ns+2) :: theta_extended_sss, theta_extended_ssss
    real(kind=pr),dimension(0:ns-1) :: theta_s, theta_ss, theta_sss, theta_ssss, T_s, p_s
    real(kind=pr) :: alpha, alpha_t, alpha_tt
    real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)


end subroutine RHS_beameqn


!-------------------------------------------------------------------------------
! Tension
! Compute tension inside the beam for a given theta, theta_dot and external loads
! equation is valid for non-const eta, mu, we updated the code for readability:
! as in the analoguous matlab script, we just solve A*x=b (linear system)
! and A includes both BC (in previous versions, we insisted on a triag system)
! INPUT:
!   theta, theta_dot: angles and angular velocity of beam at current state
!   pressure, tau_beam: external forces
!   beam_solid: compatibility
! OUTPUT:
!   T, T_s: tension and tension derivative (note: BOTH WITH GHOSTNODES!)
!-------------------------------------------------------------------------------
subroutine Tension ( time, T, T_s, theta, theta_dot, pressure, tau_beam, beam_solid)
  implicit none
  real(kind=pr), intent (in) :: time
  real(kind=pr),dimension(0:ns-1) :: theta_s, theta_ss, theta_sss
  real(kind=pr),dimension(0:ns-1) :: tau_beam_s
  real(kind=pr),dimension(0:ns-1), intent (in) ::  theta, theta_dot, pressure, tau_beam
  type(solid) :: beam_solid
  real(kind=pr):: K1,K2
  real(kind=pr):: alpha, alpha_t, alpha_tt
  real(kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real(kind=pr),dimension(-1:ns-1,-1:ns-1) :: A_matrix
  real(kind=pr),dimension(-1:ns-1) :: rhs
  real(kind=pr),dimension(-1:ns-1), intent (out) :: T, T_s
  real(kind=pr),dimension(0:ns-1) :: ipiv !used only for the MKL lib
  real(kind=pr),dimension(0:ns-1) :: zeta, zeta_s, zeta_ss, zeta_sss, mu, mu_s, mu_star
  integer :: i, info

  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beam_solid )

  call Differentiate1D (theta, theta_s, ns, ds, 1)
  call Differentiate1D (theta, theta_ss, ns, ds, 2)
  call Differentiate1D (theta, theta_sss, ns, ds, 3)
  call Differentiate1D (tau_beam, tau_beam_s, ns, ds, 1)

  zeta     = beam_solid%zeta(0:ns-1)
  zeta_s   = beam_solid%zeta_s(0:ns-1)
  zeta_ss  = beam_solid%zeta_ss(0:ns-1)
  zeta_sss = beam_solid%zeta_sss(0:ns-1)
  mu      = beam_solid%mu(0:ns-1)
  mu_s    = beam_solid%mu_s(0:ns-1)
  mu_star = beam_solid%mu_star(0:ns-1)

  ! constants from the boundary conditions
  K1 = mu(0)*( LeadingEdge(5)*dcos(alpha)+LeadingEdge(6)*dsin(alpha)+grav*dsin(alpha)) - tau_beam(0)
  K2 = mu(0)*(-LeadingEdge(5)*dsin(alpha)+LeadingEdge(6)*dcos(alpha)+grav*dcos(alpha)) + pressure(0)

  !----------------------------------------------------------------
  ! first step: construct the opertator matrix A
  !----------------------------------------------------------------
  A_matrix = 0.d0

  ! note we skip first line (-1), we'll set that later (Neumann BC)
  ! note point (ns-1) is a dirichlet BC so we skip that as well
  do i = 0, ns-2
    A_matrix(i,i-1) = +1.d0/ds**2 - 1.d0/(2.d0*ds) * mu_star(i)
    A_matrix(i,i)   = -2.d0/ds**2 - theta_s(i)**2
    A_matrix(i,i+1) = +1.d0/ds**2 + 1.d0/(2.d0*ds) * mu_star(i)
  enddo

  ! Dirichlet condition (trailing edge):
  A_matrix(ns-1,ns-1) = 1.d0

  ! neumann condition (leading edge)
  A_matrix(-1,-1) = -1.d0/(2.d0*ds)
  A_matrix(-1,1) = +1.d0/(2.d0*ds)

  !----------------------------------------------------------------
  ! construct RHS for T-equation
  !----------------------------------------------------------------
  rhs(0:ns-1) = -pressure*theta_s - 2.d0*zeta*theta_s*theta_sss &
              - zeta*theta_ss**2  - mu*(alpha_t+theta_dot)**2 - tau_beam_s &
              -2.d0 * theta_s**2 * zeta_ss &
              -5.d0 * theta_s * theta_ss * zeta_s &
              +mu_star * (theta_s*theta_ss*zeta + tau_beam + zeta_s*theta_s**2)

  ! Dirichlet (overwrite)
  rhs(ns-1) = 0.d0
  ! Neumann (eqn B.4 of PHD thesis)
  rhs(-1)   = (K1 - zeta(0)*theta_ss(0)*theta_s(0) - zeta_s(0)*theta_s(0)**2);

  ! solve system LU decomposition of matrix
  call Solve_linear_system ( A_matrix, rhs, T)


  call Differentiate1D ( T, T_s, ns+1, ds, 1)

end subroutine Tension


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

  ! beam(0:ns-1,1) = beam_solid%x(0:ns-1)
  ! beam(0:ns-1,2) = beam_solid%y(0:ns-1)
  ! beam(0:ns-1,3) = beam_solid%vx(0:ns-1)
  ! beam(0:ns-1,4) = beam_solid%vy(0:ns-1)
  ! beam(0:ns-1,5) = beam_solid%theta(0:ns-1)
  ! beam(0:ns-1,6) = beam_solid%theta_dot(0:ns-1)
  !
  ! call RK4( time, dt, beam, beam_solid%pressure_old(0:ns-1), T, beam_solid%tau_old(0:ns-1), &
  !           beam_solid )
  !
  ! beam_solid%x(0:ns-1)= beam(0:ns-1,1)
  ! beam_solid%y(0:ns-1) = beam(0:ns-1,2)
  ! beam_solid%vx(0:ns-1) = beam(0:ns-1,3)
  ! beam_solid%vy(0:ns-1) = beam(0:ns-1,4)
  ! beam_solid%theta(0:ns-1) = beam(0:ns-1,5)
  ! beam_solid%theta_dot(0:ns-1) = beam(0:ns-1,6)
  !
  ! call integrate_position( time, beam_solid )
end subroutine RK4_wrapper




!-------------------------------------------------------------------------------
subroutine RK4 (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid)
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
 !  beam_old=beam
 !  theta = beam(:,5)
 !  theta_dot = beam(:,6)
 !
 !  theta_1 = 0.d0
 !  theta_2 = 0.d0
 !  theta_3 = 0.d0
 !  theta_4 = 0.d0
 !  theta_dot_1 = 0.d0
 !  theta_dot_2 = 0.d0
 !  theta_dot_3 = 0.d0
 !  theta_dot_4 = 0.d0
 !
 !  T1 = T
 !  T2 = T
 !  T3 = T
 !  T4 = T
 !
 !
 ! !is just a runge-kutta 4th order
 ! !subroutine RHS_beameqn (time, theta , theta_dot, pressure_beam )
 !    theta_1 = theta
 !    theta_dot_1 = theta_dot
 !  call RHS_beameqn (time, theta_1 , theta_dot_1, pressure_beam, T1, tau_beam, beam_solid)
 !    theta_2     = theta + 0.5d0*dt_beam*theta_1
 !    theta_dot_2 = theta_dot + 0.5d0*dt_beam*theta_dot_1
 !  call RHS_beameqn (time+0.5d0*dt_beam, theta_2 , theta_dot_2, pressure_beam, T2, tau_beam, beam_solid)
 !    theta_3     = theta + 0.5d0*dt_beam*theta_2
 !    theta_dot_3 = theta_dot + 0.5d0*dt_beam*theta_dot_2
 !  call RHS_beameqn (time+0.5d0*dt_beam, theta_3 , theta_dot_3, pressure_beam, T3, tau_beam, beam_solid)
 !    theta_4     = theta + dt_beam*theta_3
 !    theta_dot_4 = theta_dot + dt_beam*theta_dot_3
 !  call RHS_beameqn (time+dt_beam, theta_4 , theta_dot_4, pressure_beam, T4, tau_beam, beam_solid)
 !    theta       = theta + dt_beam * (theta_1 + 2.d0*theta_2 + 2.d0*theta_3 + theta_4 )/6.0d0
 !    theta_dot   = theta_dot + dt_beam * (theta_dot_1 + 2.d0*theta_dot_2 + 2.d0*theta_dot_3 + theta_dot_4 )/6.0d0
 !
 !  beam(:,5) = theta
 !  beam(:,6) = theta_dot
 !
 !  T=T1
 !
 !  !---------------------------------------------------------------------
 !  !     emergency brake
 !  !---------------------------------------------------------------------
 !  if (maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0 ) then
 !    write (*,'(A)') "!!! rk4-Solver: I found maxval(abs(beam(:,6)-beam_old(:,6) ))>100.d0"
 !    write (*,'(A)') "possible instability"
 !    write (*,'("time=",es11.4, " dt=",es11.4)') time, dt_beam
 !  endif

end subroutine RK4


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

  ! beam(:,1) = beam_solid%x(0:ns-1)
  ! beam(:,2) = beam_solid%y(0:ns-1)
  ! beam(:,3) = beam_solid%vx(0:ns-1)
  ! beam(:,4) = beam_solid%vy(0:ns-1)
  ! beam(:,5) = beam_solid%theta(0:ns-1)
  ! beam(:,6) = beam_solid%theta_dot(0:ns-1)
  !
  ! call EE1( time, dt, beam, beam_solid%pressure_old(0:ns-1), T, &
  !      beam_solid%tau_old(0:ns-1), beam_solid )
  !
  ! beam_solid%x(0:ns-1) = beam(:,1)
  ! beam_solid%y(0:ns-1) = beam(:,2)
  ! beam_solid%vx(0:ns-1) = beam(:,3)
  ! beam_solid%vy(0:ns-1) = beam(:,4)
  ! beam_solid%theta(0:ns-1) = beam(:,5)
  ! beam_solid%theta_dot(0:ns-1) = beam(:,6)
  !
  ! call integrate_position( time, beam_solid )
end subroutine EE1_wrapper


!-------------------------------------------------------------------------------

subroutine EE1 (time, dt_beam, beam, pressure_beam, T, tau_beam, beam_solid)
  implicit none
  real(kind=pr), intent (in) :: dt_beam, time
  real(kind=pr),dimension(0:ns-1, 1:6), intent (inout) :: beam
  real(kind=pr),dimension(0:ns-1), intent(out) :: T
  real(kind=pr),dimension(0:ns-1), intent (in) :: pressure_beam, tau_beam
  type (solid) :: beam_solid
  real(kind=pr),dimension(0:ns-1) :: theta, theta_dot
  real(kind=pr),dimension(0:ns-1) :: theta_1
  real(kind=pr),dimension(0:ns-1) :: theta_dot_1

  ! theta     = beam(:,5)
  ! theta_dot = beam(:,6)
  !
  ! theta_1     = theta
  ! theta_dot_1 = theta_dot
  !
  ! call RHS_beameqn (time, theta_1 , theta_dot_1, pressure_beam, T, tau_beam, beam_solid)
  !
  ! beam(:,5) = theta + dt_beam * theta_1
  ! beam(:,6) = theta_dot + dt_beam * theta_dot_1

end subroutine EE1




!-------------------------------------------------------------------------------
!  NUMERIC COMPUTATION OF THE JACOBIAN
!  -> I think the row/column ordering in J (in IBES) is ackward; pay attention.
!-------------------------------------------------------------------------------
subroutine Jacobi_num (time, dt, J, x, beam_solid)
  implicit none
  type (solid), intent(in) :: beam_solid
  real(kind=pr), dimension(1:2*ns+4,1:2*ns+4), intent (out) :: J
  real(kind=pr), dimension(1:2*ns+4), intent(in) :: x
  real(kind=pr), dimension(1:2*ns+4) :: F1,F2, x1,x2
  real(kind=pr), intent (in) :: time, dt
  integer :: i

  do i=1,2*ns+4
    x1 = x
    x2 = x
    x1(i) = x1(i) - 1.0d-6
    x2(i) = x2(i) + 1.0d-6

    call F_nonlinear (time, dt, F1, x1, beam_solid)
    call F_nonlinear (time, dt, F2, x2, beam_solid)

    J(i,:)= (F2-F1)/(2.0d-6)
  enddo

end subroutine Jacobi_num


!-------------------------------------------------------------------------------
! initial guess
! return a full vector x with all values for theta, T and all ghost nodes
! to start the newton iteration.
! in the current version, and unlike previous versions, we prefer "guessing"
! the old solution at time t^n. we just have to compute the ghost nodes here.
! INPUT
!   beam: the full beam at time t. we take theta and compute T, and ghost nodes
! OUTPUT
!   x: a full guess [theta(-1) theta(0) ... theta(ns+1) T(-1) ... T(ns-1)] for
!      a state at the new time level n+1
!-------------------------------------------------------------------------------
subroutine initial_guess( time, dt, beam, x_guess )
  use vars
  implicit none
  real(kind=pr),intent (in)   :: dt, time
  type (solid), intent(inout) :: beam
  real(kind=pr),dimension(1:2*ns+4),intent(out) :: x_guess
  real(kind=pr),dimension(-1:ns-1) :: T, T_s
  real(kind=pr),dimension(-1:ns+1) :: theta_guess
  real(kind=pr),dimension(0:ns-1) :: theta
  real(kind=pr) :: K1, K2
  real(kind=pr) :: alpha, alpha_t, alpha_tt
  real(kind=pr), dimension(1:6) :: LeadingEdge ! LeadingEdge: x, y, vx, vy, ax, ay (Array)


  call mouvement(time+dt, alpha, alpha_t, alpha_tt, LeadingEdge, beam )

  ! constants from the boundary conditions
  K1 = beam%mu(0)*( LeadingEdge(5)*dcos(alpha)+LeadingEdge(6)*dsin(alpha)+grav*dsin(alpha)) - beam%tau_old(0)
  K2 = beam%mu(0)*(-LeadingEdge(5)*dsin(alpha)+LeadingEdge(6)*dcos(alpha)+grav*dcos(alpha)) + beam%pressure_old(0)


  ! get tension including ghost node T(-1)
  call Tension ( time, T, T_s, beam%theta, beam%theta_dot, beam%pressure_old, beam%tau_old, beam)

  ! copy interior values of theta
  theta_guess(0:ns-1) = beam%theta(0:ns-1)

  ! theta ghost nodes
  theta = beam%theta(0:ns-1)

  ! ghost node leading edge
  theta_guess(-1) = (10.d0*beam%zeta(0)*theta(0) - 12.d0*beam%zeta(0)*theta(1) + 6.d0*beam%zeta(0)*theta(2) &
                  - beam%zeta(0)*theta(3) + 2.d0*K2*ds**3 - 8.d0*ds*beam%zeta_s(0)*theta(0) + 4.d0*ds*beam%zeta_s(0)*theta(1) &
                  - T(0)*ds**2*theta(1) + ds**2*beam%zeta_ss(0)*theta(1))/(3.d0*beam%zeta(0) - 4.d0*ds*beam%zeta_s(0) &
                  - T(0)*ds**2 + ds**2*beam%zeta_ss(0))

  ! ghost nodes trailing edge
  theta_guess(ns) = (15.d0*beam%zeta(ns-1)*theta(ns-1) - 12.d0*beam%zeta(ns-1)*theta(ns-2) &
                  + beam%zeta(ns-1)*theta(ns-3) + 3.d0*ds*beam%zeta_s(ns-1)*theta(ns-2))&
                  / (4.d0*beam%zeta(ns-1) + 3.d0*ds*beam%zeta_s(ns-1))

  theta_guess(ns+1) = (120.d0*beam%zeta(ns-1)*theta(ns-1) - 128.d0*beam%zeta(ns-1)*theta(ns-2) &
                    + 12.d0*beam%zeta(ns-1)*theta(ns-3) + 3.d0*ds*beam%zeta_s(ns-1)*theta(ns-3))&
                    / (4.d0*beam%zeta(ns-1) + 3.d0*ds*beam%zeta_s(ns-1));

  !-- rearrange the guess in the vector of unknowns x
  x_guess(1:ns+3)      = theta_guess(-1:ns+1)
  x_guess(ns+4:2*ns+4) = T(-1:ns-1)
end subroutine initial_guess

end module solid_model
