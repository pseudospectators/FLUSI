!-------------------------------------------------------------------------------
! Flexible solid time stepping routines
! BDF2 method with Euler startup
!-------------------------------------------------------------------------------

subroutine flexible_solid_time_step(time, dt0, dt1, it, wing)
    use vars
    use mpi
    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr) :: c1, c2, c3
    real(kind=pr) :: t0
    integer :: i,itri,j,np

    !-----------------------------------------------------------------------------
    ! Startup time step: use EE1 as first step of BDF2
    !-----------------------------------------------------------------------------
    if (wing%StartupStep) then  ! this is the first time step
      wing%StartupStep = .false.  ! we're about to do the first step
      if (TimeMethodFlexibleSolid=="BDF2") then  ! if we deal with BDF2
        ActuallyBDF2 = .true.    ! Remember to switch back to BDF2 (at the end of the step)
      endif
    endif

    ! select scheme
    if (TimeMethodFlexibleSolid=='EI1' .or. ((ActuallyBDF2) .and. (TimeMethodFlexibleSolid=="BDF2"))) then

        !call translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)
        call flexible_wing_motions ( time, wing )

        ! EULER startup scheme
        ! compute position and velocity at the first time step
        call flexible_solid_solver_euler(time, dt1, it, wing)
        !call moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,Wings)

        !-----------------------------------------------------------------------------
        ! Change marching scheme, if necessary
        !-----------------------------------------------------------------------------
        if (TimeMethodFlexibleSolid=="BDF2") then   ! Remember to switch back to BDF2
          ActuallyBDF2 = .false.
        endif

    elseif (TimeMethodFlexibleSolid=='BDF2') then

        !call translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)
        call flexible_wing_motions ( time, wing )


        ! BDF2 scheme
        ! compute position and velocity at new time step
        ! (updates wings%u_new using wings%u_old and wings%u_old)
        call flexible_solid_solver_BDF2(time, dt0, dt1, it, wing)
        !call moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,Wings)

    endif

      !-----------------------------------------------------------------------------
      ! emergency brake
      !-----------------------------------------------------------------------------
      if ((maxval(abs(wing%u_new(1:6*wing%np) - &
          wing%u_old(1:6*wing%np)))>100.d0 ).and.(root)) then
        write (*,'(A)') "!!! Flexible-Solid-Solver: "
        write(*,'("Running on ",i5," CPUs")') mpisize
        write(*,'("I found maxval(abs(wings(",i5,")%u_new - wings(",i5,")%u_old))>100.d0")') i, i
        write (*,'(A)') "    That indicates a possible instability."
        write (*,'("time=",es11.4)') time
      endif

    ! TODO Adding periodization: if the insect moves out of the box for free flight

end subroutine

subroutine flexible_solid_solver_euler(time, dt1, it, wing)

    use vars

    implicit none

    real(kind=pr),intent(in) :: time, dt1
    integer,intent (in) :: it
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr) :: du, err, err_rel, coef=1.0
    integer :: i, iter, i_NAN, j_NAN, iJ, jJ, np
    logical :: iterate
    real(kind=pr) :: t0

    iterate = .true.
    err     = 1.0d0
    err_rel = 1.0d0
    iter    = 0

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wing%u_new = wing%u_old

    ! Begin the Newton-Raphson iterations
    do while (iterate)

      iter = iter + 1

      ! Get total number of mass points
      np = wing%np

      t0 = MPI_wtime()
      ! Construct the external force vector
      call external_forces_construction(time,it,wing)
      call toc("Flexible wing (external_forces_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wing)
      call toc("Flexible wing (internal_forces_construction)", MPI_wtime() - t0)

      if (Vector_isNAN(wing%Fint(1:3*wing%np))) then
        if (root) write(*,*) "FlexibleSolidSolver: Internal force vector contains NaNs"
        call abort(9836,"The internal force vector for the solid solver contains NaN..abort")
      endif

      ! Calculate RHS vector
      call RHS_for_NR_method(dt1, dt1, it, wing)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the external force vector
      call external_forces_derivatives_construction(wing)
      call toc("Flexible wing (external_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(wing)
      call toc("Flexible wing (internal_forces_derivatives_construction)", MPI_wtime() - t0)

      if (Matrix_isNAN(wing%FJ(1:3*np,1:3*np),i_NAN, j_NAN)) then
        if (root) write(*,*) "FlexibleSolidSolver: Jacobian matrix contains NaNs"
        call abort(9835,"The Jacobian matrix for the solid solver contains NaN..abort")
      endif

      t0 = MPI_wtime()
      ! Solve for the step of NR method
      call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                      wing%FJ(1:3*np,1:3*np),wing%FJ_ext(1:3*np,1:3*np), &
                                                      wing%m(1:np), wing%c(1:np), &
                                                      dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), coef)
      call toc("Flexible wing (solve_linear_system_using_schur_complement)", MPI_wtime() - t0)

      wing%u_new(1:6*wing%np) = wing%u_new(1:6*wing%np) - wing%du(1:6*wing%np)
      err            = dsqrt(sum(wing%du**2))
      err_rel        = abs(dsqrt(sum(wing%du**2)) / dsqrt(sum(wing%u_new**2)))


      ! convergence test
      if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
        iterate = .false.
        
      endif

      ! emergency brake
      if (iter>100) then
        call abort(309214,"!!! ERROR: The flexible solid solver performed like 100 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wing%u_oldold = wing%u_old
    wing%u_old = wing%u_new

    call update_solver_solutions_to_create_wing_mask(wing,time)

end subroutine

subroutine flexible_solid_solver_BDF2(time, dt1, dt0, it, wing)

    use vars

    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr) :: r, coef
    real(kind=pr) :: du, err, err_rel
    integer :: i, iter, np
    logical :: iterate
    real(kind=pr) :: t0

    ! Calculate the coefficient for time stepping scheme
    r = dt1/dt0
    coef = (1+r)/(1+2*r)

    ! Get total number of mass points
    np = wing%np

    iterate = .true.
    err     = 1.0d0
    err_rel = 1.0d0
    iter    = 0

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wing%u_new = wing%u_old

    do while (iterate)

      iter = iter + 1

      t0 = MPI_wtime()
      ! Construct the external force vector
      call external_forces_construction(time,it,wing)
      call toc("Flexible wing (external_forces_construction)", MPI_wtime() - t0)


      t0 = MPI_wtime()
      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wing)
      call toc("Flexible wing (internal_forces_construction)", MPI_wtime() - t0)

      ! Calculate RHS vector
      call RHS_for_NR_method(dt1, dt0, it, wing)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the external force vector
      call external_forces_derivatives_construction(wing)
      call toc("Flexible wing (external_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(wing)
      call toc("Flexible wing (internal_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Solve for the step of NR method
      call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                      wing%FJ(1:3*np,1:3*np),wing%FJ_ext(1:3*np,1:3*np), &
                                                      wing%m(1:np), wing%c(1:np), &
                                                      dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), coef)
      call toc("Flexible wing (solve_linear_system_using_schur_complement)", MPI_wtime() - t0)

      wing%u_new = wing%u_new - wing%du
      err            = dsqrt(sum(wing%du**2))
      err_rel        = abs(dsqrt(sum(wing%du**2)) / dsqrt(sum(wing%u_new**2)))

      ! convergence test
      if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
        iterate = .false.
      
      endif

      ! emergency brake
      if (iter>100) then
        call abort(337214,"!!! ERROR: The flexible solid solver performed like 100 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wing%u_oldold = wing%u_old
    wing%u_old = wing%u_new

    call update_solver_solutions_to_create_wing_mask(wing,time)

end subroutine

subroutine RHS_for_NR_method(dt1, dt0, it, wing)

! Calculate the right hand side (RHS) F for the Newton Raphson method which is
!                         J*dx = F
! where J is the Jacobian matrix and dx is the step needed for searching
! the root of the equations

  implicit none
  integer, intent(in) :: it
  real(kind=pr),intent(in) :: dt1, dt0
  type(flexible_wing), intent (inout) :: wing
  real(kind=pr) :: c1, c2, c3, r
  integer :: np

  !Get total number of points for the simplification of writing array bounds
  np = wing%np


    if (it==0) then  ! this is the first time step

        wing%RHS_a(1:3*np) = (/wing%m(1:np), wing%m(1:np), &
                      wing%m(1:np)/)*(wing%u_new(3*np+1:6*np) - wing%u_old(3*np+1:6*np)) - &
                         dt1*(wing%Fext(1:3*np) + wing%Fint(1:3*np) - &
                         (/wing%c(1:np), wing%c(1:np), wing%c(1:np)/)*(wing%u_new(3*np+1:6*np) + &
                         wing%u_old(3*np+1:6*np)))

        wing%RHS_b(1:3*np) = wing%u_new(1:3*np) - wing%u_old(1:3*np) - &
                      (wing%u_new(3*np+1:6*np))*dt1

    else

      !Calculate scheme coefficients
      r = dt1/dt0
      c1=(1+r)**2/(1+2*r)
      c2=r**2/(1+2*r)
      c3=(1+r)/(1+2*r)

      wing%RHS_a(1:3*np) = (/wing%m(1:np), wing%m(1:np), &
                  wing%m(1:np)/)*(wing%u_new(3*np+1:6*np) - c1*wing%u_old(3*np+1:6*np) + &
                       c2*wing%u_oldold(3*np+1:6*np)) - &
                       c3*dt1*(wing%Fext(1:3*np) + wing%Fint(1:3*np) - &
                       (/wing%c(1:np), wing%c(1:np), &
                       wing%c(1:np)/)*wing%u_new(3*np+1:6*np))

    wing%RHS_b(1:3*np) = wing%u_new(1:3*np) - c1*wing%u_old(1:3*np) + &
                       c2*wing%u_oldold(1:3*np) - c3*(wing%u_new(3*np+1:6*np))*dt1

    endif


end subroutine

subroutine solve_linear_system_using_schur_complement(du,np,FJ,FJ_ext,m,c,dt,a,b,coef)

! This subroutine solves a linear system under the form
!                 [    |      ] [    ]   [   ]
!                 [ FJ |   M  ] [ dx ]   [ a ]
!                 [----|------] [----] = [---]
!                 [  I |  c.I ] [ dv ]   [ b ]
!                 [    |      ] [    ]   [   ]
! where A is an arbitrary matrix deriving from Jacobian matrix getting from
! taking derivatives of all the internal forces, M is the mass and damping
! matrix who MUST be DIAGONAL where M = diag(m) + coef*dt*diag(c),
! I is IDENTICAL matrix, D is IDENTICAL matrix times a CONSTANT COEFFICIENT coef
! determined by the discretization scheme
! ATTENTION to the structure of matrices M, I and D
! a and b are the RHS getting from evaluating the state of the system from
! previous time step.
! x and v and position and velocity vector respectively, and the phase vector is
! defined as u = [x ; v]^T
!
! (FJ-M*c**(-1))*dx = a - M*c**(-1)*b
! du = 1/(c*dt)*(dx - b)

implicit none

real(kind=pr), intent(in) :: dt, coef
real(kind=pr), intent(in) :: m(1:), c(1:), a(1:), b(1:)
real(kind=pr), intent(in) :: FJ(1:,1:), FJ_ext(1:,1:)
integer, intent(in) :: np
real(kind=pr), intent(inout) :: du(1:)
real(kind=pr), allocatable :: y(:), m_array3D(:), c_array3D(:)
real(kind=pr), allocatable :: F(:,:)
real(kind=pr) :: t0
integer :: i, j

allocate(y(1:3*np),m_array3D(1:3*np),c_array3D(1:3*np))
allocate(F(1:3*np,1:3*np))

! Transform the mass array m() into an array (/m,m,m/) corresponding to three dimensions
! x,y and z
m_array3D = (/m,m,m/)
c_array3D = (/c,c,c/)

do i=1,3*np
    do j=1,3*np
      if (i .eq. j) then
        F(i,i) = coef*dt*(FJ(i,i)+FJ_ext(i,i)) + (1/(coef*dt))*m_array3D(i) + c_array3D(i)
          y(i) = a(i) + ((1/(coef*dt))*m_array3D(i) + c_array3D(i))*b(i)
      else
        F(i,j) = coef*dt*(FJ(i,j)+FJ_ext(i,i))
      endif
    enddo
enddo

call solve_linear_system_wing ( F, y, du(1:3*np) )

do i=1,3*np

  du(i + 3*np) = (1/(coef*dt))*(du(i) - b(i))

enddo

deallocate(y)
deallocate(F)

end subroutine

subroutine update_solver_solutions_to_create_wing_mask(wing,time)

implicit none

real(kind=pr), intent(in) :: time
type(flexible_wing), intent (inout) :: wing
real(kind=pr), dimension(1:3) :: u, v, w
integer :: j
real(kind=pr) :: tau, vr

if ((wing%Motion .eq. "stationary") .or. (wing%Motion .eq. "simple_harmonic") &
    .or. (wing%Motion .eq. "harmonic_ocsillation")) then
  wing%x(1:wing%np)  = wing%u_old(1:wing%np)
  wing%y(1:wing%np)  = wing%u_old(wing%np+1:2*wing%np)
  wing%z(1:wing%np)  = wing%u_old(2*wing%np+1:3*wing%np)
  wing%vx(1:wing%np) = wing%u_old(3*wing%np+1:4*wing%np)
  wing%vy(1:wing%np) = wing%u_old(4*wing%np+1:5*wing%np)
  wing%vz(1:wing%np) = wing%u_old(5*wing%np+1:6*wing%np)
elseif (wing%Motion .eq. "revolving_wing") then

  call rotate_wing(wing)

  call translate_wing(wing)

endif


end subroutine


! y_(n+1) - C1*y_(n) + C2*y_(n-1) - C3*dt*f(t_n+1,y_n+1) = 0
!
! e = dt1/dt0
! c1=(1+e)**2/(1+2*e); c2=e**2/(1+2*e); c3=(1+e)/(1+2*e);
