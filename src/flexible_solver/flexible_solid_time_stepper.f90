!-------------------------------------------------------------------------------
! Flexible solid time stepping routines
! BDF2 method with Euler startup
!-------------------------------------------------------------------------------

subroutine flexible_solid_time_step(time, dt0, dt1, it, wing, Insect)
    use vars
    use mpi
    use module_insects
    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(diptera), intent(inout) :: Insect
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

    !-----------------------------------------------------------------------------
    ! fetch current motion state
    !-----------------------------------------------------------------------------
    call BodyMotion (time, Insect)
    call StrokePlane (time, Insect)
    if (wing%ID == "left") then
      call Flexible_wing_motions ( time, wing, Insect%kine_wing_l )
    elseif (wing%ID == "right") then
      call Flexible_wing_motions ( time, wing, Insect%kine_wing_r )
    endif

    !-----------------------------------------------------------------------------
    ! define the rotation matrices to change between coordinate systems
    !-----------------------------------------------------------------------------
    call body_rotation_matrix( Insect, Insect%M_body )
    Insect%M_body_inv = transpose(Insect%M_body)
    call MSM_solver_rotation_matrix( Wing, wing%M_solver )
    Wing%M_solver_inv = transpose(Wing%M_solver)
    call flexible_wing_rotation_matrix( Wing, Insect, Wing%M_wing )
    Wing%M_wing_inv = transpose(Wing%M_wing)

    ! rel+abs wing angular velocities in the w/b/g coordinate system
    call flexible_wing_angular_velocities (time, Wing, Insect, Insect%M_body )
    call flexible_wing_angular_accel( time, wing, Insect )

    ! select scheme
    if (TimeMethodFlexibleSolid=='EI1' .or. ((ActuallyBDF2) .and. (TimeMethodFlexibleSolid=="BDF2"))) then

        ! EULER startup scheme
        ! compute position and velocity at the first time step
        call flexible_solid_solver_euler(time, dt1, it, wing, Insect)

        !-----------------------------------------------------------------------------
        ! Change marching scheme, if necessary
        !-----------------------------------------------------------------------------
        if (TimeMethodFlexibleSolid=="BDF2") then   ! Remember to switch back to BDF2
          ActuallyBDF2 = .false.
        endif

    elseif (TimeMethodFlexibleSolid=='BDF2') then


        ! BDF2 scheme
        ! compute position and velocity at new time step using BDF2 scheme
        call flexible_solid_solver_BDF2(time, dt0, dt1, it, wing, Insect)
    elseif (TimeMethodFlexibleSolid=='Verlet') then


          ! Verlet explicit scheme
          ! compute position and velocity at new time step using Verlet scheme
          call flexible_solid_solver_Verlet(time, dt1, it, wing, Insect)

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

subroutine flexible_solid_solver_euler(time, dt1, it, wing, Insect)

    use vars

    implicit none

    real(kind=pr),intent(in) :: time, dt1
    integer,intent (in) :: it
    type(diptera), intent(inout) :: Insect
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr), allocatable :: FJ(:,:)
    real(kind=pr) :: du, coef=1.0
    integer :: i, i_NAN, j_NAN, iJ, jJ, np
    logical :: iterate
    real(kind=pr) :: t0

    iterate      = .true.
    wing%err_abs = 1.0d0
    wing%err_rel = 1.0d0
    wing%iter    = 0
    wing%coef    = 1.0

    ! Get total number of mass points
    np = wing%np

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wing%u_new = wing%u_old

    if (.not.allocated(FJ)) allocate(FJ(1:3*np,1:3*np))

    ! Begin the Newton-Raphson iterations
    do while (iterate)

      wing%iter = wing%iter + 1

      t0 = MPI_wtime()

      ! Construct the external force vector
      call external_forces_construction(time,it,wing,Insect)
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
      FJ = 0.0
      call toc("Flexible wing (FJ_initialize)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the external force vector
      call external_forces_derivatives_construction(FJ,wing,dt1)
      call toc("Flexible wing (external_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(FJ,wing,dt1)
      call toc("Flexible wing (internal_forces_derivatives_construction)", MPI_wtime() - t0)

      if (Matrix_isNAN(FJ(1:3*np,1:3*np),i_NAN, j_NAN)) then
        if (root) write(*,*) "FlexibleSolidSolver: Jacobian matrix contains NaNs"
        call abort(9835,"The Jacobian matrix for the solid solver contains NaN..abort")
      endif

      t0 = MPI_wtime()
      ! Solve for the step of NR method
      if (wing%HB_matrix_given) then
        call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                        FJ(1:3*np,1:3*np), &
                                                        wing%m(1:np), wing%c(1:np), &
                                                        dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), &
                                                        wing%coef,  wing%SparseSolver, it, wing%ID, &
                                                        wing%colptr,wing%rowind)
      else
        call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                        FJ(1:3*np,1:3*np), &
                                                        wing%m(1:np), wing%c(1:np), &
                                                        dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), &
                                                        wing%coef,  wing%SparseSolver, it, wing%ID)
      endif
      call toc("Flexible wing (solve_linear_system_using_schur_complement)", MPI_wtime() - t0)

      wing%u_new(1:6*wing%np) = wing%u_new(1:6*wing%np) - wing%du(1:6*wing%np)
      wing%err_abs        = dsqrt(sum(wing%du**2))
      wing%err_rel        = abs(dsqrt(sum(wing%du**2)) / dsqrt(sum(wing%u_new**2)))


      ! convergence test
      if ( (((wing%err_abs<error_stop) .or. (wing%err_rel<error_stop)).and.(wing%iter>2))) then
        iterate = .false.

      endif

      ! emergency brake
      if (wing%iter>100) then
        call abort(309214,"!!! ERROR: The flexible solid solver performed like 200 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wing%u_oldold = wing%u_old
    wing%u_old = wing%u_new

    call update_solver_solutions_to_create_wing_mask(time,wing,Insect)

    if (allocated(FJ)) deallocate(FJ)

end subroutine

subroutine flexible_solid_solver_BDF2(time, dt1, dt0, it, wing, Insect)

    use vars

    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(diptera), intent(inout) :: Insect
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr), allocatable :: FJ(:,:)
    real(kind=pr) :: r
    real(kind=pr) :: du
    integer :: i, np, nnz_BDF2
    logical :: iterate
    real(kind=pr) :: t0

    ! Calculate the coefficient for time stepping scheme
    r = dt1/dt0
    wing%coef = (1+r)/(1+2*r)

    ! Get total number of mass points
    np = wing%np

    iterate = .true.
    wing%err_abs     = 1.0d0
    wing%err_rel = 1.0d0
    wing%iter    = 0

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wing%u_new = wing%u_old

    if (.not.allocated(FJ)) allocate(FJ(1:3*np,1:3*np))

    do while (iterate)

      wing%iter = wing%iter + 1

      t0 = MPI_wtime()
      ! Construct the external force vector
      call external_forces_construction(time,it,wing,Insect)
      call toc("Flexible wing (external_forces_construction)", MPI_wtime() - t0)


      t0 = MPI_wtime()
      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wing)
      call toc("Flexible wing (internal_forces_construction)", MPI_wtime() - t0)

      ! Calculate RHS vector
      call RHS_for_NR_method(dt1, dt0, it, wing)

      t0 = MPI_wtime()
      FJ = 0.0
      call toc("Flexible wing (FJ_initialize)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the external force vector
      call external_forces_derivatives_construction(FJ,wing,dt1)
      call toc("Flexible wing (external_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(FJ,wing,dt1)
      call toc("Flexible wing (internal_forces_derivatives_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Solve for the step of NR method
      if (wing%HB_matrix_given) then
        call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                      FJ(1:3*np,1:3*np), &
                                                      wing%m(1:np), wing%c(1:np), &
                                                      dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), &
                                                      wing%coef, wing%SparseSolver, it, wing%ID, &
                                                      wing%colptr, wing%rowind)
      else
        call solve_linear_system_using_schur_complement(wing%du(1:6*np), np, &
                                                      FJ(1:3*np,1:3*np), &
                                                      wing%m(1:np), wing%c(1:np), &
                                                      dt1, wing%RHS_a(1:3*np), wing%RHS_b(1:3*np), &
                                                      wing%coef, wing%SparseSolver, it, wing%ID)
      endif
      call toc("Flexible wing (solve_linear_system_using_schur_complement)", MPI_wtime() - t0)

      wing%u_new = wing%u_new - wing%du
      wing%err_abs           = dsqrt(sum(wing%du**2))
      wing%err_rel        = abs(dsqrt(sum(wing%du**2)) / dsqrt(sum(wing%u_new**2)))

      ! convergence test
      if ( (((wing%err_abs<error_stop) .or. (wing%err_rel<error_stop)).and.(wing%iter>2))) then
        iterate = .false.

      endif

      ! emergency brake
      if (wing%iter>200) then
        call abort(337214,"!!! ERROR: The flexible solid solver performed like 200 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wing%u_oldold = wing%u_old
    wing%u_old = wing%u_new

    call update_solver_solutions_to_create_wing_mask(time,wing,Insect)

    if (allocated(FJ)) deallocate(FJ)

end subroutine

subroutine flexible_solid_solver_verlet(time, dt1, it, wing, Insect)

    use vars

    implicit none

    real(kind=pr),intent(in) :: time, dt1
    integer,intent (in) :: it
    type(diptera), intent(inout) :: Insect
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr) :: du, coef=1.0
    real(kind=pr), allocatable :: point_forces(:)
    integer :: i, i_NAN, j_NAN, iJ, jJ, nop
    logical :: iterate
    real(kind=pr) :: t0


    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wing%u_new = wing%u_old

      ! Get total number of mass points
      nop = wing%np

      allocate(point_forces(1:3*nop))

      ! Update position from previous velovity and acceleration
      do i=1,nop

        wing%u_new(i) = wing%u_old(i) + wing%u_old(i+3*nop)*dt1 + wing%ax_old(i)*0.5*dt1**2
        wing%u_new(i+nop) = wing%u_old(i+nop) + wing%u_old(i+4*nop)*dt1 + wing%ay_old(i)*0.5*dt1**2
        wing%u_new(i+2*nop) = wing%u_old(i+2*nop) + wing%u_old(i+5*nop)*dt1 + wing%az_old(i)*0.5*dt1**2

      enddo

      ! Update forces from the new position and velocity

      ! Construct the external force vector
      t0 = MPI_wtime()
      call external_forces_construction(time,it,wing,Insect)
      call toc("Flexible wing (external_forces_construction)", MPI_wtime() - t0)

      t0 = MPI_wtime()
      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wing)
      call toc("Flexible wing (internal_forces_construction)", MPI_wtime() - t0)

      if (Vector_isNAN(wing%Fint(1:3*wing%np))) then
        if (root) write(*,*) "FlexibleSolidSolver: Internal force vector contains NaNs"
        call abort(9836,"The internal force vector for the solid solver contains NaN..abort")
      endif

      !Update acceleration from the new forces
        point_forces(1:3*nop) = (wing%Fext(1:3*nop) + wing%Fint(1:3*nop) - &
                                (/wing%c(1:nop), wing%c(1:nop), wing%c(1:nop)/)*wing%u_new(3*nop+1:6*nop))

      do i=1,nop

          wing%ax_new(i) = point_forces(i)/wing%m(i)
          wing%ay_new(i) = point_forces(i+nop)/wing%m(i)
          wing%az_new(i) = point_forces(i+2*nop)/wing%m(i)

      enddo

      !Update the velocity from the new acceleration

      do i=1,nop

          wing%u_new(i+3*nop) = wing%u_old(i+3*nop) + 0.5*dt1*(wing%ax_old(i) + wing%ax_new(i))
          wing%u_new(i+4*nop) = wing%u_old(i+4*nop) + 0.5*dt1*(wing%ay_old(i) + wing%ay_new(i))
          wing%u_new(i+5*nop) = wing%u_old(i+5*nop) + 0.5*dt1*(wing%az_old(i) + wing%az_new(i))

      enddo

      ! emergency brake
      if (maxval(wing%u_new)>1e6) then
        call abort(309214,"!!! ERROR: The flexible solid solver returns large value. This is not normal!")
      endif

    ! Update results
    wing%ax_old = wing%ax_new
    wing%ay_old = wing%ay_new
    wing%az_old = wing%az_new
    wing%u_old = wing%u_new

    call update_solver_solutions_to_create_wing_mask(time,wing,Insect)

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

subroutine solve_linear_system_using_schur_complement(du,np,FJ,m,c,dt,a,b,&
  coef,SparseSolver,it,wing_ID,colptr_in,rowind_in)

use vars

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
real(kind=pr), intent(inout) :: FJ(1:,1:)
character, intent(in) :: SparseSolver
integer, intent(in) :: np, it
character(len=strlen), intent(in) :: wing_ID
real(kind=pr), intent(inout) :: du(1:)
integer, intent(in),optional :: rowind_in(1:,1:), colptr_in(1:)
real(kind=pr), allocatable :: y(:)
real(kind=pr) :: t0
integer :: i, j, iopt, nnz, info, nrhs_lu, nnz_FJ, nnz_FJext, nnz_pre
integer :: p, col
integer :: factors(8)
integer, allocatable :: rowind(:), colptr(:)
real(kind=pr), allocatable :: values(:)
character(len=3)  :: wingstr


allocate(y(1:3*np))

!Initialize F
y = 0.0d0

! Transform the mass array m() into an array (/m,m,m/) corresponding to three dimensions
! x,y and z
!m_array3D = (/m,m,m/)
!c_array3D = (/c,c,c/)

!do i=1,3*np
!    do j=1,3*np
!      if (i .eq. j) then
!        F(i,i) = coef*dt*(FJ(i,i)+FJ_ext(i,i)) + (1/(coef*dt))*m_array3D(i) + c_array3D(i)
!          y(i) = a(i) + ((1/(coef*dt))*m_array3D(i) + c_array3D(i))*b(i)
!      else
!        F(i,j) = coef*dt*(FJ(i,j)+FJ_ext(i,j))
!      endif
!    enddo
!enddo

write (wingstr,'(i3)') it


t0 = MPI_wtime()
do i=1,np
  do j=0,2
  FJ(i+j*np,i+j*np) = FJ(i+j*np,i+j*np) + (1/(coef*dt))*m(i) + c(i)
          y(i+j*np) = a(i+j*np) + ((1/(coef*dt))*m(i) + c(i))*b(i+j*np)
  enddo
enddo
call toc("Flexible wing (Construct global jacobian:diagonal)", MPI_wtime() - t0)

t0 = MPI_wtime()
call solve_linear_system_wing ( FJ(1:3*np,1:3*np), y(1:3*np), du(1:3*np) )
call toc("Flexible wing (solve linear system using lu factorization)", MPI_wtime() - t0)

do i=1,3*np
  du(i + 3*np) = (1/(coef*dt))*(du(i) - b(i))
enddo

deallocate(y)


end subroutine

subroutine update_solver_solutions_to_create_wing_mask(time, wing, Insect)

implicit none

real(kind=pr), intent(in) :: time
type(diptera), intent(inout) :: Insect
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
elseif ((wing%Motion .eq. "from_file") .or. (wing%Motion .eq. "revolving_wing")) then

  !call rotate_wing(wing)

  !call translate_wing(wing)

  call rotate_and_translate_wing_into_global_system(wing, Insect)

  call construct_total_velocity(wing,Insect%M_body,Insect%M_body_inv)

endif


end subroutine
