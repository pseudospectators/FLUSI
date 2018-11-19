!-------------------------------------------------------------------------------
! Flexible solid time stepping routines
! BDF2 method with Euler startup
!-------------------------------------------------------------------------------

subroutine flexible_solid_time_step(time, dt0, dt1, it, wings)
    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(flexible_wing), dimension(1:nWings), intent (inout) :: wings
    real(kind=pr) :: c1, c2, c3
    integer :: i

    ! select scheme
    if (it == 0) then

        call translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)

        ! Construct the external force vector
        write(*,*) wings(1)%at_inertia
        call external_forces_construction(time,dt0,dt1, it,wings)
        ! EULER startup scheme
        ! compute position and velocity at new time step
        ! (updates wings%u_new using wings%u_old)
        call flexible_solid_solver_euler(time, dt1, it, wings)
        !call moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,Wings)
    else

        call translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)

        write(*,*) wings(1)%at_inertia
        ! Construct the external force vector
        call external_forces_construction(time,dt0,dt1, it,wings)

        ! BDF2 scheme
        ! compute position and velocity at new time step
        ! (updates wings%u_new using wings%u_old and wings%u_old)
        call flexible_solid_solver_BDF2(time, dt0, dt1, it, wings)
        !call moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,Wings)

    endif

    do i=1,nWings
      !-----------------------------------------------------------------------------
      ! emergency brake
      !-----------------------------------------------------------------------------
      if ((maxval(abs(wings(i)%u_new(1:6*wings(i)%np) - &
          wings(i)%u_old(1:6*wings(i)%np)))>100.d0 ).and.(root)) then
        write (*,'(A)') "!!! Flexible-Solid-Solver: "
        write(*,'("Running on ",i5," CPUs")') mpisize
        write(*,'("I found maxval(abs(wings(",i5,")%u_new - wings(",i5,")%u_old))>100.d0")') i, i
        write (*,'(A)') "    That indicates a possible instability."
        write (*,'("time=",es11.4)') time
      endif
    enddo

    ! TODO Adding periodization: if the insect moves out of the box for free flight

end subroutine

subroutine flexible_solid_solver_euler(time, dt1, it, wings)

    implicit none

    real(kind=pr),intent(in) :: time, dt1
    integer,intent (in) :: it
    type(flexible_wing), dimension(1:nWings), intent (inout) :: wings
    real(kind=pr) :: du, err, err_rel, coef=1.0
    integer :: i, iter, i_NAN, j_NAN, iJ,jJ,np
    logical :: iterate

  do i=1,nWings

    iterate = .true.
    err     = 1.0d0
    err_rel = 1.0d0
    iter    = 0

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wings(i)%u_new = wings(i)%u_old

    ! Begin the Newton-Raphson iterations
    do while (iterate)

      iter = iter + 1

      ! Get total number of mass points
      np = wings(i)%np

      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wings(i))

      if (Vector_isNAN(wings(i)%Fint(1:3*wings(i)%np))) then
        if (root) write(*,*) "FlexibleSolidSolver: Internal force vector contains NaNs"
        call abort(9836,"The internal force vector for the solid solver contains NaN..abort")
      endif

      ! Calculate RHS vector
      call RHS_for_NR_method(dt1, dt1, it, wings(i))

      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(wings(i))

      if (Matrix_isNAN(wings(i)%FJ(1:3*wings(i)%np,1:3*wings(i)%np),i_NAN, j_NAN)) then
        if (root) write(*,*) "FlexibleSolidSolver: Jacobian matrix contains NaNs"
        call abort(9835,"The Jacobian matrix for the solid solver contains NaN..abort")
      endif

      ! Solve for the step of NR method
      call solve_linear_system_using_schur_complement(wings(i)%du(1:6*np), np, &
                                                      wings(i)%FJ(1:3*np,1:3*np), &
                                                      wings(i)%m(1:np), wings(i)%c(1:np), &
                                                      dt1, wings(i)%RHS_a(1:3*np), wings(i)%RHS_b(1:3*np), coef)


      wings(i)%u_new(1:6*wings(i)%np) = wings(i)%u_new(1:6*wings(i)%np) - wings(i)%du(1:6*wings(i)%np)
      err            = dsqrt(sum(wings(i)%du**2))
      err_rel        = abs(dsqrt(sum(wings(i)%du**2)) / dsqrt(sum(wings(i)%u_new**2)))


      ! convergence test
      if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
        iterate = .false.
      endif

      ! emergency brake
      if (iter>20) then
        call abort(309214,"!!! ERROR: The flexible solid solver performed like 50 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wings(i)%u_oldold = wings(i)%u_old
    wings(i)%u_old = wings(i)%u_new
    wings(i)%x(1:wings(i)%np)  = wings(i)%u_old(1:wings(i)%np)
    wings(i)%y(1:wings(i)%np)  = wings(i)%u_old(wings(i)%np+1:2*wings(i)%np)
    wings(i)%z(1:wings(i)%np)  = wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np)
    wings(i)%vx(1:wings(i)%np) = wings(i)%u_old(3*wings(i)%np+1:4*wings(i)%np)
    wings(i)%vy(1:wings(i)%np) = wings(i)%u_old(4*wings(i)%np+1:5*wings(i)%np)
    wings(i)%vz(1:wings(i)%np) = wings(i)%u_old(5*wings(i)%np+1:6*wings(i)%np)
  enddo

end subroutine

subroutine flexible_solid_solver_BDF2(time, dt1, dt0, it, wings)

    implicit none

    real(kind=pr),intent(in) :: time, dt1, dt0
    integer,intent (in) :: it
    type(flexible_wing), dimension(1:nWings), intent (inout) :: wings
    real(kind=pr) :: r, coef
    real(kind=pr) :: du, err, err_rel
    integer :: i, iter, np
    logical :: iterate

    ! Calculate the coefficient for time stepping scheme
    r = dt1/dt0
    coef = (1+r)/(1+2*r)

  do i=1,nWings

    ! Get total number of mass points
    np = wings(i)%np

    iterate = .true.
    err     = 1.0d0
    err_rel = 1.0d0
    iter    = 0

    ! Assign the initial guess value for the Newton-Raphson method equal to the
    ! value of the state vector of the previous time step
    wings(i)%u_new = wings(i)%u_old

    do while (iterate)

      iter = iter + 1


      ! Calculate internal force vector from the new state vector u_new
      call internal_forces_construction(wings(i))

      ! Calculate RHS vector
      call RHS_for_NR_method(dt1, dt0, it, wings(i))

      ! Calculate the Jacobian matrix of the internal force vector
      call internal_forces_derivatives_construction(wings(i))

      ! Solve for the step of NR method
      call solve_linear_system_using_schur_complement(wings(i)%du(1:6*np), np, &
                                                      wings(i)%FJ(1:3*np,1:3*np), &
                                                      wings(i)%m(1:np), wings(i)%c(1:np), &
                                                      dt1, wings(i)%RHS_a(1:3*np), wings(i)%RHS_b(1:3*np), coef)

      wings(i)%u_new = wings(i)%u_new - wings(i)%du
      err            = dsqrt(sum(wings(i)%du**2))
      err_rel        = abs(dsqrt(sum(wings(i)%du**2)) / dsqrt(sum(wings(i)%u_new**2)))

      ! convergence test
      if ( (((err<error_stop) .or. (err_rel<error_stop)).and.(iter>2))) then
        iterate = .false.
        if (root) then
          write(*,*) minval((/err, err_rel/))
          write(*,*) iter
        endif
      endif


      ! emergency brake
      if (iter>100) then
        call abort(337214,"!!! ERROR: The flexible solid solver performed like 50 iterations. This is not normal!")
      endif

    enddo

    ! Update results
    wings(i)%u_oldold = wings(i)%u_old
    wings(i)%u_old = wings(i)%u_new
    wings(i)%x(1:wings(i)%np)  = wings(i)%u_old(1:wings(i)%np)
    wings(i)%y(1:wings(i)%np)  = wings(i)%u_old(wings(i)%np+1:2*wings(i)%np)
    wings(i)%z(1:wings(i)%np)  = wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np)
    wings(i)%vx(1:wings(i)%np) = wings(i)%u_old(3*wings(i)%np+1:4*wings(i)%np)
    wings(i)%vy(1:wings(i)%np) = wings(i)%u_old(4*wings(i)%np+1:5*wings(i)%np)
    wings(i)%vz(1:wings(i)%np) = wings(i)%u_old(5*wings(i)%np+1:6*wings(i)%np)

  enddo

end subroutine

subroutine RHS_for_NR_method(dt1, dt0, it, wings)

! Calculate the right hand side (RHS) F for the Newton Raphson method which is
!                         J*dx = F
! where J is the Jacobian matrix and dx is the step needed for searching
! the root of the equations

  implicit none
  integer, intent(in) :: it
  real(kind=pr),intent(in) :: dt1, dt0
  type(flexible_wing), intent (inout) :: wings
  real(kind=pr) :: c1, c2, c3, r
  integer :: np

  !Get total number of points for the simplification of writing array bounds
  np = wings%np


    if (it==0) then  ! this is the first time step

        wings%RHS_a(1:3*np) = (/wings%m(1:np), wings%m(1:np), &
                      wings%m(1:np)/)*(wings%u_new(3*np+1:6*np) - wings%u_old(3*np+1:6*np)) - &
                         dt1*(wings%Fext(1:3*np) + wings%Fint(1:3*np) - &
                         (/wings%c(1:np), wings%c(1:np), wings%c(1:np)/)*(wings%u_new(3*np+1:6*np) + &
                         wings%u_old(3*np+1:6*np)))

        wings%RHS_b(1:3*np) = wings%u_new(1:3*np) - wings%u_old(1:3*np) - &
                      (wings%u_new(3*np+1:6*np))*dt1

    else

      !Calculate scheme coefficients
      r = dt1/dt0
      c1=(1+r)**2/(1+2*r)
      c2=r**2/(1+2*r)
      c3=(1+r)/(1+2*r)

      wings%RHS_a(1:3*np) = (/wings%m(1:np), wings%m(1:np), &
                  wings%m(1:np)/)*(wings%u_new(3*np+1:6*np) - c1*wings%u_old(3*np+1:6*np) + &
                       c2*wings%u_oldold(3*np+1:6*np)) - &
                       c3*dt1*(wings%Fext(1:3*np) + wings%Fint(1:3*np) - &
                       (/wings%c(1:np), wings%c(1:np), &
                       wings%c(1:np)/)*wings%u_new(3*np+1:6*np))

    wings%RHS_b(1:3*np) = wings%u_new(1:3*np) - c1*wings%u_old(1:3*np) + &
                       c2*wings%u_oldold(1:3*np) - c3*(wings%u_new(3*np+1:6*np))*dt1

    endif


end subroutine

subroutine solve_linear_system_using_schur_complement(du,np,FJ,m,c,dt,a,b,coef)

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
real(kind=pr), intent(in) :: FJ(1:,1:)
integer, intent(in) :: np
real(kind=pr), intent(inout) :: du(1:)
real(kind=pr), allocatable :: y(:), m_array3D(:), c_array3D(:)
real(kind=pr), allocatable :: F(:,:)
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
        F(i,i) = coef*dt*FJ(i,i) + (1/(coef*dt))*m_array3D(i) + c_array3D(i)
          y(i) = a(i) + ((1/(coef*dt))*m_array3D(i) + c_array3D(i))*b(i)
      else
        F(i,j) = coef*dt*FJ(i,j)
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

! y_(n+1) - C1*y_(n) + C2*y_(n-1) - C3*dt*f(t_n+1,y_n+1) = 0
!
! e = dt1/dt0
! c1=(1+e)**2/(1+2*e); c2=e**2/(1+2*e); c3=(1+e)/(1+2*e);
