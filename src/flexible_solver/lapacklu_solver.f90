module lapacklu_solver

  use vars
  implicit none

contains

subroutine solve_linear_system_using_schur_complement(du,np,FJ,FJ_ext,m,c,dt,a,b,coef,SparseSolver)


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
character, intent(in) :: SparseSolver
integer, intent(in) :: np
real(kind=pr), intent(inout) :: du(1:)
real(kind=pr), allocatable :: y(:), m_array3D(:), c_array3D(:)
real(kind=pr), allocatable :: F(:,:)
real(kind=pr) :: t0
integer :: i, j, iopt, nnz, info, nrhs_lu, nnz_FJ, nnz_FJext, nnz_pre
integer :: factors(8)
integer, allocatable :: rowind(:), colptr(:)
real(kind=pr), allocatable :: values(:)


allocate(y(1:3*np),m_array3D(1:3*np),c_array3D(1:3*np))
allocate(F(1:3*np,1:3*np))

!Initialize F
F = 0.0d0
y = 0.0d0

! Transform the mass array m() into an array (/m,m,m/) corresponding to three dimensions
! x,y and z
m_array3D = (/m,m,m/)
c_array3D = (/c,c,c/)

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


t0 = MPI_wtime()
 F = coef*dt*(FJ+FJ_ext)
call toc("Flexible wing (Construct global jacobian matrix:the remainings)", MPI_wtime() - t0)

t0 = MPI_wtime()
do i=1,3*np
  F(i,i) = F(i,i) + (1/(coef*dt))*m_array3D(i) + c_array3D(i)
    y(i) = a(i) + ((1/(coef*dt))*m_array3D(i) + c_array3D(i))*b(i)
enddo
call toc("Flexible wing (Construct global jacobian:diagonal)", MPI_wtime() - t0)


  if (root) then
    write(*,*) "Warning! You are solving the dynamics of the mass-spring system without &
                SuperLU solver. This can take forever to run the simulation. It is recommanded &
                to model flexible wing with SuperLU solver!"
  endif


  t0 = MPI_wtime()
  call solve_linear_system_wing ( F, y, du(1:3*np) )
  call toc("Flexible wing (solve linear system using lu factorization)", MPI_wtime() - t0)




do i=1,3*np

  du(i + 3*np) = (1/(coef*dt))*(du(i) - b(i))

enddo

deallocate(y)
deallocate(F)

end subroutine


end module lapacklu_solver
