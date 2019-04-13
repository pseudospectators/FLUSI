
subroutine cal_objectiveFunctional_test(ux,J,u0)
  use mpi 

  implicit none
  real(kind=pr),intent(in)  :: ux(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(in)  :: u0(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(out) :: J 

  !local variables
  integer                   :: ierr

  J= sum(ux - u0)**2._pr 

  call MPI_ALLREDUCE(J,J,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  

end subroutine cal_objectiveFunctional_test


subroutine add_adjointSourceTerm_test(rhsx)
  use vars_adjoint

  implicit none
  real(kind=pr),intent(inout):: rhsx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  !local variables
  logical,save               :: alreadyCalled

  if (.not.alreadyCalled) then
    rhsx = rhsx + 2._pr*u_forward   
    alreadyCalled = .true.
  endif

end subroutine add_adjointSourceTerm_test


subroutine cal_gradient_test(u_adj,gradient,i,j,k,ivar)
  implicit none

  real(kind=pr),intent(in)                             :: u_adj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),dimension(:),allocatable,intent(inout) :: gradient
  integer,      intent(in)                             :: i,j,k,ivar

  !local variables
  integer                                              :: ierr

  if (.not.allocated(gradient)) allocate(gradient(1))

  if (i.ge.ra(1).and.i.le.rb(1).and.j.ge.ra(2).and.j.le.rb(2).and.k.ge.ra(3).and.k.le.rb(3)) then
    gradient = u_adj(i,j,k,ivar)
  else
    gradient = 0._pr
  endif

  call MPI_ALLREDUCE(gradient(1),gradient(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

end subroutine cal_gradient_test



subroutine add_forwardSourceTerm_test(rhsx,eps,i,j,k,ivar)  
  implicit none

  real(kind=pr),intent(inout):: rhsx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(in)   :: eps
  integer,      intent(in)   :: i,j,k,ivar
   
  integer                    :: ii

  if (i.ge.ra(1).and.i.le.rb(1).and.j.ge.ra(2).and.j.le.rb(2).and.k.ge.ra(3).and.k.le.rb(3)) then
    rhsx (i,j,k,ivar) = rhsx (i,j,k,ivar) + eps
  endif


end subroutine add_forwardSourceTerm_test
