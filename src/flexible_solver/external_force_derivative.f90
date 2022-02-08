!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct Jacobian matrix of internal force vector for Newton-Raphson method
! since we use implicit scheme for time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_forces_derivatives_construction(FJ,Wing,dt1)

type(flexible_wing), intent(in)  :: Wing
real(kind=pr), intent(inout) :: FJ(1:,1:)
real(kind=pr), intent(in) :: dt1
real(kind=pr) :: t0
integer :: np,j

! Get the number of mass points for the sake of simplicity in coding
np = wing%np

! Initialize
!t0 = MPI_wtime()
!wing%FJ = 0.d0
!call toc("Flexible wing (FJ_initialize)", MPI_wtime() - t0)

  !Construct external force derivative matrix aka Jacobian matrix
  do j=1,np
      FJ(j,j)           = FJ(j,j) - &
                              wing%coef*dt1*wing%m(j)*(wing%vr0(2)**2 + wing%vr0(3)**2)
      FJ(j,j+np)        = FJ(j,j+np) - &
                              wing%coef*dt1*wing%m(j)*( - wing%vr0(1)*wing%vr0(2) + wing%ar0(3))
      FJ(j,j+2*np)      = FJ(j,j+2*np) + &
                              wing%coef*dt1*wing%m(j)*(wing%vr0(1)*wing%vr0(3) + wing%ar0(2))

      FJ(j+np,j)        = FJ(j+np,j) + &
                                   wing%coef*dt1*wing%m(j)*(wing%vr0(1)*wing%vr0(2) + wing%ar0(3))
      FJ(j+np,j+np)     = FJ(j+np,j+np) - &
                                   wing%coef*dt1*wing%m(j)*(wing%vr0(1)**2 + wing%vr0(3)**2)
      FJ(j+np,j+2*np)   = FJ(j+np,j+2*np) - &
                                   wing%coef*dt1*wing%m(j)*( - wing%vr0(2)*wing%vr0(3) + wing%ar0(1))

      FJ(j+2*np,j)      = FJ(j+2*np,j) - &
                                   wing%coef*dt1*wing%m(j)*( - wing%vr0(1)*wing%vr0(3) + wing%ar0(2))
      FJ(j+2*np,j+np)   = FJ(j+2*np,j+np) + &
                                   wing%coef*dt1*wing%m(j)*(wing%vr0(2)*wing%vr0(3) + wing%ar0(1))
      FJ(j+2*np,j+2*np) = FJ(j+2*np,j+2*np) - &
                                   wing%coef*dt1*wing%m(j)*(wing%vr0(1)**2 + wing%vr0(2)**2)
  enddo

end subroutine
