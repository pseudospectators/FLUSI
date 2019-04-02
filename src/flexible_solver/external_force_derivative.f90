!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct Jacobian matrix of internal force vector for Newton-Raphson method
! since we use implicit scheme for time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_forces_derivatives_construction(Wing)

type(flexible_wing), intent(inout)  :: Wing
integer :: np,j

! Get the number of mass points for the sake of simplicity in coding
np = wing%np

! Initialize
wing%FJ_ext = 0.d0

  !Construct external force derivative matrix aka Jacobian matrix
  do j=1,np
      wing%FJ_ext(j,j)           = wing%FJ_ext(j,j) - &
                                   wing%m(j)*(wing%vr0(2)**2 + wing%vr0(3)**2)
      wing%FJ_ext(j,j+np)        = wing%FJ_ext(j,j+np) - &
                                   wing%m(j)*( - wing%vr0(1)*wing%vr0(2) + wing%ar0(3))
      wing%FJ_ext(j,j+2*np)      = wing%FJ_ext(j,j+2*np) + &
                                   wing%m(j)*(wing%vr0(1)*wing%vr0(3) + wing%ar0(2))

      wing%FJ_ext(j+np,j)        = wing%FJ_ext(j+np,j) + &
                                   wing%m(j)*(wing%vr0(1)*wing%vr0(2) + wing%ar0(3))
      wing%FJ_ext(j+np,j+np)     = wing%FJ_ext(j+np,j+np) - &
                                   wing%m(j)*(wing%vr0(1)**2 + wing%vr0(3)**2)
      wing%FJ_ext(j+np,j+2*np)   = wing%FJ_ext(j+np,j+2*np) - &
                                   wing%m(j)*( - wing%vr0(2)*wing%vr0(3) + wing%ar0(1))

      wing%FJ_ext(j+2*np,j)      = wing%FJ_ext(j+2*np,j) - &
                                   wing%m(j)*( - wing%vr0(1)*wing%vr0(3) + wing%ar0(2))
      wing%FJ_ext(j+2*np,j+np)   = wing%FJ_ext(j+2*np,j+np) + &
                                   wing%m(j)*(wing%vr0(2)*wing%vr0(3) + wing%ar0(1))
      wing%FJ_ext(j+2*np,j+2*np) = wing%FJ_ext(j+2*np,j+2*np) - &
                                   wing%m(j)*(wing%vr0(1)**2 + wing%vr0(2)**2)
  enddo

end subroutine
