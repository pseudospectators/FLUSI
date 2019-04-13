module module_adjoint
 
  use vars 

  implicit none

  private

  ! public functions
  public :: cal_objectiveFunctional
  public :: add_adjointSourceTerm
  public :: cal_gradient
  public :: add_forwardSourceTerm

contains

  !---------------------------------------
  ! note these include files also have to be specified as dependencies in the
  ! Makefile for make to check if one of them changed
  include "case_test.f90"
  !---------------------------------------

  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the objective functional. 
  !-------------------------------------------------------------------------------
  subroutine cal_objectiveFunctional(ux,J,u0,time)
    implicit none
    real(kind=pr),intent(in)           :: ux(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),intent(out)          :: J
    real(kind=pr),intent(in), optional :: u0(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),intent(in), optional :: time

    select case (adjoint_case)
      case("test")
        if (.not.present(u0)) then
          call abort(1904131314,"When the adjoint case is 'test' the ini condition is necessary &
                     to calculate the objective functional.")
        endif
        call cal_objectiveFunctional_test(ux,J,u0)  
      case default 
        call abort(1904051545,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select

  end subroutine cal_objectiveFunctional

  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the adjoint source term.
  !-------------------------------------------------------------------------------
  subroutine add_adjointSourceTerm(rhsx,time)
    implicit none 
    real(kind=pr),intent(inout)        :: rhsx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),intent(in), optional :: time

    select case (adjoint_case)
      case("test")
        call add_adjointSourceTerm_test(rhsx)  
      case default 
        call abort(1904051546,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select


  end subroutine add_adjointSourceTerm


  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the gradient. 
  !-------------------------------------------------------------------------------
  subroutine cal_gradient(u_adj,gradient,time,i,j,k,ivar)
    implicit none
    
    real(kind=pr),intent(in)                             :: u_adj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),dimension(:),allocatable,intent(inout) :: gradient
    real(kind=pr),intent(in), optional                   :: time
    integer      ,intent(in), optional                   :: i,j,k,ivar
    

    select case (adjoint_case)
      case("test")
        if (.not.present(i).or..not.present(j).or..not.present(k).or..not.present(ivar)) then
          call abort(1904131343,"When the adjoint case is 'test' the position of the source term is necessary &
                     to calculate the adjoint gradient.")
        endif
        call cal_gradient_test(u_adj,gradient,i,j,k,ivar)  
      case default 
        call abort(1904051547,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select

  end subroutine cal_gradient


  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the source term in the forward calculation.
  ! only needed for testing!!! 
  !-------------------------------------------------------------------------------
  subroutine add_forwardSourceTerm(rhsx,eps,i,j,k,ivar)
    implicit none 
    
    real(kind=pr),intent(inout)        :: rhsx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),intent(in)           :: eps
    integer      ,intent(in)           :: i,j,k,ivar

    select case (adjoint_case)
      case("test")
        call cal_forwardSourceTerm_test(rhsx,eps,i,j,k,ivar)  
      case default 
        call abort(1904051548,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select
  end subroutine add_forwardSourceTerm 


end module module_adjoint
