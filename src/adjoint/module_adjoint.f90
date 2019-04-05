module module_adjoint
 
  use vars 

  implicit none

  private

  ! public functions
  public :: cal_objectiveFunctional
  public :: cal_adjointSourceTerm
  public :: cal_gradient
  public :: cal_forwardSourceTerm

contains

  !---------------------------------------
  ! note these include files also have to be specified as dependencies in the
  ! Makefile for make to check if one of them changed
  include "case_test.f90"
  !---------------------------------------

  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the objective functional. 
  !-------------------------------------------------------------------------------
  subroutine cal_objectiveFunctional()
    implicit none

    select case (adjoint_case)
      case("test")
        call cal_objectiveFunctional_test()  
      case default 
        call abort(1,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select

  end subroutine cal_objectiveFunctional

  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the adjoint source term.
  !-------------------------------------------------------------------------------
  subroutine cal_adjointSourceTerm()
    implicit none 

    select case (adjoint_case)
      case("test")
        call cal_adjointSourceTerm_test()  
      case default 
        call abort(1,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select


  end subroutine cal_adjointSourceTerm


  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the gradient. 
  !-------------------------------------------------------------------------------
  subroutine cal_gradient()
    implicit none

    select case (adjoint_case)
      case("test")
        call cal_gradient_test()  
      case default 
        call abort(1,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select



  end subroutine cal_gradient


  !-------------------------------------------------------------------------------
  ! Wrapper for the calculation of the source term in the forward calculation.
  ! only needed for testing!!! 
  !-------------------------------------------------------------------------------
  subroutine cal_forwardSourceTerm()
    implicit none 

    select case (adjoint_case)
      case("test")
        call cal_forwardSourceTerm_test()  
      case default 
        call abort(1,"Please set a valid variable for 'case' in the Adjoint section. &
               Possible are: 'test' ") 
    end select
  end subroutine cal_forwardSourceTerm 


end module module_adjoint
