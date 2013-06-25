! Wrapper for init_fields
subroutine init_fields(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  use mpi_header
  use vars
  implicit none

  integer,intent (inout) :: n1,it
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex (kind=pr), intent(inout) :: &
       uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout) :: &
       work_nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr), intent(inout) ::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)

  if (mpirank == 0) write(*,*) "Initializating fields:"
 
  select case(method)
  case("fsi") 
     call init_fields_fsi(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  case("mhd")
     call init_fields_mhd(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in init_fields."
     call abort
  end select
end subroutine init_fields
