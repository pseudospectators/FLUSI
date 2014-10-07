! Wrapper for init_fields
subroutine init_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
 
  if (mpirank==0) write(*,*) "Set up initial conditions...."
  
  !-----------------------------------------------------------------------------
  ! initialize fields, possibly read in backup file
  !----------------------------------------------------------------------------- 
  select case(method)
  case("fsi") 
     call init_fields_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  case("mhd")
     call abort()
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in init_fields."
     call abort()
  end select

  n0 = 1-n1 !important to do this now in case we're retaking a backp
  
  !-----------------------------------------------------------------------------
  ! create startup mask function
  !-----------------------------------------------------------------------------
!   if (mpirank==0) write(*,'("Creating startup mask...time=",es12.4)') time
!   call create_mask(time,Insect,beams)
  
  !-----------------------------------------------------------------------------
  ! save initial conditions (if not resuming a backup)
  !-----------------------------------------------------------------------------
  if (index(inicond,'backup::')==0) then
    if (mpirank==0) write(*,*) "Saving initial conditions to disk..."
    call save_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  endif
end subroutine init_fields