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
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
 
  if (mpirank==0) write(*,*) "Set up initial conditions...."
  
  !-----------------------------------------------------------------------------
  ! initialize fields, possibly read in backup file
  !----------------------------------------------------------------------------- 
  call init_fields_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)

  !-----------------------------------------------------------------------------
  ! create startup mask function
  !-----------------------------------------------------------------------------
  if (mpirank==0) write(*,'("Creating startup mask...time=",es12.4)') time%time
  call create_mask( time%time, mask, mask_color, us, Insect, beams )
  
  !-----------------------------------------------------------------------------
  ! save initial conditions (if not resuming a backup)
  !-----------------------------------------------------------------------------
  if (index(inicond,'backup::')==0) then
    if (mpirank==0) write(*,*) "Saving initial conditions to disk..."
    call save_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  endif
end subroutine init_fields