! checks if a given file ("fname") exists. if not, code is stopped brutally
subroutine check_file_exists(fname)
  use vars
  use mpi
  implicit none

  character (len=*), intent(in) :: fname
  logical :: exist1
  integer :: mpicode

  inquire ( file=fname, exist=exist1 )
  if ( exist1 .eqv. .false.) then
    if (root) write (*,'("ERROR! file: ",A," not found")') trim(adjustl(fname))
    call abort(223,"file not found")
  endif

end subroutine check_file_exists



! overwrite and initialize file
subroutine init_empty_file( fname )
  use vars
  implicit none
  character (len=*), intent(in) :: fname

  if (root) then
    open (15, file=fname, status='replace')
    close(15)
  endif
  
end subroutine
