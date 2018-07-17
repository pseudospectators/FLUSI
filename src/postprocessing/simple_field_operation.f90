!-------------------------------------------------------------------------------
! ./flusi -p --simple-field-operation field1.h5 OP field2.h5 output.h5
!-------------------------------------------------------------------------------
! load two fields and perform a simple operation ( + - / * )
! example: ./flusi -p --simple-field-operation vorabs_0000.h5 * mask_0000.h5 vor2_0000.h5
! this gives just the product vor*mask
subroutine simple_field_operation(help)
  use vars
  use basic_operators
  use p3dfft_wrapper
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: file1, file2, file_out, operation
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --simple-field-operation field1.h5 OP field2.h5 output.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! load two fields and perform a simple operation ( + - / * )"
    write(*,*) "! example: ./flusi -p --simple-field-operation vorabs_0000.h5 * mask_0000.h5 vor2_0000.h5"
    write(*,*) "! this gives just the product vor*mask"
    write(*,*) "! NOTE: USE AN ESCAPE CHARACTER FOR OPERATION!"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,file1)
  call get_command_argument(4,operation)
  call get_command_argument(5,file2)
  call get_command_argument(6,file_out)


  call check_file_exists( file1 )
  call check_file_exists( file2 )

  if (mpirank==0) write(*,*) "Performing simple operation using the fields"
  if (mpirank==0) write(*,*) trim(adjustl(file1))//" "//trim(adjustl(file2))
  if (mpirank==0) write(*,*) "Operation="//trim(adjustl(operation))
  if (mpirank==0) write(*,*) "Outfile="//trim(adjustl(file_out))

  call fetch_attributes( file1, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))

  call read_single_file ( file1, u(:,:,:,1) )
  call read_single_file ( file2, u(:,:,:,2) )

  select case (operation)
  case ("*")
    if(mpirank==0) write(*,*) "multiplication"
    u(:,:,:,3) = u(:,:,:,1) * u(:,:,:,2)
  case ("/")
    if(mpirank==0) write(*,*) "division"
    u(:,:,:,3) = u(:,:,:,1) / u(:,:,:,2)
  case ("+")
    if(mpirank==0) write(*,*) "addition"
    u(:,:,:,3) = u(:,:,:,1) + u(:,:,:,2)
  case ("-")
    if(mpirank==0) write(*,*) "substraction"
    u(:,:,:,3) = u(:,:,:,1) - u(:,:,:,2)
  case default
    write(*,*) "error operation not supported::"//operation
  end select

  call save_field_hdf5 ( time, file_out, u(:,:,:,3) )

  if (mpirank==0) write(*,*) "Wrote result to "//trim(file_out)

  deallocate (u)
  call fft_free()

end subroutine simple_field_operation
