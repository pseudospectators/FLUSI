!-------------------------------------------------------------------------------
! ./flusi --postprocess --time-avg file_list.txt avgx_0000.h5
! Reads in a list of files from a file, then loads one file after the other and
! computes the average field, which is then stored in the specified file.
!-------------------------------------------------------------------------------
subroutine time_avg_HDF5(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_this, fname_avg
  real(kind=pr), dimension(:,:,:), allocatable :: field_avg, field
  integer :: ix, iy ,iz, io_error=0, i=0
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --time-avg file_list.txt avgx_0000.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Reads in a list of files from a file, then loads one file after the other and"
    write(*,*) " computes the average field, which is then stored in the specified file."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: Yes"
    return
  endif

  call get_command_argument(3,fname)
  call get_command_argument(4,fname_avg)

  !-----------------------------------------------------------------------------
  ! check if input file exists, the file contains the list of h5 files to be avg
  !-----------------------------------------------------------------------------
  call check_file_exists ( fname )
  if (root) write(*,*) "Reading list of files from "//fname

  !-----------------------------------------------------------------------------
  ! read in the file, loop over lines
  !-----------------------------------------------------------------------------
  open( unit=14, file=fname, action='read', status='old')
  do while (io_error==0)
    ! fetch current filename
    read (14,'(A)', iostat=io_error) fname_this

    if (io_error == 0) then
      call check_file_exists ( fname_this )

      ! initialization is done after first read.
      if (i==0) then
        ! get file size etc
        call fetch_attributes( fname_this, nx, ny, nz, xl, yl, zl, time, nu, origin )
        ! initialization parallel module (no FFTS)
        call decomposition_initialize()
        ! memory
        allocate(field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
        allocate(field_avg(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

        field_avg = 0.d0
      endif

      ! read the field from file
      call read_single_file( fname_this, field )

      ! add it to the avg field
      field_avg = field_avg + field

      i = i+1
    endif
  enddo
  close (14)

  ! we're done an can free this array.
  deallocate(field)

  ! compute average
  field_avg = field_avg / dble(i)

  ! save and exit
  call save_field_hdf5(0.d0, fname_avg, field_avg)

  deallocate(field_avg)
end subroutine time_avg_HDF5
