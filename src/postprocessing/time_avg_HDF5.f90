!-------------------------------------------------------------------------------
! ./flusi --postprocess --time-avg file_list.txt avgx_0000.h5
! Reads in a list of files from a file, then loads one file after the other and
! computes the average field, which is then stored in the specified file.
!-------------------------------------------------------------------------------
subroutine time_avg_HDF5(help)
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_bin, fname_avg
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
    write(*,*) "Parallel: NO"
    return
  endif


  call get_command_argument(3,fname)
  call get_command_argument(4,fname_avg)

  !-----------------------------------------------------------------------------
  ! check if input file exists, the file contains the list of h5 files to be avg
  !-----------------------------------------------------------------------------
  call check_file_exists ( fname )
  write(*,*) "Reading list of files from "//fname

  if ( mpisize>1 ) then
    write (*,*) "--time-avg is currently a serial version only, run it on 1CPU"
    return
  endif

  !-----------------------------------------------------------------------------
  ! read in the file, loop over lines
  !-----------------------------------------------------------------------------
  open( unit=14,file=fname, action='read', status='old')
  do while (io_error==0)
    ! fetch current filename
    read (14,'(A)', iostat=io_error) fname_bin
    write(*,*) "read "//trim(adjustl(fname_bin))
    if (io_error == 0) then
      write(*,*) "Processing file "//trim(adjustl(fname_bin))

      call check_file_exists ( fname_bin )
      call fetch_attributes( fname_bin, nx, ny, nz, xl, yl, zl, time, nu, origin )

      ! first time? allocate then.
      if ( .not. allocated(field_avg) ) then
        ra=(/0,0,0/)
        rb=(/nx-1,ny-1,nz-1/)
        allocate(field_avg(0:nx-1,0:ny-1,0:nz-1))
        allocate(field(0:nx-1,0:ny-1,0:nz-1))
        field_avg = 0.d0
      endif

      ! read the field from file
      call read_single_file( fname_bin, field )

      field_avg = field_avg + field

      i = i+1
    endif
  enddo
  close (14)

  field_avg = field_avg / dble(i)

  call save_field_hdf5(0.d0, fname_avg, field_avg)

  deallocate(field_avg)
  deallocate(field)


end subroutine time_avg_HDF5
