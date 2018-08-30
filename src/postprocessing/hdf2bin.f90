!-------------------------------------------------------------------------------
! ./flusi --postprocess --hdf2bin ux_00000.h5 filename.bin [--double-precision]
!-------------------------------------------------------------------------------
! converts the *.h5 file to an ordinairy binary file
subroutine convert_hdf2bin(help)
  use vars
  use mpi
  use basic_operators
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_bin, pr_out_flag
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer :: ix, iy ,iz
  real(kind=pr) :: time


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --hdf2bin ux_00000.h5 filename.bin [--double-precision]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "convert the given HDF5 file to a FORTRAN binary file"
    write(*,*) "Ordering:"
    write(*,*) "write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)"
    write(*,*) "LITTLE ENDIAN"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: NO"
    return
  endif


  call get_command_argument(3,fname)
  call get_command_argument(4,fname_bin)
  call get_command_argument(5,pr_out_flag)

  ! check if input file exists
  call check_file_exists ( fname )

  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return
  endif

  call fetch_attributes( fname, nx, ny, nz, xl, yl, zl, time, nu, origin )

  write (*,'("Converting ",A," to ",A," Resolution is",3(i4,1x))') &
  trim(fname), trim(fname_bin), nx,ny,nz
  write (*,'("time=",es12.4," xl=",es12.4," yl=",es12.4," zl=",es12.4)') &
  time, xl, yl, zl

  ra=(/0,0,0/)
  rb=(/nx-1,ny-1,nz-1/)
  allocate ( field(0:nx-1,0:ny-1,0:nz-1) )

  ! read field from hdf file
  call read_single_file (fname, field)

  write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)

  ! dump binary file (this file will be called ux_00100.h5.binary)
  open (12, file = trim(fname_bin), form='unformatted', status='replace',&
  convert="little_endian")

  ! set the output precision and write into a file
  if (pr_out_flag == "--double-precision") then
    write (12) ((( real(field(ix,iy,iz),kind=8), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
    write(*,*) "DOUBLE PRECISION, LITTLE ENDIAN"
  else
    write (12) ((( real(field(ix,iy,iz),kind=4), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
    write(*,*) "SINGLE PRECISION, LITTLE ENDIAN"
  endif

  close (12)

  deallocate (field)
end subroutine convert_hdf2bin
