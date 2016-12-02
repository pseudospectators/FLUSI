!-------------------------------------------------------------------------------
! ./flusi --postprocess --hdf2bin ux_00000.h5 filename.bin
!-------------------------------------------------------------------------------
! converts the *.h5 file to an ordinairy binary file
subroutine convert_hdf2bin(help)
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_bin
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer, parameter :: pr_out = 4
  integer :: ix, iy ,iz
  real(kind=pr_out), dimension(:,:,:), allocatable :: field_out ! single precision
  real(kind=pr) :: time


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --hdf2bin ux_00000.h5 filename.bin"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "convert the given HDF5 file to a FORTRAN binary file"
    write(*,*) "Ordering:"
    write(*,*) "write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)"
    write(*,*) "SINGLE PRECISION, LITTLE_ENDIAN"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: NO"
    return
  endif


  call get_command_argument(3,fname)
  call get_command_argument(4,fname_bin)

  ! check if input file exists
  call check_file_exists ( fname )

  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return
  endif

  call fetch_attributes( fname, nx, ny, nz, xl, yl, zl, time, nu )

  write (*,'("Converting ",A," to ",A," Resolution is",3(i4,1x))') &
  trim(fname), trim(fname_bin), nx,ny,nz
  write (*,'("time=",es12.4," xl=",es12.4," yl=",es12.4," zl=",es12.4)') &
  time, xl, yl, zl

  ra=(/0,0,0/)
  rb=(/nx-1,ny-1,nz-1/)
  allocate ( field(0:nx-1,0:ny-1,0:nz-1),field_out(0:nx-1,0:ny-1,0:nz-1) )

  ! read field from hdf file
  call read_single_file (fname, field)

  ! convert to single precision
  field_out = real(field, kind=pr_out)

  write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field_out),minval(field_out)

  ! dump binary file (this file will be called ux_00100.h5.binary)
  open (12, file = trim(fname_bin), form='unformatted', status='replace',&
  convert="little_endian")
  write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
  !  write(12) field_out
  close (12)

  deallocate (field, field_out)
end subroutine convert_hdf2bin
