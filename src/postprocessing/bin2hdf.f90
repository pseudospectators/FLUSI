!-------------------------------------------------------------------------------
! ./flusi --postprocess --bin2hdf [file_bin] [file_hdf5] [nx] [ny] [nz] [xl] [yl] [zl] [time]
! ./flusi --postprocess --bin2hdf ux_file.binary ux_00000.h5 128 128 384 3.5 2.5 10.0 0.0
!-------------------------------------------------------------------------------
! converts the given binary file into an HDF5 file following flusi's conventions
subroutine convert_bin2hdf(help)
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_bin,fname_hdf,tmp
  ! we read a single precision field from binary
  real, dimension(:,:,:), allocatable :: field
  ! and convert it to a double precision field (however: hdf file is single again)
  real(kind=pr), dimension(:,:,:), allocatable :: field_out
  integer, parameter :: pr_out = 4
  integer :: ix, iy ,iz, i,j,k
  integer(kind=8) :: record_length
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --bin2hdf [file_bin] [file_hdf5] [nx] [ny] [nz] [xl] [yl] [zl] [time]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "converts the given binary file into an HDF5 file following flusi's conventions"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: Nope"
    return
  endif


  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return
  endif


  ! binary file name
  call get_command_argument(3,fname_bin)
  ! hdf5 file name
  call get_command_argument(4,fname_hdf)
  call get_command_argument(5,tmp)
  read (tmp,*) nx
  call get_command_argument(6,tmp)
  read (tmp,*) ny
  call get_command_argument(7,tmp)
  read (tmp,*) nz
  call get_command_argument(8,tmp)
  read (tmp,*) xl
  call get_command_argument(9,tmp)
  read (tmp,*) yl
  call get_command_argument(10,tmp)
  read (tmp,*) zl
  call get_command_argument(11,tmp)
  read (tmp,*) time

  write(*,'("converting ",A," into ",A," resolution: ",3(i4,1x)," box size: ",&
  &3(es15.8,1x)," time=",es15.8)') trim(adjustl(fname_bin)), &
  trim(adjustl(fname_hdf)), nx,ny,nz, xl,yl,zl,time

  !-----------------------------------------------------------------------------
  ! read in the binary field to be converted
  !-----------------------------------------------------------------------------
  allocate ( field(0:nx-1,0:ny-1,0:nz-1), field_out(0:nx-1,0:ny-1,0:nz-1) )

  !   inquire (iolength=record_length) field
  !   open(11, file=fname_bin, form='unformatted', &
  !   access='direct', recl=record_length, convert="little_endian")
  !   read (11,rec=1) field
  !   close (11)
  !   write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)


  OPEN(10,FILE=fname_bin,FORM='unformatted',STATUS='OLD', convert='LITTLE_ENDIAN')
  read(10) (((field(i,j,k),i=0,nx-1),j=0,ny-1),k=0,nz-1)
  CLOSE(10)
  write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)
  !-----------------------------------------------------------------------------
  ! write the field data to an HDF file
  !-----------------------------------------------------------------------------
  ! initializes serial domain decomposition:
  ra=(/0, 0, 0/)
  rb=(/nx-1, ny-1, nz-1/)
  field_out = real(field,kind=pr)
  striding = 1
  call save_field_hdf5(time,fname_hdf,field_out)

  deallocate (field,field_out)
end subroutine convert_bin2hdf
