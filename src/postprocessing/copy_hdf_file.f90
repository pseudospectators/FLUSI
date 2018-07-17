!-------------------------------------------------------------------------------
! copy one hdf5 file to another one, with different name. Not you cannot do this
! in terminal with "cp" since "cp" does not touch the dataset name in the file,
! which is then not conformal to flusi naming convention.
! call:
! ./flusi --postprocess --cp ux_00000.h5 new_00000.h5
!
!
! since I learned about h5copy tool, this subroutine is deprecated
!
! h5copy -i mask_000000.h5 -s mask -o hallo.h5 -d test
!
!-------------------------------------------------------------------------------
subroutine copy_hdf_file(help)
  use mpi
  use vars
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_in, fname_out
  real(kind=pr)::time
  ! input field
  real(kind=pr), dimension(:,:,:), allocatable :: field_in

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "/flusi --postprocess --cp ux_00000.h5 new_00000.h5 "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! copy one hdf5 file to another one, with different name. Not you cannot do this"
    write(*,*) "! in terminal with cp since cp does not touch the dataset name in the file,"
    write(*,*) "! which is then not conformal to flusi naming convention."
    write(*,*) "!"
    write(*,*) "! since I learned about h5copy tool, this subroutine is deprecated"
    write(*,*) "!"
    write(*,*) "! h5copy -i mask_000000.h5 -s mask -o hallo.h5 -d test"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: no"
    return
  endif


  if (mpisize/=1) then
    write(*,*)
    call abort(111, "./flusi --postprocess --cp is a SERIAL routine, use 1CPU only")
  endif



  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  ! get filename to save file to
  call get_command_argument(4,fname_out)

  call fetch_attributes( fname_in, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ra=0
  rb=(/nx-1,ny-1,nz-1/)
  write(*,*) "copying ",trim(adjustl(fname_in)), " to ", trim(adjustl(fname_out))

  allocate ( field_in(0:nx-1,0:ny-1,0:nz-1) )
  call read_single_file(fname_in,field_in)

  call save_field_hdf5 ( time, fname_out, field_in )

  deallocate (field_in)
end subroutine copy_hdf_file
