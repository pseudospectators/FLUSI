!-------------------------------------------------------------------------------
! ./flusi -p --magnitude ux_00.h5 uy_00.h5 uz_00.h5 outfile_00.h5
!-------------------------------------------------------------------------------
! load the vector components from file and compute & save the magnitude to
! another HDF5 file.
subroutine magnitude_post(help)
  use vars
  use basic_operators
  use p3dfft_wrapper
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfile
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --magnitude ux_00.h5 uy_00.h5 uz_00.h5 outfile_00.h5 "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! load the vector components from file and compute & save the magnitude to"
    write(*,*) "! another HDF5 file. mag(u) = sqrt(ux^2+uy^2+uz^2)  "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,outfile)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if (mpirank==0) then
    write(*,*) "Computing magnitude of vector from these files: "
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Outfile="//trim(adjustl(outfile))
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  work = dsqrt( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
  deallocate( u )

  call save_field_hdf5 ( time, outfile, work )
  if (mpirank==0) write(*,*) "Wrote magnitude to "//trim(outfile)

  deallocate (work)
  call fft_free()

end subroutine magnitude_post
