!-------------------------------------------------------------------------------
! ./flusi -p --energy ux_00.h5 uy_00.h5 uz_00.h5 outfile_00.h5
!-------------------------------------------------------------------------------
! load the vector components from file and compute & save the energy to
! another HDF5 file. (energy = (ux^2+uy^2+uz^2)/2)
subroutine energy_post(help)
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
    write(*,*) "./flusi -p --energy ux_00.h5 uy_00.h5 uz_00.h5 outfile_00.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! load the vector components from file and compute & save the energy to"
    write(*,*) "! another HDF5 file. (energy = (ux^2+uy^2+uz^2)/2)"
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

  if (mpirank==0) write(*,*) "Computing energy of vector from these files: "
  if (mpirank==0) write(*,*) trim(adjustl(fname_ux))//" "//trim(adjustl(fname_uy))//" "//trim(adjustl(fname_uz))
  if (mpirank==0) write(*,*) "Outfile="//trim(adjustl(outfile))

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  work = 0.5d0 * ( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )

  call save_field_hdf5 ( time, outfile, work )
  if (mpirank==0) write(*,*) "Wrote energy to "//trim(outfile)

  deallocate (u,work)
  call fft_free()

end subroutine energy_post
