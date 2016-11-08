!-------------------------------------------------------------------------------
! ./flusi --postprocess --vor-abs ux_00000.h5 uy_00000.h5 uz_00000.h5 --second-order
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! directly compute the absolute value of vorticity, do not save components
! can be done in parallel
subroutine convert_abs_vorticity(help)
  use vars
  use p3dfft_wrapper
  use mpi
  use helpers
  use basic_operators
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, fname_vor, order
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --vor-abs ux_00000.h5 uy_00000.h5 uz_00000.h5 vorabs_0000.h5 [--second-order]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " load the velocity components from file and compute & save the vorticity"
    write(*,*) " directly compute the absolute value of vorticity, do not save components"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,fname_vor)
  call get_command_argument(7,order)

  call check_file_exists(fname_ux)
  call check_file_exists(fname_uy)
  call check_file_exists(fname_uz)

  if (mpirank == 0) then
    write(*,*) "Compute magnitude(vorticity) from velocity files:"
    write(*,*) "ux="//trim(adjustl(fname_ux))
    write(*,*) "uy="//trim(adjustl(fname_uy))
    write(*,*) "uz="//trim(adjustl(fname_uz))
    write(*,*) "output="//trim(adjustl(fname_vor))
    write(*,*) "order flag is: "//trim(adjustl(order))
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu )

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  if (mpirank==0) write (*,*) "Done fft_initialize"

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  if (mpirank==0) write (*,*) "Allocated memory"

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft (uk(:,:,:,1),u(:,:,:,1))
  call fft (uk(:,:,:,2),u(:,:,:,2))
  call fft (uk(:,:,:,3),u(:,:,:,3))

  if (order=="--second-order") then
    if (mpirank==0) write(*,*) "using second order!"
    call curl_2nd(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  else
    if (mpirank==0) write(*,*) "using spectral accuracy!"
    call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  endif

  call ifft (u(:,:,:,1),uk(:,:,:,1))
  call ifft (u(:,:,:,2),uk(:,:,:,2))
  call ifft (u(:,:,:,3),uk(:,:,:,3))


  if (mpirank == 0) then
    write (*,'("Writing mag(vor) to file: ",A)') trim(fname_vor)
  endif

  ! compute absolute vorticity:
  u(:,:,:,1) = dsqrt(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)

  call save_field_hdf5 ( time,fname_vor,u(:,:,:,1) )

  deallocate (u)
  deallocate (uk)
  call fft_free()

end subroutine convert_abs_vorticity
