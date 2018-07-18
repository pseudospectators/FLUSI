!-------------------------------------------------------------------------------
! ./flusi -p --TKE-mean ekinavg_00.h5 uavgx_00.h5 uavgy_00.h5 uavgz_00.h5 tkeavg_000.h5
! From the time-avg kinetic energy field and the components of the time avg
! velocity field, compute the time averaged turbulent kinetic energy.
! See TKE note 18 feb 2015 (Thomas) and 13 feb 2015 (Dmitry)
! Can be done in parallel
!-------------------------------------------------------------------------------
subroutine TKE_mean(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers
  
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, fname_ekin, outfile, mode
  real(kind=pr),dimension(:,:,: ),allocatable :: ekin
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --TKE-mean ekinavg_00.h5 uavgx_00.h5 uavgy_00.h5 uavgz_00.h5 tkeavg_000.h5 --urms"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! From the time-avg kinetic energy field and the components of the time avg"
    write(*,*) "! velocity field, compute the time averaged turbulent kinetic energy."
    write(*,*) "! See TKE note 18 feb 2015 (Thomas) and 13 feb 2015 (Dmitry)"
    write(*,*) "! TKE_avg = ekin_avg - 0.5d0*(ux_avg^2 + uy_avg^2 + uz_avg^2)"
    write(*,*) "! If the --urms switich is set, we save the RMS velocity instead"
    write(*,*) "! URMS = sqrt(2*TKE/3)"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,fname_ux)
  call get_command_argument(5,fname_uy)
  call get_command_argument(6,fname_uz)
  call get_command_argument(7,outfile)
  call get_command_argument(8,mode)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )
  call check_file_exists( fname_ekin )

  if (mpirank == 0) then
    write(*,*) "Processing "//trim(adjustl(fname_ux))//" "//trim(adjustl(fname_uy))//&
    &" "//trim(adjustl(fname_uz))//" and "//fname_ekin
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  call decomposition_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(ekin(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )
  call read_single_file ( fname_ekin, ekin )

  ekin = ekin - 0.5d0*(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)
  if ( mode == "--urms" ) then
      if (root) write(*,*) "computing rms velocity instead of TKE"
      ekin = dsqrt( 2.0d0*ekin/3.0d0)
  endif

  if (mpirank==0) write(*,*) "Wrote to "//trim(adjustl(outfile))
  call save_field_hdf5 ( time,outfile,ekin)


  deallocate (u,ekin)
end subroutine tke_mean
