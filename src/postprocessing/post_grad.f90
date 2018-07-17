!-------------------------------------------------------------------------------
! compute gradient of a scalar field
!
!-------------------------------------------------------------------------------
subroutine post_grad(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_outx, fname_outy, fname_outz, pr_out_flag
  complex(kind=pr),dimension(:,:,:,:),allocatable :: vecfldk
  real(kind=pr),dimension(:,:,:,:),allocatable :: vecfld
  complex(kind=pr),dimension(:,:,:),allocatable :: scfldk
  real(kind=pr),dimension(:,:,:),allocatable :: scfld
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --gradient scfld_000.h5 dscflddx_000.h5 dscflddy_000.h5 dscflddz_000.h5 [--double-precision]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! compute gradient of a scalar field"
    write(*,*) "! with spectral accuracy"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  call get_command_argument(3,fname)
  call get_command_argument(4,fname_outx)
  call get_command_argument(5,fname_outy)
  call get_command_argument(6,fname_outz)
  call get_command_argument(7,pr_out_flag)

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "gradient computation"
    write(*,*) "Computing gradient of a scalar field given in this file: "
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(fname_outx))
    write(*,*) trim(adjustl(fname_outy))
    write(*,*) trim(adjustl(fname_outz))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname )

  call fetch_attributes( fname, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  ! allocate memory
  allocate(scfld(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(scfldk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  ! read scalar field from file
  call read_single_file ( fname, scfld(:,:,:) )

  ! Fourier transform
  call fft( inx=scfld, outk=scfldk )
  deallocate (scfld)

  ! gradient in Fourier space
  allocate(vecfldk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  call gradient(scfldk,vecfldk(:,:,:,1),vecfldk(:,:,:,2),vecfldk(:,:,:,3))
  deallocate (scfldk)

  ! inverse Fourier transform
  allocate(vecfld(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  call ifft3( ink=vecfldk, outx=vecfld )
  deallocate (vecfldk)

  ! set the output precision and write into a file
  if (pr_out_flag == "--double-precision") then
    field_precision = "double"
    if (mpirank==0) then
      write(*,*) "DOUBLE PRECISION"
    endif
  else
    field_precision = "single"
    if (mpirank==0) then
      write(*,*) "SINGLE PRECISION"
    endif
  endif

  ! save  in physical space
  call save_field_hdf5 ( time,fname_outx,vecfld(:,:,:,1) )
  call save_field_hdf5 ( time,fname_outy,vecfld(:,:,:,2) )
  call save_field_hdf5 ( time,fname_outz,vecfld(:,:,:,3) )

  if (mpirank==0) then
    write(*,*) "Done writing output!"
    write(*,'(80("-"))')
  endif

  deallocate (vecfld)
  call fft_free()
end subroutine post_grad
