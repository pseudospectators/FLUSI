!-------------------------------------------------------------------------------
! compute divergence of a scalar field
!
!-------------------------------------------------------------------------------
subroutine post_div(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy ,fname_uz, pr_out_flag, fname_out

  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  complex(kind=pr),dimension(:,:,:),allocatable :: divuk
  real(kind=pr),dimension(:,:,:),allocatable :: divu

  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --divergence ux_00.h5 uy_00.h5 uz_00.h5 div_00.h5 [--double-precision]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! compute divergence of a scalar field"
    write(*,*) "! with spectral accuracy"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,fname_out)
  call get_command_argument(7,pr_out_flag)

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "divergence computation"
    write(*,*) "Computing divergence of a vector field given in this file: "
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Writing divergence to:"
    write(*,*) trim(adjustl(fname_out))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  ! allocate memory
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))

  ! read vector field from file
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  ! Fourier transform
  call fft3( inx=u, outk=uk )
  deallocate (u)

  ! divergence in Fourier space
  allocate(divuk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  call divergence(uk, divuk)
  deallocate(uk)

  ! inverse Fourier transform
  allocate(divu(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call ifft( ink=divuk, outx=divu )
  deallocate(divuk)

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
  call save_field_hdf5 ( time, fname_out, divu )

  if (mpirank==0) then
    write(*,*) "Done writing output!"
    write(*,'(80("-"))')
  endif

  deallocate (divu)
  call fft_free()
end subroutine post_div
