subroutine post_CVE(help)
  use vars
  use p3dfft_wrapper
  use flusi_wavelet_lib
  use module_helpers

  implicit none

  logical, intent(in) :: help
  integer :: save_coh, save_inc, iterations
  real(kind=pr) :: time
  real(kind=pr), dimension(:,:,:,:), allocatable :: work1, work2, work3
  type(orth_wavelet) :: wavelet
  character(len=strlen) :: mode, fname_ux, fname_uy, fname_uz


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --coherent-vortex-extraction vorx_00.h5 vory_00.h5 vorz_00.h5 save_coh save_inc iterations"
    write(*,*) "./flusi -p --CVE vorx_00.h5 vory_00.h5 vorz_00.h5 save_coh save_inc iterations"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --coherent-scalar-extraction phi_00.h5 save_coh save_inc iterations"
    write(*,*) "./flusi -p --CSE phi_00.h5 save_coh save_inc iterations"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! Coherent vortex extraction, split given vorticity field into coherent and incoherent "
    write(*,*) "! part. "
    write(*,*) "! "
    write(*,*) "! Saveflag: [COH,INC]"
    write(*,*) "! Iterations: maximum number of iterattions performed"
    write(*,*) "! "
    write(*,*) "! OUTPUT: coh*_*.h5 and inc*_*.h5 if saveflag is 11"
    write(*,*) "! "
    write(*,*) "! Coiflet12 is used (orthogonal, periodized)"
    write(*,*) "! -------------------------------------------------------------------------"
    write(*,*) "! Algorithm:"
    write(*,*) "!   In general, the coherent part is what remains after denoising. In order"
    write(*,*) "!   to denoise, one has to know the variance of the noise (=incoherent"
    write(*,*) "!   part), but that is a priori unknown. Therefore, an iterative method is"
    write(*,*) "!   employed. We start assuming the entire signal is just noise, denoise,"
    write(*,*) "!   and take the variance of what remains as better estimate."
    write(*,*) "! -------------------------------------------------------------------------"
    write(*,*) "! Original version for 1D signals by Antonio Merulla, Michele Caldoro"
    write(*,*) "! Modified for 3D CVE, 11/2015, T.Engels"
    write(*,*) "! Ported to flusi (from matlab) 07/2017, T. Engels"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  ! we can apply this method to vector or scalar problems.
  call get_command_argument(2,mode)
  if (mode=="--coherent-scalar-extraction" .or. mode == "--CSE") then
    ! scalar mode
    nd = 1
    call get_command_argument(3,fname_ux)
    call check_file_exists( fname_ux )

    call get_command_argument(4,mode)
    read(mode,*) save_coh
    call get_command_argument(5,mode)
    read(mode,*) save_inc
    call get_command_argument(6,mode)
    read(mode,*) iterations
  else
    ! vector mode
    nd = 3
    call get_command_argument(3,fname_ux)
    call get_command_argument(4,fname_uy)
    call get_command_argument(5,fname_uz)

    call check_file_exists( fname_ux )
    call check_file_exists( fname_uy )
    call check_file_exists( fname_uz )

    call get_command_argument(6,mode)
    read(mode,*) save_coh
    call get_command_argument(7,mode)
    read(mode,*) save_inc
    call get_command_argument(8,mode)
    read(mode,*) iterations
  endif


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()
  ! call fft_initialize()
  call setup_cart_groups()
  ! initialize wavelet transform
  call setup_coiflet_coefs( 1, wavelet )

  if (root) then
    write(*,*) "Coherent vortex/scalar extraction"
    write(*,'("We will perform i=",i2," iterations at most")') iterations
    write(*,'("save coherent=",i1," save incoherent=",i1)') save_coh, save_inc
  endif

  allocate( work1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) )
  allocate( work2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) )
  allocate( work3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) )

  if (nd == 3) then
    !-----------------------------------------------------------------------------
    ! vector mode
    !-----------------------------------------------------------------------------
    ! read input files
    call read_single_file( fname_ux, work1(:,:,:,1))
    call read_single_file( fname_uy, work1(:,:,:,2))
    call read_single_file( fname_uz, work1(:,:,:,3))
    ! perfom actual decomposition
    call coherent_vortex_extraction( wavelet, iterations, work1, work2, work3 )
    ! save coherent field, if desired
    if (save_coh == 1) then
      call save_field_hdf5( time, "coh"//fname_ux, work2(:,:,:,1))
      call save_field_hdf5( time, "coh"//fname_uy, work2(:,:,:,2))
      call save_field_hdf5( time, "coh"//fname_uz, work2(:,:,:,3))
    endif
    ! save incoherent field, if desired
    if (save_inc == 1) then
      call save_field_hdf5( time, "inc"//fname_ux, work3(:,:,:,1))
      call save_field_hdf5( time, "inc"//fname_uy, work3(:,:,:,2))
      call save_field_hdf5( time, "inc"//fname_uz, work3(:,:,:,3))
    endif

  else
    !-----------------------------------------------------------------------------
    ! scalar mode
    !-----------------------------------------------------------------------------
    ! read only one input file
    call read_single_file( fname_ux, work1(:,:,:,1))
    ! perfom actual decomposition
    call coherent_scalar_extraction( wavelet, iterations, work1, work2, work3 )
    ! save coherent field
    if (save_coh == 1) then
      call save_field_hdf5( time, "coh"//fname_ux, work2)
    endif
    ! save incoherent field
    if (save_inc == 1) then
      call save_field_hdf5( time, "inc"//fname_ux, work3)
    endif
  endif

  deallocate( work1, work2, work3 )
end subroutine post_CVE
