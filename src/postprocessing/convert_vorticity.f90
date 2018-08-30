subroutine convert_vorticity(help)
  use vars
  use p3dfft_wrapper
  use module_helpers
  use basic_operators

  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, order, mode
  character(len=strlen) :: prefix, fname_outx, fname_outy, fname_outz
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  complex(kind=pr),dimension(:,:,:),allocatable :: workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, workr, workr2
  real(kind=pr),dimension(:,:,:),allocatable :: w
  real(kind=pr) :: time, divu_max

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p [--vor-abs | --vor-abs-FD | --vorticity | --vorticity-FD | --Q | --Q-FD]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Read velocity components from file and compute their curl or Q-criterion."
    write(*,*) " Optionally second order or non-periodic using finite differences."
    write(*,*) " "
    write(*,*) " Compute vorticity vector (fourier)"
    write(*,*) " ./flusi -p --vorticity ux_00.h5 uy_00.h5 uz_00.h5 vorx_00.h5 vory_00.h5 vorz_00.h5 [--second-order]"
    write(*,*) " "
    write(*,*) " Directly compute vorticity magnitude: (fourier)"
    write(*,*) " ./flusi -p --vor-abs ux_00.h5 uy_00.h5 uz_00.h5 vorabs_00.h5 [--second-order]"
    write(*,*) " "
    write(*,*) " Finite difference versions: (periodic, 2nd and 4th order or non-periodic 2nd order)"
    write(*,*) " ./flusi -p --vor-abs-FD ux_00.h5 uy_00.h5 uz_00.h5 vorabs_00.h5 [ORDER]"
    write(*,*) " ./flusi -p --vorticity-FD ux_00.h5 uy_00.h5 uz_00.h5 vorx_00.h5 vory_00.h5 vorz_00.h5 [ORDER]"
    write(*,*) " where ORDER is one of: "
    write(*,*) " --second-order"
    write(*,*) " --fourth-order"
    write(*,*) " --second-order-nonper"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "  Q-criterion"
    write(*,*) " ./flusi -p --Q-FD ux_00.h5 uy_00.h5 uz_00.h5 Q_00.h5 [ORDER]"
    write(*,*) ""
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes, fully"
    return
  endif

  call get_command_argument(2,mode)
  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)

  select case (mode)
  case ("--vorticity","--vorticity-FD")
    call get_command_argument(6,fname_outx)
    call get_command_argument(7,fname_outy)
    call get_command_argument(8,fname_outz)
    call get_command_argument(9,order)

  case ("--vor-abs", "--vor-abs-FD", "--Q", "--Q-FD")
    call get_command_argument(6,fname_outx)
    call get_command_argument(7,order)
    fname_outy = ""
    fname_outz = ""

  case default
    call abort(123321,"unkown in convert_vorticity.f90")
  end select

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "Computing vorticity from velocity given in these files: "
    write(*,'(80("-"))')
    write(*,*) "mode is: "//trim(adjustl(mode))
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(fname_outx))
    write(*,*) trim(adjustl(fname_outy))
    write(*,*) trim(adjustl(fname_outz))
    write(*,*) "Using order flag: "//trim(adjustl(order))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )

  !-----------------------------------------------------------------------------
  ! initialize code and scaling factors for derivatives, also domain decomposition
  !-----------------------------------------------------------------------------
  select case (mode)
    !*** Fourier
  case ("--vorticity", "--vor-abs", "--Q")
    ! spectral version
    if (root) write(*,*) "using spectral code (FOURIER)"

    ! initialize FFTS
    call fft_initialize()

    ! alloc memory for fourier transform
    allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

    !*** Finite-differences
  case ("--vorticity-FD", "--vor-abs-FD", "--Q-FD")
    ng = 2
    if (order=="--fourth-order") order="centered_4th"
    if (order=="--second-order") order="centered_2nd"
    ! in the non-periodic case, do not use ghost nodes at all (we use MPI_TRANSPOSES)
    if (order=="--second-order-nonper") ng = 0
    ! finite differences
    if (root) write(*,*) "using finite differences code"
    ! set number of ghost points
    if (root) write(*,'("Set up ng=",i1," ghost points")') ng

    ! initialize domain decomposition, but do not use FFTs
    if (ng == 0) then
      ! if non-periodic, 1d decomposition, if possible, is more efficient
      call decomposition_initialize()
    else
      ! force code to use 2d decomposition if using ghost nodes -> higher
      ! efficiency
      call decomposition_initialize( .true. )
    endif

    call setup_cart_groups()

  end select

  if (root) write(*,*) "order flag is: "//trim(adjustl(order))

  !-----------------------------------------------------------------------------
  ! read input data
  !-----------------------------------------------------------------------------
  ! allocate memory for input data
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  call read_single_file( fname_ux, u(:,:,:,1) )
  call read_single_file( fname_uy, u(:,:,:,2) )
  call read_single_file( fname_uz, u(:,:,:,3) )


  !-----------------------------------------------------------------------------
  ! compute curl, possibly with second order filter
  !-----------------------------------------------------------------------------
  select case (mode)
  case ("--vorticity","--vor-abs")
    !*****************************************************************************
    !*** Fourier
    !*****************************************************************************
    ! to Fourier space
    call fft3 (inx=u, outk=uk)

    ! actually compute curl, possibly with 2nd order filter (PERIODIC)
    if (order=="--second-order") then
      if (root) write(*,*) "using SECOND ORDER accuracy"
      call curl_2nd(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
    else
      if (root) write(*,*) "using SPECTRAL accuracy"
      call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
    endif

    ! back to x-space
    call ifft3(ink=uk, outx=u)

    ! free memory
    deallocate(uk)

    if (mode=="--vor-abs") then
      ! vorticity magnitude
      u(:,:,:,1) = sqrt( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)
      call save_field_hdf5( time, fname_outx, u(:,:,:,1) )
    else
      ! now u contains the vorticity in physical space
      call save_field_hdf5( time, fname_outx, u(:,:,:,1) )
      call save_field_hdf5( time, fname_outy, u(:,:,:,2) )
      call save_field_hdf5( time, fname_outz, u(:,:,:,3) )
    endif

    ! free memory
    deallocate (u)

    ! bye bye
    call fft_free()

  case ("--vorticity-FD", "--vor-abs-FD" )
    !*****************************************************************************
    !*** Finite-differences
    !*****************************************************************************
    ! finite differences
    if (root) write(*,*) "using finite differences code"

    ! allocate work array WITH ghost nodes. Note in the case "order==--second-order-nonper"
    ! ng=0 and thus ra=ga and rb=gb (no ghost nodes used, we use MPI_TRANSPOSES)
    allocate(workr(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3))
    allocate(workr2(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3))

    ! fill interior of work array wth input data. ghost node synching is done,
    ! if the periodic case is used, in the curl-subroutine.
    workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),:) = u
    deallocate(u)

    if (order=="--second-order-nonper") then
      if (root) write(*,*) "NON-PERIODIC DISCRETIZATION USED"
      call curl_FD_nonper( workr, workr2, order )
    else
      if (root) write(*,*) "PERIODIC DISCRETIZATION USED"
      call curl_FD( workr, workr2, order )
    endif

    ! as we have the curl, we can free the data array
    deallocate(workr)

    ! save either vorticity vector in 3 files or just the magnitude in one file
    if (mode=="--vor-abs-FD") then
      ! vorticity magnitude
      workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) = sqrt( workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1)**2 &
      + workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2)**2 &
      + workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3)**2)
      call save_field_hdf5( time, fname_outx, workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
    else
      ! now work2 contains the vorticity in physical space
      call save_field_hdf5( time, fname_outx, workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
      call save_field_hdf5( time, fname_outy, workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2) )
      call save_field_hdf5( time, fname_outz, workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3) )
    endif

    ! free memory
    deallocate( workr2 )

  case ("--Q-FD" )
    !*****************************************************************************
    !*** Finite-differences, Q-criterion
    !*****************************************************************************
    ! finite differences
    if (root) write(*,*) "using finite differences code"

    ! allocate work array WITH ghost nodes. Note in the case "order==--second-order-nonper"
    ! ng=0 and thus ra=ga and rb=gb (no ghost nodes used, we use MPI_TRANSPOSES)
    allocate(workr(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3))
    allocate(workr2(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1))

    ! fill interior of work array wth input data. ghost node synching is done,
    ! if the periodic case is used, in the curl-subroutine.
    workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),:) = u
    deallocate(u)

    if (order=="--second-order-nonper") then
      if (root) write(*,*) "NON-PERIODIC DISCRETIZATION USED"
      call Q_FD_nonper( workr, workr2(:,:,:,1), order )
    else
      if (root) write(*,*) "PERIODIC DISCRETIZATION USED"
      call Q_FD( workr, workr2(:,:,:,1), order )
    endif

    ! as we have the curl, we can free the data array
    deallocate(workr)

    ! save data
    call save_field_hdf5( time, fname_outx, workr2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )

    ! free memory
    deallocate( workr2 )

  end select
end subroutine convert_vorticity
