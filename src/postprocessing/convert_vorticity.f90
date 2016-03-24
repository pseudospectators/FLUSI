!-------------------------------------------------------------------------------
! Read velocity components from file and compute their curl. Optionally second order
!
! You can write to standard output: (vorx_0000.h5 in the example:)
! ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 [--second-order]
!
! Or specify a prefix for the output files: (writes to curl_0000.h5 in example:)
! ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --outputprefix curl [--second-order]
!
! Or specifiy the ouput files directly:
! ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --outfiles vx_00.h5 vy_00.h5 vz_00.h5 [--second-order]
!
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! can be done in parallel. the flag --second order can be used for filtering
subroutine convert_vorticity(help)
  use vars
  use p3dfft_wrapper
  use helpers
  use basic_operators
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, order
  character(len=strlen) :: prefix, fname_outx, fname_outy, fname_outz
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  complex(kind=pr),dimension(:,:,:),allocatable :: workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr),dimension(:,:,:),allocatable :: workr
  real(kind=pr) :: time, divu_max

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Read velocity components from file and compute their curl. Optionally second order"
    write(*,*) " "
    write(*,*) " You can write to standard output: (vorx_0000.h5 in the example:)"
    write(*,*) " ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 [--second-order]"
    write(*,*) " "
    write(*,*) " Or specify a prefix for the output files: (writes to curl_0000.h5 in example:)"
    write(*,*) " ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --outputprefix curl [--second-order]"
    write(*,*) " "
    write(*,*) " Or specifiy the ouput files directly:"
    write(*,*) " ./flusi -p --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --outfiles vx_00.h5 vy_00.h5 vz_00.h5 [--second-order]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,order)

  if (order == "--outputprefix") then
    call get_command_argument(7,prefix)
    call get_command_argument(8,order)
    fname_outx=trim(adjustl(prefix))//"x"//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
    fname_outy=trim(adjustl(prefix))//"y"//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
    fname_outz=trim(adjustl(prefix))//"z"//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)

  elseif (order == "--outfiles") then
    call get_command_argument(7,fname_outx)
    call get_command_argument(8,fname_outy)
    call get_command_argument(9,fname_outz)
    call get_command_argument(10,order)

  else
    fname_outx='vorx'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
    fname_outy='vory'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
    fname_outz='vorz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)
  endif

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "vor2u (Biot-Savart)"
    write(*,*) "Computing vorticity from velocity given in these files: "
    write(*,'(80("-"))')
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

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu )

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  ! to Fourier space
  call fft3 (inx=u, outk=uk)

  !-----------------------------------------------------------------------------
  ! compute divergence of input fields and show the maximum value
  !-----------------------------------------------------------------------------
  allocate( workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)) )
  call divergence( uk, workc)
  call ifft( ink=workc, outx=workr )
  divu_max = fieldmax(workr)
  if(mpirank==0) write(*,'("maximum divergence in input field=",es12.4)') divu_max
  deallocate(workr,workc)

  !-----------------------------------------------------------------------------
  ! compute curl, possibly with second order filter
  !-----------------------------------------------------------------------------
  if (order=="--second-order") then
    if (mpirank==0) write(*,*) "using SECOND ORDER accuracy"
    call curl_2nd(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  else
    if (mpirank==0) write(*,*) "using SPECTRAL accuracy"
    call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  endif

  call ifft3 (ink=uk, outx=u)

  ! now u contains the vorticity in physical space
  call save_field_hdf5 ( time,fname_outx,u(:,:,:,1) )
  call save_field_hdf5 ( time,fname_outy,u(:,:,:,2) )
  call save_field_hdf5 ( time,fname_outz,u(:,:,:,3) )

  deallocate (u)
  deallocate (uk)
  call fft_free()

end subroutine convert_vorticity
