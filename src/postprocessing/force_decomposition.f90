!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine force_decomposition(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfile
  character(len=strlen) :: fname_potx, fname_poty, fname_potz, order
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, upot, vor
  real(kind=pr) :: time, ux,uy,uz,upotx,upoty,upotz,vorx,vory,vorz
  integer :: ix,iy,iz, mpicode

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --force-decomp ux_00.h5 uy_00.h5 uz_00.h5 upotx_00.h5 upoty_00.h5 upotz_00.h5"
    write(*,*) " output.h5 --second-order"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "!"
    write(*,*) "!"
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
  call get_command_argument(6,fname_potx)
  call get_command_argument(7,fname_poty)
  call get_command_argument(8,fname_potz)
  call get_command_argument(9,outfile)
  call get_command_argument(10,order)


  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "Force decomposition [1]"
    write(*,*) "[1] C.C. Chang: Potential flow and forces for incompressible viscous &
               &flow. Proc. R. Soc. Lond. A (1992), 437, 517--525"
    write(*,'(80("-"))')
    write(*,*) "DNS flow field data read from:"
    write(*,'("udnsx=",A)') trim(adjustl(fname_ux))
    write(*,'("udnsy=",A)') trim(adjustl(fname_uy))
    write(*,'("udnsz=",A)') trim(adjustl(fname_uz))
    write(*,*) "Potential flow read from:"
    write(*,'("upotx=",A)') trim(adjustl(fname_potx))
    write(*,'("upoty=",A)') trim(adjustl(fname_poty))
    write(*,'("upotz=",A)') trim(adjustl(fname_potz))
    write(*,'("Order flag:",A)') trim(adjustl(order))
    write(*,'("Output=",A)') trim(adjustl(outfile))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )
  call check_file_exists( fname_potx )
  call check_file_exists( fname_poty )
  call check_file_exists( fname_potz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )

  scalex = 2.d0*pi/xl
  scaley = 2.d0*pi/yl
  scalez = 2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  neq=3
  nd=3
  ncw=3
  nrw=3

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(upot(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  ! read vorticity from files to u
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call read_single_file ( fname_potx, upot(:,:,:,1) )
  call read_single_file ( fname_poty, upot(:,:,:,2) )
  call read_single_file ( fname_potz, upot(:,:,:,3) )

  if (root) write(*,*) "done reading input files."

  ! we need to compute the curl of the DNS data. this assumes by the way that the
  ! data is periodic (so one can run into trouble when using subsets of a field)
  call fft3( inx=u,outk=uk ) ! uk is now vork
  if (order=="--second-order") then
    if (root) write(*,*) "using second order accuracy!"
    call curl_2nd( uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3) )
  else
    if (root) write(*,*) "using spectral accuracy!"
    call curl( uk )
  endif
  call ifft3( ink=uk, outx=vor )
  deallocate (uk)

  if (root) write(*,*) "done computing curl(u_dns)."

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ! local loop variables
        ux   = u(ix,iy,iz,1)
        uy   = u(ix,iy,iz,2)
        uz   = u(ix,iy,iz,3)

        vorx = vor(ix,iy,iz,1)
        vory = vor(ix,iy,iz,2)
        vorz = vor(ix,iy,iz,3)

        upotx = upot(ix,iy,iz,1)
        upoty = upot(ix,iy,iz,2)
        upotz = upot(ix,iy,iz,3)

        upot(ix,iy,iz,1) = (vorz*uy-vory*uz)*upotx &
                         + (vorx*uz-vorz*ux)*upoty &
                         + (vory*ux-vorx*uy)*upotz
        ! invert sign, according to Chang's work
        upot(ix,iy,iz,1) = - upot(ix,iy,iz,1)
      enddo
    enddo
  enddo

  if (root) write(*,*) "done computing force decomposition."
  if (root) write(*,*) "writing output"
  call save_field_hdf5 ( time, outfile, upot(:,:,:,1) )

  deallocate (u)
  deallocate (upot, vor)
  call fft_free()

end subroutine force_decomposition
