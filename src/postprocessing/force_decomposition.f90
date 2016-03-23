!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine force_decomposition(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfiles_given
  character(len=strlen) :: fname_outx, fname_outy, fname_outz, prefix
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, upot, vor
  real(kind=pr) :: time, ux,uy,uz,upotx,upoty,upotz,vorx,vory,vorz
  integer :: ix,iy,iz, mpicode

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --force-decomp ux_00.h5 uy_00.h5 uz_00.h5 upotx_00.h5 upoty_00.h5 upotz_00.h5"
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

  call get_command_argument(6,fname_outx)
  call get_command_argument(7,fname_outy)
  call get_command_argument(8,fname_outz)


  ! ! header and information
  ! if (mpirank==0) then
  !   write(*,'(80("-"))')
  !   write(*,*) "vor2u (Biot-Savart)"
  !   write(*,*) "Computing velocity from vorticity given in these files: "
  !   write(*,'(80("-"))')
  !   write(*,*) trim(adjustl(fname_ux))
  !   write(*,*) trim(adjustl(fname_uy))
  !   write(*,*) trim(adjustl(fname_uz))
  !   write(*,*) "Writing to:"
  !   write(*,*) trim(adjustl(fname_outx))
  !   write(*,*) trim(adjustl(fname_outy))
  !   write(*,*) trim(adjustl(fname_outz))
  !   write(*,'(80("-"))')
  !   if ((fname_ux(1:4).ne."vorx").or.(fname_uy(1:4).ne."vory").or.(fname_uz(1:4).ne."vorz")) then
  !     write (*,*) "WARNING in arguments, files do not start with vorx vory and vorz"
  !     write (*,*) "note files have to be in the right order"
  !   endif
  ! endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )
  call check_file_exists( fname_outx )
  call check_file_exists( fname_outy )
  call check_file_exists( fname_outz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu )

  pi = 4.d0 *datan(1.d0)
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

  call read_single_file ( fname_outx, upot(:,:,:,1) )
  call read_single_file ( fname_outy, upot(:,:,:,2) )
  call read_single_file ( fname_outz, upot(:,:,:,3) )

  call fft3( inx=u,outk=uk ) ! uk is now vork
  call curl( uk )
  call ifft3( ink=uk, outx=vor )
  deallocate (uk)

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
      enddo
    enddo
  enddo

  call save_field_hdf5 ( time, 'result_00.h5', upot(:,:,:,1) )

  deallocate (u)
  deallocate (upot, vor)
  call fft_free()

end subroutine force_decomposition
