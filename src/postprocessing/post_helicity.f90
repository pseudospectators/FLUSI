!-------------------------------------------------------------------------------
! compute helicity from velocity field with spectral accuracy
!
!-------------------------------------------------------------------------------
subroutine post_helicity(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfile, normalized
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk, vork
  complex(kind=pr),dimension(:,:,:),allocatable :: helk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, vor
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, divu_max, errx, erry, errz
  integer :: ix,iy,iz, mpicode, k, nk
  real(kind=pr) :: kx, ky, kz, kreal, kmax, dk
  real(kind=pr), dimension(:),allocatable :: SpecHel, kvec
  real(kind=pr), dimension(:),allocatable  :: SpecHel_loc

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --helicity ux_000.h5 uy_000.h5 uz_000.h5 helicity_000.h5 [--normalized]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! compute (normalized) helicity from given velocity field"
    write(*,*) "! employs spectral precision"
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
  call get_command_argument(6,outfile)
  call get_command_argument(7,normalized)


  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "helicity computation"
    write(*,*) "Computing helicit from velocity given in these files: "
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(outfile))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  neq=3
  nd=3
  ncw=3
  nrw=3

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(vork(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(helk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  ! read vorticity from files to u
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )


  call fft3( inx=u, outk=uk )
  call curl( ink=uk, outk=vork )

  !--------------SpecHel
  helk = uk(:,:,:,1)*conjg(vork(:,:,:,1)) &
       + uk(:,:,:,2)*conjg(vork(:,:,:,2)) &
       + uk(:,:,:,3)*conjg(vork(:,:,:,3))

  allocate( SpecHel(0:nx-1),SpecHel_loc(0:nx-1),kvec(0:nx-1) )
  SpecHel_loc=0.d0

  kmax = minval( (/scalex*dble(nx/2),scaley*dble(ny/2),scalez*dble(nz/2)/) )
  dk = minval( (/scalex,scaley,scalez/) )
  nk = nint(kmax / dk)

  do iz=ca(1),cb(1)
    kz = wave_z(iz)
    do iy=ca(2),cb(2)
      ky = wave_y(iy)
      do ix=ca(3),cb(3)
        kx = wave_x(ix)
        ! compute 2-norm of wavenumber, scaled to our actual domain size
        kreal = dsqrt( (kx*kx)+(ky*ky)+(kz*kz))
        ! note spectrum omits parts of the data (circle in square problem)
        if (kreal <= kmax) then
          k = nint( kreal/dk )
          SpecHel_loc(k)  = SpecHel_loc(k)  + dble(real(helk(iz,iy,ix)))
        endif
      enddo
    enddo
  enddo

  ! the returned wavenumber is the middle of the bin: (and not kreal!)
  do k = 0, nk
    kvec(k) = dble(k)*dk
  enddo

  call MPI_ALLREDUCE(SpecHel_loc,SpecHel,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

  if(root) then
    call init_empty_file(trim(adjustl(outfile))//'heli_spectrum.txt')
    open(14,file=trim(adjustl(outfile))//'heli_spectrum.txt',status='unknown',position='append')
    do ix = 0, nx-1
      write (14,'(2(es15.8,1x))') kvec(ix), SpecHel(ix)
    enddo
    close(14)
  endif
  deallocate(SpecHel, SpecHel_loc)
!--------------SpecHel

  call ifft3( ink=vork, outx=vor )

  if (normalized == "--normalized" ) then
    if (mpirank==0) write(*,*) "computing normalized helicity"
    call helicity_norm( u, vor, work)
  else
    if (mpirank==0) write(*,*) "computing absolute helicity"
    call helicity( u, vor, work)
  endif

  ! now u contains the velocity in physical space
  call save_field_hdf5 ( time, outfile, work )

  deallocate (u,uk,vor,work,helk,vork)
  call fft_free()
end subroutine post_helicity
