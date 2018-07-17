!-------------------------------------------------------------------------------
! compute helicity from velocity field with spectral accuracy
!
!-------------------------------------------------------------------------------
subroutine post_helicity(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfile, normalized
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk, vork
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, vor
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, divu_max, errx, erry, errz
  integer :: ix,iy,iz, mpicode, k, nk
  real(kind=pr) :: kx, ky, kz, kreal, kmax, dk
  real(kind=pr), dimension(:),allocatable :: SpecHel, kvec
  real(kind=pr), dimension(:),allocatable  :: SpecHel_loc
  real(kind=pr) :: mean_helicity

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
    write(*,*) "Computing helicity from velocity given in these files: "
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

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
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

  ! read vorticity from files to u
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )


  call fft3( inx=u, outk=uk )
  call curl( ink=uk, outk=vork )
  call ifft3( ink=vork, outx=vor )

  !--------------SpecHel
  allocate( SpecHel(0:nx-1),SpecHel_loc(0:nx-1),kvec(0:nx-1) )
  SpecHel_loc=0.d0

  ! for details, see comments on energy spectrum
  ! kmax = minval( (/scalex*dble(nx/2),scaley*dble(ny/2),scalez*dble(nz/2)/) )
  kmax = norm2( (/scalex*dble(nx/2),scaley*dble(ny/2),scalez*dble(nz/2)/) )
  dk = minval( (/scalex,scaley,scalez/) )
  nk = nint(kmax / dk)

  if (nx <= nk) then
    call abort(1212,"for some reason, we have too many wavenumbers in spectrum..")
  endif

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
          if ( ix==0 .or. ix==nx/2 ) then
            SpecHel_loc(k)  = SpecHel_loc(k)  +(dble(real( uk(iz,iy,ix,1)*conjg(vork(iz,iy,ix,1)) )) &
                                              + dble(real( uk(iz,iy,ix,2)*conjg(vork(iz,iy,ix,2)) )) &
                                              + dble(real( uk(iz,iy,ix,3)*conjg(vork(iz,iy,ix,3)) )) ) /2.d0
          else
            SpecHel_loc(k)  = SpecHel_loc(k) + dble(real( uk(iz,iy,ix,1)*conjg(vork(iz,iy,ix,1)) )) &
                                             + dble(real( uk(iz,iy,ix,2)*conjg(vork(iz,iy,ix,2)) )) &
                                             + dble(real( uk(iz,iy,ix,3)*conjg(vork(iz,iy,ix,3)) ))
          endif
         endif
      enddo
    enddo
  enddo

  ! NOTE
  ! The integral helicity, as represented by the sum of specHel, is a factor of 2
  ! too small to give you the enire mean helicity as it is integrated in physical space
  ! Therefore, to compare plancherell, you have to multiply the helicity spectrum by 2.
  ! (we simply skip the negative wavenumber part)
  ! On the other hand, if you check abs(H(k)) / k*E(k), you have to omit the factor 2
  ! since you have to sum from -intfy to +infty

  ! the returned wavenumber is the middle of the bin: (and not kreal!)
  kvec=0.d0
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


!--------------SpecHel

  if (normalized == "--normalized" ) then
    if (mpirank==0) write(*,*) "computing normalized helicity"
    call helicity_norm( u, vor, work)
  else
    if (mpirank==0) write(*,*) "computing absolute helicity"
    call helicity( u, vor, work)
  endif

  mean_helicity = mpisum(sum(work)) / dble(nx*ny*nz)

  if (root) then
    write(*,*) "spectrum sum: ", sum(SpecHel)
    write(*,*) "mean helicity: ", mean_helicity
    write(*,*) "cactoe: ", mean_helicity / sum(SpecHel)
  endif


  ! now u contains the velocity in physical space
  call save_field_hdf5 ( time, outfile, work )

  deallocate(SpecHel, SpecHel_loc)
  deallocate (u,uk,vor,work,vork)
  call fft_free()
end subroutine post_helicity
