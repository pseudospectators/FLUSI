!-------------------------------------------------------------------------------
! ./flusi --postprocess --spectrum ux_00000.h5 uy_00000.h5 uz_00000.h5 spectrum.dat
!-------------------------------------------------------------------------------
! NOTE: I actually did not figure out what happens if xl=yl=zl/=2*pi
! which is a rare case in all isotropic turbulence situtations, and neither
! of the corresponding routines have been tested for that case.
subroutine post_spectrum(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use helpers

  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, spectrum_file
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time, E
  real(kind=pr), dimension(:), allocatable :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec
  integer :: mpicode, k

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --spectrum ux_00000.h5 uy_00000.h5 uz_00000.h5 spectrum.dat"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,spectrum_file)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  if (mpirank==0) then
    write(*,'("Computing spectrum of ",A,1x,A,1x,A)') &
    trim(adjustl(fname_ux)), trim(adjustl(fname_uy)), trim(adjustl(fname_uz))
    write(*,'("Spectrum will be written to ",A)') trim(adjustl(spectrum_file))
    write(*,'("Resolution is ",i4,1x,i4,1x,i4)') nx, ny, nz
    write(*,'("Domain size is", es12.4,1x,es12.4,1x,es12.4)') xl, yl ,zl
  endif

  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("mpirank=",i5," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
  &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
  mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)


  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(S_Ekinx(0:nx-1),S_Ekiny(0:nx-1),S_Ekinz(0:nx-1),S_Ekin(0:nx-1),kvec(0:nx-1))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft3( inx=u, outk=uk)

  ! compute the actual spectrum
  call compute_spectrum( time,kvec,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

  ! compute energy in physical space (to check parsevals identity)
  u(:,:,:,1) = u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2
  E = mpisum( sum(u(:,:,:,1))/2.d0 )

  ! on root, write it to disk
  if (mpirank == 0) then
    open(10,file=spectrum_file,status='replace')
    write(10,'(5(A15,1x))') '%   K ','E_u(K)','E_ux(K)','E_uy(K)','E_uz(K)'
    do k=0,nx-1
      write(10,'(5(1x,es15.8))') kvec(k),S_Ekin(k),S_Ekinx(k),S_Ekiny(k),S_Ekinz(k)
      write(*,'(5(1x,es15.8))') kvec(k),S_Ekin(k),S_Ekinx(k),S_Ekiny(k),S_Ekinz(k)
    enddo

    write(10,*) '% parsevals identity:'
    write(10,*) '% Energy integrated from spectrum: ',sum(S_ekin)
    write(10,*) '% Energy in phys. space:', E/(dble(nx)*dble(ny)*dble(nz))
    write(10,*) '% sum = ', sum(S_Ekin), sum(S_Ekinx), sum(S_Ekiny), sum(S_Ekinz)
    write(10,*) '% sum*dk = integral = ', sum(S_Ekin)*(kvec(2)-kvec(1))
    write(10,*) '% time = ',time
    close(10)
  endif

  deallocate (u)
  deallocate (uk)
  deallocate(S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec)
  call fft_free()

end subroutine post_spectrum
