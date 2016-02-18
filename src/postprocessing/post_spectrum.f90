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
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, spectrum_file
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time, sum_u
  real(kind=pr), dimension(:), allocatable :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin
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

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
    write (*,*) "Error in arguments, files do not start with ux uy and uz"
    write (*,*) "note files have to be in the right order"
    call abort()
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  if (mpirank==0) then
    write(*,'("Computing spectrum of ",A,1x,A,1x,A)') &
    trim(adjustl(fname_ux)), trim(adjustl(fname_uy)), trim(adjustl(fname_uz))
    write(*,'("Spectrum will be written to ",A)') trim(adjustl(spectrum_file))
    write(*,'("Resolution is ",i4,1x,i4,1x,i4)') nx, ny, nz
    write(*,'("Domain size is", es12.4,1x,es12.4,1x,es12.4)') xl, yl ,zl
  endif

  call fft_initialize() ! also initializes the domain decomp


  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("mpirank=",i5," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
  &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
  mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)


  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(S_Ekinx(0:nx-1),S_Ekiny(0:nx-1),S_Ekinz(0:nx-1),S_Ekin(0:nx-1))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft (uk(:,:,:,1),u(:,:,:,1))
  call fft (uk(:,:,:,2),u(:,:,:,2))
  call fft (uk(:,:,:,3),u(:,:,:,3))

  ! compute the actual spectrum
  call compute_spectrum( time,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

  ! on root, write it to disk
  if (mpirank == 0) then
    open(10,file=spectrum_file,status='replace')
    write(10,'(5(A15,1x))') '%   K ','E_u(K)','E_ux(K)','E_uy(K)','E_uz(K)'
    do k=0,nx-1
      write(10,'(5(1x,es15.8))') dble(k),S_Ekin(k),S_Ekinx(k),S_Ekiny(k),S_Ekinz(k)
    enddo

    sum_u=0.0d0
    do k=1,nx-1
      sum_u=sum_u +S_Ekin(k)
    enddo
    write(10,*) '% Etot = ',sum_u
    write(10,*) '% time = ',time
    close(10)
  endif

  deallocate (u)
  deallocate (uk)
  deallocate(S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
  call fft_free()

end subroutine post_spectrum
