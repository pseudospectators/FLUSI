!-------------------------------------------------------------------------------
! ./flusi --postprocess --turbulence-analysis ux_00000.h5 uy_00000.h5 uz_00000.h5 nu outfile.dat
! NOTE: I actually did not figure out what happens if xl=yl=zl/=2*pi
! which is a rare case in all isotropic turbulence situtations, and neither
! of the corresponding routines have been tested for that case.
!-------------------------------------------------------------------------------
subroutine turbulence_analysis(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers
  
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, viscosity, outfile
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk,vork
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, vor
  real(kind=pr) :: time, epsilon_loc, epsilon, fact, E, u_rms,lambda_macro,lambda_micro
  real(kind=pr), dimension(:), allocatable :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec
  integer :: ix,iy,iz, mpicode
  real(kind=pr)::kx,ky,kz,kreal

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --turbulence-analysis ux_00000.h5 uy_00000.h5 uz_00000.h5 nu outfile.dat"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "compute a bunch of values relevant to Homogeneous isotropic turbulence and write them to outfile"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  neq = 3
  nd = 3

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,viscosity)
  call get_command_argument(7,outfile)
  read(viscosity,*) nu

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    call postprocessing_ascii_header(17)
    write(17,'(A)') "-----------------------------------"
    write(17,'(A)') "FLUSI turbulence analysis"
    write(17,'("call: ./flusi -p --turbulence-analysis ",5(A,1x))') trim(adjustl(fname_ux)),&
    trim(adjustl(fname_uy)),trim(adjustl(fname_uz)),trim(adjustl(viscosity)),trim(adjustl(outfile))
    write(17,'(A)') "-----------------------------------"
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(vork(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(S_Ekinx(0:nx-1),S_Ekiny(0:nx-1),S_Ekinz(0:nx-1),S_Ekin(0:nx-1),kvec(0:nx-1))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft3 (inx=u,outk=uk)
  call curl (uk,vork)
  call ifft3 (ink=vork,outx=vor)

  ! compute spectrum
  call compute_spectrum( time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

  !-----------------------------------------------------------------------------
  ! dissipation rate from velocity in Fourier space
  !-----------------------------------------------------------------------------
  call dissipation_rate(uk, epsilon)
  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') epsilon, "Dissipation rate from velocity in Fourier space"
  endif

  !-----------------------------------------------------------------------------
  ! dissipation rate from vorticty
  !-----------------------------------------------------------------------------
  epsilon_loc = nu * sum(vor(:,:,:,1)**2+vor(:,:,:,2)**2+vor(:,:,:,3)**2)
  call MPI_REDUCE(epsilon_loc,epsilon,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
  MPI_COMM_WORLD,mpicode)

  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') epsilon/(dble(nx)*dble(ny)*dble(nz)), "Dissipation rate from vorticity"
  endif

  !-----------------------------------------------------------------------------
  ! dissipation rate from spectrum; see Ishihara, Kaneda "High
  ! resolution DNS of incompressible Homogeneous forced turbulence -time dependence
  ! of the statistics" or my thesis
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    epsilon=0.0
    do ix = 0,nx-1
      epsilon = epsilon + 2.d0 * nu * dble(ix**2) * S_Ekin(ix)
    enddo
    write(17,'(g15.8,5x,A)') epsilon, "Dissipation rate from spectrum"
  endif

  !-----------------------------------------------------------------------------
  ! energy from velocity
  !-----------------------------------------------------------------------------
  epsilon_loc = 0.5d0*sum(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2)
  call MPI_REDUCE(epsilon_loc,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
  MPI_COMM_WORLD,mpicode)

  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') E/(dble(nx)*dble(ny)*dble(nz)), "energy from velocity"
  endif

  !-----------------------------------------------------------------------------
  ! energy from spectrum
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    E=0.0
    do ix = 0,nx-1
      E = E + S_Ekin(ix)
    enddo
    write(17,'(g15.8,5x,A)') E, "energy from spectrum"
  endif

  u_rms=dsqrt(2.d0*E/3.d0)
  !-----------------------------------------------------------------------------
  ! kolmogrov scales
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') (nu**3 / epsilon)**(0.25d0), "kolmogorov length scale"
    write(17,'(g15.8,5x,A)') (nu / epsilon)**(0.5d0), "kolmogorov time scale"
    write(17,'(g15.8,5x,A)') (nu*epsilon)**(0.25d0), "kolmogorov velocity scale"
    write(17,'(g15.8,5x,A)') u_rms, "RMS velocity"
  endif

  !-----------------------------------------------------------------------------
  ! taylor scales
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    lambda_micro = (15.d0*nu*u_rms**2 / epsilon)**(0.5d0)
    lambda_macro=0.0
    do ix = 1,nx-1
      lambda_macro = lambda_macro + pi/(2.d0*u_rms**2) * S_Ekin(ix) / dble(ix)
    enddo
    write(17,'(g15.8,5x,A)') lambda_micro, "Taylor micro scale"
    write(17,'(g15.8,5x,A)') lambda_macro, "Taylor macro scale"
    write(17,'(g15.8,5x,A)') u_rms*lambda_macro/nu, "Reynolds Taylor macro scale"
    write(17,'(g15.8,5x,A)') u_rms*lambda_micro/nu, "Reynolds Taylor micro scale"
    write(17,'(g15.8,5x,A)') lambda_macro/u_rms, "eddy turnover time"
    write(17,'(g15.8,5x,A)') (2./3.)*(dble(nx/2-1)), "kmax"
    write(17,'(g15.8,5x,A)') (2./3.)*(dble(nx/2-1))*(nu**3 / epsilon)**(0.25d0), "kmax*eta"
  endif

  if(mpirank==0) close(17)

  deallocate (u,vor,vork)
  deallocate (uk)
  deallocate (S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec)
  call fft_free()

end subroutine turbulence_analysis
