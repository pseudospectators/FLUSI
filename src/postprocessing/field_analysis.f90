
!-------------------------------------------------------------------------------
! ./flusi -p --field-analysis ux_00000.h5 uy_00000.h5 uz_00000.h5 outfile.dat
!-------------------------------------------------------------------------------
subroutine field_analysis(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfile
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time, epsilon_loc, epsilon, fact, E, u_rms,lambda_macro,lambda_micro
  integer :: ix,iy,iz, mpicode
  real(kind=pr) :: Z_loc,Z_tot,nu2

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --field-analysis ux_00000.h5 uy_00000.h5 uz_00000.h5 outfile.dat"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "FLUSI field analysis. From a given vector field u (three files, one for each component)"
    write(*,*) "we'll compute the kinetic energy E=(ux^2 + uy^2 + uz^2)/2, then the curl vor = curl(u)"
    write(*,*) "and the enstrophy Z=(vorx^2 + vory^2 + vorz^2)/2 as well as the dissipation rate "
    write(*,*) "epsilon=nu*Z, where the viscosity is read from the files. (note you might have to add it"
    write(*,*) "using --set-hdf5-attribute to the file if it is missing) We print the output (integral "
    write(*,*) "and mean) to an ascii file given in the call, as well as on the screen"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,outfile)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if (mpirank==0) then
  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
    write (*,*) "Warning in arguments, files do not start with ux uy and uz"
    write (*,*) "note files have to be in the right order"
  endif
  endif

  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    call postprocessing_ascii_header(17)
    write(17,'(A)') "-----------------------------------"
    write(17,'(A)') "FLUSI field analysis"
    write(17,'(A)') "-----------------------------------"
  endif

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  ! field energy
  epsilon_loc = 0.5d0*sum(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2)*(dx*dy*dz)
  call MPI_REDUCE(epsilon_loc,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)

  ! compute vorticity
  if (mpirank==0) write(*,*) "Computing vorticity.."
  call fft3 (inx=u,outk=uk)
  call curl3_inplace (uk)
  call ifft3 (ink=uk,outx=u)

  ! compute enstrophy
  Z_loc = 0.5d0*sum(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2)*(dx*dy*dz)
  call MPI_REDUCE(Z_loc,Z_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)

  if(mpirank==0) then
    write(*,'("viscosity=",es15.8)') nu
    write(*,'("Total kinetic energy =",es15.8)') E
    write(*,'("Total Enstrophy =",es15.8)') Z_tot
    write(*,'("Total Dissipation =",es15.8)') nu*Z_tot
    write(*,'("Mean kinetic energy =",es15.8)') E / (xl*yl*zl)
    write(*,'("Mean Enstrophy =",es15.8)') Z_tot / (xl*yl*zl)
    write(*,'("Mean Dissipation =",es15.8)') nu*Z_tot / (xl*yl*zl)

    write(17,'("viscosity=",es15.8)') nu
    write(17,'("Total kinetic energy =",es15.8)') E
    write(17,'("Total Enstrophy =",es15.8)') Z_tot
    write(17,'("Total Dissipation =",es15.8)') nu*Z_tot
    write(17,'("Mean kinetic energy =",es15.8)') E / (xl*yl*zl)
    write(17,'("Mean Enstrophy =",es15.8)') Z_tot / (xl*yl*zl)
    write(17,'("Mean Dissipation =",es15.8)') nu*Z_tot / (xl*yl*zl)
  endif

  deallocate (u,uk)
  call fft_free()
end subroutine field_analysis
