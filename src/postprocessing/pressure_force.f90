subroutine pressure_force(help)
  use vars
  use p3dfft_wrapper
  use module_helpers
  use basic_operators
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_p, fname_mask, outfile
  complex(kind=pr),dimension(:,:,:),allocatable :: pk
  complex(kind=pr),dimension(:,:,:,:),allocatable :: workc
  real(kind=pr),dimension(:,:,:),allocatable :: p, chi
  real(kind=pr),dimension(:,:,:),allocatable :: workr
  real(kind=pr) :: time, divu_max
  real(kind=pr),dimension(1:3) :: force
  integer :: i


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pressure-force p_0000.h5 mask_0000.h5 outfile.txt"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Compute the integral force contribution of the pressure. "
    write(*,*) " "
    write(*,*) " Note that F = \int ( sigma.dA) = \int (div(sigma)) dV"
    write(*,*) " therefore, we can set sigma = -p*delta_ij and we are left with integration of grad(p)"
    write(*,*) " over the volume of the obstacle."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_p)
  call get_command_argument(4,fname_mask)
  call get_command_argument(5,outfile)


  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "Compute pressure force acting on obstacle"
    write(*,*) "Pressure and mask files:"
    write(*,*) trim(adjustl(fname_p))
    write(*,*) trim(adjustl(fname_mask))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(outfile))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_p )
  call check_file_exists( fname_mask )

  call fetch_attributes( fname_p, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  allocate(p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(chi(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  call read_single_file ( fname_p, p )
  call read_single_file ( fname_mask, chi )

  ! to Fourier space
  call fft (inx=p, outk=pk)

  deallocate( p )
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3) )
  allocate(workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call gradient(pk,workc(:,:,:,1),workc(:,:,:,2),workc(:,:,:,3))

  do i=1,3
    ! to x-space
    call ifft( ink=workc(:,:,:,i), outx=workr )
    ! apply mask
    workr = workr * chi
    ! sum gives the pressure force (volume integral inside the obstacle)
    force(i) = -1.d0*mpisum(sum(workr))*dx*dy*dz
  end do

  if (mpirank==0) then
    open  (14, file = outfile, status = 'replace')
    write(14,'(A)') "% pressure force computed from "//trim(adjustl(fname_p))//" "//trim(adjustl(fname_mask))
    write (14,'(4(es15.8,1x))') time, force
    write (* ,'(4(es15.8,1x))') time, force
    close (14)
  end if

  deallocate (chi,pk,workc,workr)
  call fft_free()

end subroutine pressure_force
