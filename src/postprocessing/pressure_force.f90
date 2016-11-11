subroutine pressure_force(help)
  use vars
  use p3dfft_wrapper
  use helpers
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


  call get_command_argument(3,fname_p)
  call get_command_argument(4,fname_mask)
  call get_command_argument(5,outfile)


  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "Compute pressure force acting on obstacle"
    write(*,*) "Pressure and mask files:"
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname_p))
    write(*,*) trim(adjustl(fname_mask))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(outfile))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_p )
  call check_file_exists( fname_mask )
  call check_file_exists( outfile )

  call fetch_attributes( fname_p, nx, ny, nz, xl, yl, zl, time, nu )

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
    force(i) = sum(workr)*dx*dy*dz
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
