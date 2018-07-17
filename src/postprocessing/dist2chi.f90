!-------------------------------------------------------------------------------
! ./flusi -p --energy ux_00.h5 uy_00.h5 uz_00.h5 outfile_00.h5
!-------------------------------------------------------------------------------
! load the vector components from file and compute & save the energy to
! another HDF5 file. (energy = (ux^2+uy^2+uz^2)/2)
subroutine dist2chi(help)
  use vars
  use basic_operators
  use p3dfft_wrapper
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: distfile,maskfile,dummy
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, sign, tmp
  integer :: ix,iy,iz

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --dist2chi dist_000.h5 mask_000.h5 1"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! read in signed distance and convert it to a mask function suitable for simulations"
    write(*,*) "! Attention: it depends on the sign what is considered outside and what is inside"
    write(*,*) "! in some cases, the sign has to be inverted, call with -1 then"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,distfile)
  call get_command_argument(4,maskfile)
  call get_command_argument(5,dummy)
  read(dummy,*) sign

  if (root) then
    write(*,*) "reading signed distance from file="//trim(adjustl(distfile))
    write(*,*) "writing mask function to file="//trim(adjustl(maskfile))
    write(*,*) "sign is=", sign
  endif

  call check_file_exists( distfile )
  call fetch_attributes( distfile, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file ( distfile, work )

  !-----------------------------------------------------------------------------
  ! convert signed distance function to mask function chi
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        call smoothstep(tmp, sign*work(ix,iy,iz), 0.d0, 2.d0*dz)
        work(ix,iy,iz) = tmp !steps( sign*work(ix,iy,iz), 0.d0 )
      enddo
    enddo
  enddo

  ! save output
  call save_field_hdf5 ( time, maskfile, work )

  deallocate (work)
  call fft_free()

end subroutine dist2chi
