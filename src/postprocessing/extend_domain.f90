!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine extend_domain(help)
  use vars
  use p3dfft_wrapper
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_in, fname_out, tmp
  integer :: nx_new, ny_new, nz_new
  integer :: nx_org, ny_org, nz_org
  real(kind=pr) :: time
  real(kind=pr),dimension(:,:,:),allocatable :: u_org, u_new

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --extend-domain source.h5 target.h5 256 256 512"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Load an existing field with given resolution, then copy its content to the bottom lower"
    write(*,*) "corner of a larger field with the given size. We keep the same resolution: that means,"
    write(*,*) "we are not doing upsampling, but we extend the domain."
    write(*,*) ""
    write(*,*) "As a mpi-decomposed field cannot simply be copied locally, but instead involves communication"
    write(*,*) "this routine has to be run in serial, unfortunately. at least, it does not allocate any extra memory."
    write(*,*) ""
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: nope"
    return
  endif

  if (mpisize/=1) then
    write(*,*) "./flusi --postprocess --extend-domain is a SERIAL routine, use 1CPU only"
    call abort()
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  call get_command_argument(4,fname_out)

  ! read target resolution from command line
  call get_command_argument(5,tmp)
  read (tmp,*) nx_new
  call get_command_argument(6,tmp)
  read (tmp,*) ny_new
  call get_command_argument(7,tmp)
  read (tmp,*) nz_new

  call fetch_attributes( fname_in, nx_org, ny_org, nz_org, xl, yl, zl, time, nu )

  write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  write(*,*) "reading file "//trim(adjustl(fname_in))//" and writing to "//trim(adjustl(fname_out))
  write(*,'("Target field size= ",3(i4,1x))') nx_new, ny_new, nz_new
  write(*,'("Source field size= ",3(i4,1x))') nx_org, ny_org, nz_org
  write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  if ((nx_new<nx_org).or.(ny_new<ny_org).or.(nz_new<nz_org)) then
    call abort("new resolution may not be smaller than old resolution")
  endif

  write(*,*) "Initializing SMALL field..."
  nx = nx_org
  ny = ny_org
  nz = nz_org
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  ra = (/0,0,0/)
  rb = (/nx-1,ny-1,nz-1/)
  !-----------------------

  allocate(u_org(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call read_single_file(fname_in,u_org)

  !-----------------------
  write(*,*) "Initializing BIG field..."
  nx = nx_new
  ny = ny_new
  nz = nz_new
  xl = dx*dble(nx_new)
  yl = dy*dble(ny_new)
  zl = dz*dble(nz_new)
  ra = (/0,0,0/)
  rb = (/nx-1,ny-1,nz-1/)

  !-----------------------
  allocate(u_new(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  write(*,*) "data allocated. we'll now copy the data from small to big field"
  ! copy data to lower bottom corner
  u_new(0:nx_org-1, 0:ny_org-1, 0:nz_org-1) = u_org
  deallocate( u_org )

  write(*,*) "Saving upsampled field to " // trim(adjustl(fname_out))
  call save_field_hdf5(time,fname_out,u_new)

  deallocate( u_new )
end subroutine extend_domain
