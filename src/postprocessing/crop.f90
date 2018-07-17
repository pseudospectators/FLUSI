subroutine crop(help)
  use p3dfft_wrapper
  use vars
  use hdf5
  use hdf5_wrapper
  use basic_operators
  use module_helpers

  implicit none

  logical, intent(in) :: help
  character(len=strlen) :: fname_in, fname_out, dsetname_in, dsetname_out
  character(len=strlen) :: xset,yset,zset
  integer :: ix,iy,iz,i
  ! reduced domain size
  integer :: nx1,nx2, ny1,ny2, nz1,nz2
  ! sizes of the new array
  integer :: nx_red, ny_red, nz_red, ix_red, iy_red, iz_red
  ! reduced domain extends
  real(kind=pr) :: xl1, yl1, zl1
  real(kind=pr), dimension(:,:,:), allocatable :: field

  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: nx_file, ny_file, nz_file, mpierror

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1


  if (help.and.root) then
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "./flusi -p --crop ux_00000.h5 sux_00000.h5 0:256 0:128 513:1024"
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "! crop field."
    write(*,'(A)') "! extracts a subset of a file to another file. Of course, you can also overwrite the source file."
    write(*,'(A)') "!---"
    write(*,'(A)') "! Note: this routine is Parallel but does not allow striding. Use --extract-subset"
    write(*,'(A)') "! for a serial version that does support striding."
    write(*,'(A)') "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A)') "Parallel: yes"
    return
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  ! get filename to save subset to
  call get_command_argument(4,fname_out)

  dsetname_in  = get_dsetname( fname_in )
  dsetname_out = get_dsetname( fname_out )

  ! fetch attributes from source file
  call Fetch_attributes( fname_in, nx_file, ny_file, nz_file, xl_file, yl_file, zl_file, time, nu, origin )

  call get_command_argument(5,xset)
  call get_command_argument(6,yset)
  call get_command_argument(7,zset)

  ! red in subset from command line. it is given in the form
  ! ixmin:xspacing:ixmax as a string.
  read (xset(1:index(xset,':')-1) ,*) nx1
  read (xset(index(xset,':',.true.)+1:len_trim(xset)),*) nx2

  read (yset(1:index(yset,':')-1) ,*) ny1
  read (yset(index(yset,':',.true.)+1:len_trim(yset)),*) ny2

  read (zset(1:index(zset,':')-1) ,*) nz1
  read (zset(index(zset,':',.true.)+1:len_trim(zset)),*) nz2

  if ( nx1<0 .or. nx2>nx_file-1 .or. ny1<0 .or. ny2>ny_file-1 .or. nz1<0 .or. nz2>nz_file-1) then
    if (root) write (*,*) "subset indices exceed array bounds....proceed, but correct mistake"
    nx1 = max(nx1,0)
    ny1 = max(ny1,0)
    nz1 = max(nz1,0)
    nx2 = min(nx_file-1,nx2)
    ny2 = min(ny_file-1,ny2)
    nz2 = min(nz_file-1,nz2)
  endif

  if (root) then
    write(*,'("Cropping field from 0:",i4," | 0:",i4," | 0:",i4,"   to subset   "&
    &,i4,":",i4," | ",i4,":",i4," | ",i4,":",i4)')&
    nx_file-1,ny_file-1,nz_file-1,nx1,nx2,ny1,ny2,nz1,nz2
  endif



  !-----------------------------------------------------------------------------
  ! compute dimensions of reduced subset:
  nx_red = nx2-nx1 + 1
  ny_red = ny2-ny1 + 1
  nz_red = nz2-nz1 + 1
  if (root) write (*,'("Size of subset is ",3(i4,1x))') nx_red, ny_red, nz_red
  !-----------------------------------------------------------------------------

  ! define new origin of grid
  origin = (/ dble(nx1)*xl_file/dble(nx_file), dble(ny1)*yl_file/dble(ny_file), dble(nz1)*zl_file/dble(nz_file) /)
  if (root) write (*,'("New origin is ",3(g12.4,1x))') origin

  ! set up dimensions in global variables, since save_field_hdf5 relies on this
  nx = nx_red
  ny = ny_red
  nz = nz_red
  dx = xl_file / dble(nx_file)
  dy = yl_file / dble(ny_file)
  dz = zl_file / dble(nz_file)
  xl = dx + dble(nx1+(nx_red-1))*dx - dble(nx1)*dx
  yl = dy + dble(ny1+(ny_red-1))*dy - dble(ny1)*dy
  zl = dz + dble(nz1+(nz_red-1))*dz - dble(nz1)*dz

  call decomposition_initialize()
  ra = ra + (/nx1,ny1,nz1/)
  rb = rb + (/nx1,ny1,nz1/)
  allocate(field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  ! NOTE: the hdf5 wrapper actually reads a subset of the total data for each mpirank
  ! so it was actually very easy to adapt
  call read_field_hdf5( fname_in, get_dsetname(fname_in), ra, rb, field )

  ra = ra - (/nx1,ny1,nz1/)
  rb = rb - (/nx1,ny1,nz1/)

  ! Done! Write extracted subset to disk and be happy with the result
  call save_field_hdf5 ( time, fname_out, field )
  call write_attrib_dble(fname_out, get_dsetname(fname_out), "origin", origin)
end subroutine crop
