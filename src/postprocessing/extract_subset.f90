!-------------------------------------------------------------------------------
! extract subset
! loads a file to memory and extracts a subset, writing to a different file. We
! assume here that you do this for visualization; in this case, one usually keeps
! the original, larger files. For simplicity, ensure all files follow FLUSI
! naming convention.
!---
! Note: through using module_helpers.f90::get_dsetname, it is finally possible to write
! to a subfolder *.h5 file directly from flusi
! (Thomas, 03/2015)
!---
! Using HDF5s hyperslab functions, we can read only a specific part into the
! memory - at no point we have to load the entire original file before we can
! downsample it. This is a good step forward. (Thomas 03/2015)
!-------------------------------------------------------------------------------
! call:
! ./flusi --postprocess --extract-subset ux_00000.h5 sux_00000.h5 128:1:256 128:2:1024 1:1:9999
!-------------------------------------------------------------------------------
subroutine extract_subset(help)
  use mpi
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
  integer :: nx1,nx2, ny1,ny2, nz1,nz2, nxs,nys,nzs
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
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --extract-subset ux_00000.h5 sux_00000.h5 128:1:256 128:2:1024 1:1:9999"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! extract subset"
    write(*,*) "! loads a file to memory and extracts a subset, writing to a different file. We"
    write(*,*) "! assume here that you do this for visualization; in this case, one usually keeps"
    write(*,*) "! the original, larger files. For simplicity, ensure all files follow FLUSI"
    write(*,*) "! naming convention."
    write(*,*) "!---"
    write(*,*) "! Note: through using module_helpers.f90::get_dsetname, it is finally possible to write"
    write(*,*) "! to a subfolder *.h5 file directly from flusi"
    write(*,*) "! (Thomas, 03/2015)"
    write(*,*) "!---"
    write(*,*) "! Using HDF5s hyperslab functions, we can read only a specific part into the"
    write(*,*) "! memory - at no point we have to load the entire original file before we can"
    write(*,*) "! downsample it. This is a good step forward. (Thomas 03/2015)"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: Nope"
    return
  endif


  if (mpisize/=1) then
    call abort(111, "./flusi --postprocess --extract-subset is a SERIAL routine, use 1CPU only")
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  ! get filename to save subset to
  call get_command_argument(4,fname_out)

  dsetname_in  = get_dsetname( fname_in )
  dsetname_out = get_dsetname( fname_out )

  write(*,'("dsetname=",A,1x,A)') trim(adjustl(dsetname_in)),trim(adjustl(dsetname_out))

  ! fetch attributes from source file
  call fetch_attributes( fname_in, nx, ny, nz, xl, yl, zl, time, nu, origin )

  call get_command_argument(5,xset)
  call get_command_argument(6,yset)
  call get_command_argument(7,zset)

  ! red in subset from command line. it is given in the form
  ! ixmin:xspacing:ixmax as a string.
  read (xset(1:index(xset,':')-1) ,*) nx1
  read (xset(index(xset,':',.true.)+1:len_trim(xset)),*) nx2
  read (xset(index(xset,':')+1:index(xset,':',.true.)-1),*) nxs

  read (yset(1:index(yset,':')-1) ,*) ny1
  read (yset(index(yset,':',.true.)+1:len_trim(yset)),*) ny2
  read (yset(index(yset,':')+1:index(yset,':',.true.)-1),*) nys

  read (zset(1:index(zset,':')-1) ,*) nz1
  read (zset(index(zset,':',.true.)+1:len_trim(zset)),*) nz2
  read (zset(index(zset,':')+1:index(zset,':',.true.)-1),*) nzs


  ! stop if subset exceeds array bounds
  if ( nx1<0 .or. nx2>nx-1 .or. ny1<0 .or. ny2>ny-1 .or. nz1<0 .or. nz2>nz-1) then
    write (*,*) "subset indices exceed array bounds....proceed, but correct mistake"
    nx1 = max(nx1,0)
    ny1 = max(ny1,0)
    nz1 = max(nz1,0)
    nx2 = min(nx-1,nx2)
    ny2 = min(ny-1,ny2)
    nz2 = min(nz-1,nz2)
  endif


  write(*,'("Cropping field from " &
  &,"0:",i4," | 0:",i4," | 0:",i4,&
  &"   to subset   "&
  &,i4,":",i2,":",i4," | ",i4,":",i2,":",i4," | ",i4,":",i2,":",i4)')&
  nx-1,ny-1,nz-1,nx1,nxs,nx2,ny1,nys,ny2,nz1,nzs,nz2



  !-----------------------------------------------------------------------------
  ! compute dimensions of reduced subset:
  nx_red = nx1 + floor( dble(nx2-nx1)/dble(nxs) )  - nx1 + 1
  ny_red = ny1 + floor( dble(ny2-ny1)/dble(nys) )  - ny1 + 1
  nz_red = nz1 + floor( dble(nz2-nz1)/dble(nzs) )  - nz1 + 1
  write (*,'("Size of subset is ",3(i4,1x))') nx_red, ny_red, nz_red
  !-----------------------------------------------------------------------------

  ! we figured out how big the subset array is
  allocate ( field(0:nx_red-1,0:ny_red-1,0:nz_red-1) )


  call Fetch_attributes( fname_in, nx_file,ny_file,nz_file,xl_file,yl_file,zl_file,time, nu, origin )

  ! define new origin of grid
  origin = (/ dble(nx1)*xl/dble(nx), dble(ny1)*yl/dble(ny), dble(nz1)*zl/dble(nz) /)
  write (*,'("New origin is ",3(g12.4,1x))') origin

  !-----------------------------------------------------------------------------
  ! load the file
  ! the basic idea is to just allocate the smaller field, and use hdf5
  ! to read just this field from the input file.
  !-----------------------------------------------------------------------------
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call h5fopen_f (fname_in, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Definition of memory distribution
  dimensions_file  = (/ nx_file, ny_file, nz_file/)
  dimensions_local = (/ nx_red , ny_red,  nz_red /)
  offset = (/ nx1, ny1, nz1 /)
  stride = (/ nxs, nys, nzs /)
  chunking_dims = 1 !min(nx_red,ny_red,nz_red)

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname_in, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, dimensions_local, &
  error, stride, count)


  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  ! actual read is the next command:
  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
  mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  ! check if we loaded crap
  call checknan(field,"recently loaded field")

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error)
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id,error)
  call H5close_f(error)

  ! set up dimensions in global variables, since save_field_hdf5 relies on this
  ra = 0
  rb(1) = nx_red-1
  rb(2) = ny_red-1
  rb(3) = nz_red-1
  nx = nx_red
  ny = ny_red
  nz = nz_red
  dx = xl_file / dble(nx_file)
  dy = yl_file / dble(ny_file)
  dz = zl_file / dble(nz_file)
  xl = dx + dble(nx1+(nx_red-1)*nxs)*dx - dble(nx1)*dx
  yl = dy + dble(ny1+(ny_red-1)*nys)*dy - dble(ny1)*dy
  zl = dz + dble(nz1+(nz_red-1)*nzs)*dz - dble(nz1)*dz

  if ( dble(nx)*dble(ny)*dble(nz) > 500.0d6) then
    write(*,'(A)') "WARNING! The subset you ordered is rather larger (>500M points). The hdf wrapper"
    write(*,'(A)') "might crash during writing. If that happens, change max_chunk to say 256 in hdf5_wrapper.f90"
    write(*,'(A)') "then recompile and retry"
  endif

  ! Done! Write extracted subset to disk and be happy with the result
  call save_field_hdf5 ( time, fname_out, field )
  call write_attrib_dble(fname_out, get_dsetname(fname_out), "origin", origin)

end subroutine extract_subset
