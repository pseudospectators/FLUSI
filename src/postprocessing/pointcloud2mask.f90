subroutine pointcloud2mask(help)
  use vars
  use module_helpers
  use module_ini_files_parser_mpi
  use p3dfft_wrapper
  use penalization, only : mask_color
  use module_insects
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: cloudfile, outfile, mode, dummy
  real(kind=pr),dimension(:,:),allocatable :: data_raw, points, normals
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, d0, smoothing
  type(inifile) :: PARAMS
  integer :: ntri, matrixlines, matrixcols, safety
  type(diptera) :: Insect

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pointcloud2mask pointcloud.txt outfile_000.h5 x0 y0 z0 xl yl zl nx ny nz d0 [--minimal-domain]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! read an point cloud file, which basically contains a list of points and their normals,"
    write(*,*) "! and compute the mask function for it. "
    write(*,*) "! "
    write(*,*) "! The ascii point cloud file is expected to have the following columns:"
    write(*,*) "! x, y, z, .... (unused) ...., nx, ny ,nz"
    write(*,*) "! The first line is assumed to be the header and skipped. At least 6 columns are required."
    write(*,*) "! Additional columns in the middle are ignored (maybe RGB values etc). Normals are used for"
    write(*,*) "! defining interior/exterior."
    write(*,*) "! "
    write(*,*) "! The origin is specified as coordinates x0, y0, z0 and added to all points in the cloud. It"
    write(*,*) "! serves thus to translate the cloud (x_cloud + x0)"
    write(*,*) "! "
    write(*,*) "! The output hdf5 field has box size xl,yl,zl and resolution nx,ny,nz."
    write(*,*) "! "
    write(*,*) "! Usually, the points are interpreted as points ON the interface, i.e. the distance on the grid"
    write(*,*) "! to a point is d = sqrt( (x-xp)**2+(y-yp)**2+(z-zp)**2) but it may be handy to make the mask thicker"
    write(*,*) "! i.e. compute d = d-d0 and then compute the mask function mask(d). This was used in uCT data handling"
    write(*,*) "! in order to use an existing segmentation as a mask for the original data."
    write(*,*) "! "
    write(*,*) "! "
    write(*,*) "! "
    write(*,*) "! NOTE: we compute the mask ONLY near the surface, NOT in the interior of the body."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,cloudfile)
  call get_command_argument(4,outfile)
  call get_command_argument(5,dummy)
  read(dummy,*) x0
  call get_command_argument(6,dummy)
  read(dummy,*) y0
  call get_command_argument(7,dummy)
  read(dummy,*) z0
  call get_command_argument(8,dummy)
  read(dummy,*) xl
  call get_command_argument(9,dummy)
  read(dummy,*) yl
  call get_command_argument(10,dummy)
  read(dummy,*) zl
  call get_command_argument(11,dummy)
  read(dummy,*) nx
  call get_command_argument(12,dummy)
  read(dummy,*) ny
  call get_command_argument(13,dummy)
  read(dummy,*) nz
  call get_command_argument(14,dummy)
  read(dummy,*) d0
  call get_command_argument(15,mode)

  if (root) then
    write(*,'(80("~"))')
    write(*,*) "Reading cloud from file="//trim(adjustl(cloudfile))
    write(*,*) "Writing mask to file="//trim(adjustl(outfile))
    write(*,'("nx=",i5," ny=",i5," nz=",i5)') nx,ny,nz
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') xl,yl,zl
    write(*,'("x0=",g12.4," y0=",g12.4," z0=",g12.4)') x0,y0,z0
    write(*,'("d0=",g12.4)') d0
    write(*,'("mode=",A)') mode
    write(*,'(80("~"))')
  endif

  !-----------------------------------------------------------------------------
  ! read point cloud data, note all mpiranks hold a copy of the entire array.
  !-----------------------------------------------------------------------------
  call check_file_exists(cloudfile)

  call count_lines_in_ascii_file_mpi(cloudfile, matrixlines, 1)
  call count_cols_in_ascii_file_mpi(cloudfile, matrixcols, 1)

  allocate( data_raw(matrixlines, matrixcols))
  allocate( normals(matrixlines,1:3) )
  allocate( points(matrixlines,1:3) )

  call read_array_from_ascii_file_mpi(cloudfile, data_raw, 1)

  ! we assume the following columns:
  ! x, y, z, .... (unused) ...., nx, ny ,nz

  normals = data_raw(:,matrixcols-2:matrixcols)
  points(:,1) = data_raw(:,1) + x0
  points(:,2) = data_raw(:,2) + y0
  points(:,3) = data_raw(:,3) + z0
  deallocate( data_raw )

  ! the safety integer defines how large the vicinity of each point is, i.e. it
  ! is mainly a computational cost parameter
  safety = 4

  ! in most cases, the origin of the grid is zero
  origin = 0.0d0

  if (mode == "--minimal-domain") then

      ! There is several solutions now:
      ! - keep the resolution specified in the call and adjust nx,ny,nz
      ! - keep the number of points in the call and end up with dx /= dy /= dz
      ! - keep the number of points the same in the best sampled direction, then adjust the other two resolutions accordingly.
      ! right now, we're doing (c)
    if (root) write(*,*) "As you set --minimal-domain, I will ignore xl,yl,zl and create"
    if (root) write(*,*) "the smallest possible cube that contains your pointcloud. I will use"
    if (root) write(*,*) "the same resolution you specified in the call."
    ! the resolution in the call was:
    dx = xl / dble(nx)

    origin = (/minval(points(:,1)),minval(points(:,2)),minval(points(:,3))/) - dble(safety)*dx
    xl = maxval(points(:,1))-origin(1) + dble(safety)*dx
    yl = maxval(points(:,2))-origin(2) + dble(safety)*dx
    zl = maxval(points(:,3))-origin(3) + dble(safety)*dx

    dx = minval( (/xl/dble(nx), yl/dble(ny), zl/dble(nz)/) )

    nx = nint(xl/dx)
    ny = nint(yl/dx)
    nz = nint(zl/dx)

    if (root) write(*,'("nx=",i5," ny=",i5," nz=",i5)') nx,ny,nz
    if (root) write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') xl,yl,zl
    if (root) write(*,'("x0=",g12.4," y0=",g12.4," z0=",g12.4)') origin

  endif

  ! check if the grid is large enough for the data
  if (root) then
    write(*,'(80("~"))')
    write(*,'("max_x=",g12.4," max_y=",g12.4," max_z=",g12.4)') maxval(points(:,1)), maxval(points(:,2)), maxval(points(:,3))
    write(*,'("min_x=",g12.4," min_y=",g12.4," min_z=",g12.4)') minval(points(:,1)), minval(points(:,2)), minval(points(:,3))
    write(*,'(80("~"))')

    if (maxval(points(:,1))>xl+origin(1) .or. minval(points(:,1))<origin(1)) then
      write(*,*) "WARNING: some points may be outside of the domain! (x)"
    endif
    if (maxval(points(:,2))>yl+origin(2) .or. minval(points(:,2))<origin(2)) then
      write(*,*) "WARNING: some points may be outside of the domain! (y)"
    endif
    if (maxval(points(:,3))>zl+origin(3) .or. minval(points(:,3))<origin(3)) then
      write(*,*) "WARNING: some points may be outside of the domain! (z)"
    endif
  endif



  !-----------------------------------------------------------------------------
  ! create mask
  !-----------------------------------------------------------------------------
  time = 0.d0
  ! we don't use ghost points right now, but maybe in the future if we want to fill
  ! the body interior with some algorithm
  ng = 0
  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()

  call insect_init( 0.d0, "params_template_fsi.ini", Insect, .false., "", (/xl,yl,zl/), nu, dx, periodic=periodic)



  ! smoothing for mask function:
  if (nx/=1) then
    smoothing = 1.5d0*max(dz,dy,dx)
  else
    smoothing = 1.5d0*max(dz,dy)
  endif

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  work = 0.0d0
  mask_color = int(0,kind=2)

  if (root) write(*,*) "...computing distance only near fluid-solid interface"
  call mask_from_pointcloud(points, normals, origin+(/dble(ra(1))*dx, dble(ra(2))*dy, dble(ra(3))*dz/), &
                            (/dx,dy,dz/), work, safety, smoothing, d0, mask_color, 0_2)

  ! save result to file
  call save_field_hdf5 ( time, outfile, work )

  ! finalize
  deallocate (work, mask_color)
  call fft_free()
end subroutine
