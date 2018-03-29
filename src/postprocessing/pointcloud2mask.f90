subroutine pointcloud2mask(help)
  use vars
  use helpers
  use ini_files_parser_mpi
  use p3dfft_wrapper
  use insect_module
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: cloudfile, outfile, mode, dummy
  real(kind=pr),dimension(:,:),allocatable :: data_raw, points, normals
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, d0
  type(inifile) :: PARAMS
  integer :: ntri, matrixlines, matrixcols, safety


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pointcloud2mask pointcloud.txt outfile_000.h5 x0 y0 z0 xl yl zl nx ny nz d0"
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

  if (root) then
    write(*,'(80("~"))')
    write(*,*) "Reading cloud from file="//trim(adjustl(cloudfile))
    write(*,*) "Writing mask to file="//trim(adjustl(outfile))
    write(*,'("nx=",i5," ny=",i5," nz=",i5)') nx,ny,nz
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') xl,yl,zl
    write(*,'("x0=",g12.4," y0=",g12.4," z0=",g12.4)') x0,y0,z0
    write(*,'("d0=",g12.4)') d0
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

  if (root) then
    write(*,'(80("~"))')
    write(*,'("max_x=",g12.4," max_y=",g12.4," max_z=",g12.4)') maxval(points(:,1)), maxval(points(:,2)), maxval(points(:,3))
    write(*,'("min_x=",g12.4," min_y=",g12.4," min_z=",g12.4)') minval(points(:,1)), minval(points(:,2)), minval(points(:,3))
    write(*,'(80("~"))')

    if (maxval(points(:,1))>xl .or. minval(points(:,1))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (x)"
    if (maxval(points(:,2))>yl .or. minval(points(:,2))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (y)"
    if (maxval(points(:,3))>zl .or. minval(points(:,3))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (z)"
  endif

  deallocate( data_raw )

  !-----------------------------------------------------------------------------
  ! create mask
  !-----------------------------------------------------------------------------
  time = 0.d0
  ! the safety integer defines how large the vicinity of each point is, i.e. it
  ! is mainly a computational cost parameter
  safety = 5
  ! we don't use ghost points right now, but maybe in the future if we want to fill
  ! the body interior with some algorithm
  ng = 0

  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()


  ! smoothing for mask function:
  if (nx/=1) then
    smoothing = 1.5d0*max(dz,dy,dx)
  else
    smoothing = 1.5d0*max(dz,dy)
  endif

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  work = 0.0d0

  if (root) write(*,*) "...computing distance only near fluid-solid interface"
  call mask_from_pointcloud(matrixlines, points, normals, work, safety, d0, int(0,kind=2) )


  ! save result to file
  call save_field_hdf5 ( time, outfile, work )

  ! finalize
  deallocate (work)
  call fft_free()
end subroutine


subroutine mask_from_pointcloud(N, points, normals, work, safety, d0, color)
  use vars
  use penalization
  use helpers
  use insect_module, only : steps
  use stl_file_reader ! for the function distance to triangle
  implicit none
  integer, intent(in) :: N,safety
  integer(kind=2), intent(in) :: color
  real(kind=pr),dimension(1:N,1:3),intent(in):: points
  real(kind=pr),dimension(1:N,1:3),intent(in):: normals
  real(kind=pr),intent(in):: d0
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix,iy,iz,i,xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=pr) :: x,y,z,tmp, s

  work = 9.0d7

  if (maxval(points(:,1))>xl .or. minval(points(:,1))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (x)"
  if (maxval(points(:,2))>yl .or. minval(points(:,2))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (y)"
  if (maxval(points(:,3))>zl .or. minval(points(:,3))<0.0d0) write(*,*) "WARNING: some points may be outside of the domain! (z)"


  do i = 1, N
    ! bounding box of the vicinity of the langrangian makert point of the cloud.
    xmin = nint(points(i,1)/dx)-safety
    xmax = nint(points(i,1)/dx)+safety

    ymin = nint(points(i,2)/dy)-safety
    ymax = nint(points(i,2)/dy)+safety

    zmin = nint(points(i,3)/dz)-safety
    zmax = nint(points(i,3)/dz)+safety

    do iz = zmin, zmax
      do iy = ymin, ymax
        do ix = xmin, xmax
          if ( on_proc((/ix,iy,iz/)) ) then
            x = dx*dble(ix)
            y = dy*dble(iy)
            z = dz*dble(iz)

            ! the distance to the current point:
            ! note this is the square of the distance (cheaper to take sqrt later!)
            tmp = (x-points(i,1))**2 + (y-points(i,2))**2 + (z-points(i,3))**2

            ! if closer (in abs value!) then use this now
            if ( dabs(tmp) < dabs(work(ix,iy,iz) ) ) then
                ! take care of sign: exterior / interior
                if ((x-points(i,1))*normals(i,1) + (y-points(i,2))*normals(i,2) + (z-points(i,3))*normals(i,3) < 0.0) then
                  tmp = -tmp
                endif
                work(ix,iy,iz)  = tmp
            endif

          endif
        enddo
      enddo
    enddo

  enddo

  !-----------------------------------------------------------------------------
  ! convert signed distance function to mask function chi
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
          tmp = work(ix,iy,iz)
          if (tmp >= 0.0d0) then
              tmp = sqrt(tmp)
          else
              tmp = -sqrt(-tmp)
          endif
        work(ix,iy,iz) = steps( tmp-d0, 0.d0 )
      enddo
    enddo
  enddo

  if (allocated(mask_color)) then
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          if (work(ix,iy,iz)>0.0) mask_color(ix,iy,iz)=color
        enddo
      enddo
    enddo
  endif

end subroutine mask_from_pointcloud
