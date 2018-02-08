subroutine pointcloud2mask(help)
  use vars
  use helpers
  use ini_files_parser_mpi
  use p3dfft_wrapper
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: stlfile, outfile, mode, dummy, flag
  real(kind=pr),dimension(:,:),allocatable :: data, points, normals
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, scale
  type(inifile) :: PARAMS
  integer :: ntri, matrixlines, matrixcols, safety


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pointcloud2mask infile.stl outfile_000.h5 NORMALIZATION x0 y0 z0 xl yl zl nx ny nz FLAG safety"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! read an *.stl file and compute the signed distance function for it"
    write(*,*) "! the NORMALIZATION flag can be either --lx, --ly, --lz or --scale=123.0"
    write(*,*) "! in the latter case, you give the normalization scale directly"
    write(*,*) "! "
    write(*,*) "! origin is specified as coordinates"
    write(*,*) "! output hdf5 field has box size xl,yl,zl and resolution nx,ny,nz"
    write(*,*) "! "
    write(*,*) "! Note *.stl files do not contain dimensions"
    write(*,*) "! "
    write(*,*) "! FLAG can be --everywhere or --surface. in the latter last flag safety is nearness to interface"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,stlfile)
  call get_command_argument(4,outfile)
  call get_command_argument(5,mode)
  call get_command_argument(6,dummy)
  read(dummy,*) x0
  call get_command_argument(7,dummy)
  read(dummy,*) y0
  call get_command_argument(8,dummy)
  read(dummy,*) z0
  call get_command_argument(9,dummy)
  read(dummy,*) xl
  call get_command_argument(10,dummy)
  read(dummy,*) yl
  call get_command_argument(11,dummy)
  read(dummy,*) zl
  call get_command_argument(12,dummy)
  read(dummy,*) nx
  call get_command_argument(13,dummy)
  read(dummy,*) ny
  call get_command_argument(14,dummy)
  read(dummy,*) nz
  call get_command_argument(15,flag)

  if (mode(1:8)=="--scale=") then
    read(mode(9:80),*) scale
    mode = "--scale"
  endif

  if (root) then
    write(*,'(80("~"))')
    write(*,*) "Reading stl file="//trim(adjustl(stlfile))
    write(*,*) "Writing signe distance to file="//trim(adjustl(outfile))
    write(*,'("nx=",i5," ny=",i5," nz=",i5)') nx,ny,nz
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') xl,yl,zl
    write(*,'("x0=",g12.4," y0=",g12.4," z0=",g12.4)') x0,y0,z0
    write(*,'(80("~"))')
  endif

  ! read stl file and normalize it according to the criteria specified in the
  ! command line call
  call check_file_exists(stlfile)

  call count_lines_in_ascii_file_mpi(stlfile, matrixlines, 1)
  call count_cols_in_ascii_file_mpi(stlfile, matrixcols, 1)

  allocate( data(matrixlines, matrixcols))
  allocate( normals(matrixlines,1:3) )
  allocate( points(matrixlines,1:3) )

  call read_array_from_ascii_file_mpi(stlfile, data, 1)

  normals = data(:,matrixcols-2:matrixcols)
  points(:,1) = data(:,1)-x0+xl/2.0d0
  points(:,2) = data(:,2)-y0+yl/2.0d0
  points(:,3) = data(:,3)-z0+zl/2.0d0

  deallocate( data )

  time=0.d0
  ! initialize code and domain decomposition, but do not use FFTs
  safety = 5
  ng = 0
  call decomposition_initialize(.false.)
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  work = 0.0d0

  if (root) write(*,*) "...computing distance only near fluid-solid interface"
  call signed_distance_from_pointcloud(matrixlines, points, normals, work, safety)


  ! save result to file
  call save_field_hdf5 ( time, outfile, work )

  ! finalize
  deallocate (work)
  call fft_free()
end subroutine


subroutine signed_distance_from_pointcloud(N, points, normals, work, safety)
  use vars
  use helpers
  use insect_module, only : steps
  use stl_file_reader ! for the function distance to triangle
  implicit none
  integer, intent(in) :: N,safety
  real(kind=pr),dimension(1:N,1:3),intent(in):: points
  real(kind=pr),dimension(1:N,1:3),intent(in):: normals
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix,iy,iz,i,xmin,xmax,ymin,ymax,zmin,zmax
  real(kind=pr) :: x,y,z,tmp, s

  work = 9.0d4

  do i = 1, N
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
            tmp = dsqrt( (x-points(i,1))**2 + (y-points(i,2))**2 + (z-points(i,3))**2 )
            if ((x-points(i,1))*normals(i,1) + (y-points(i,2))*normals(i,2) + (z-points(i,3))*normals(i,3) < 0.0) then
              tmp = -tmp
            endif
            ! if closer (in abs value!) then use this now
            if ( dabs(tmp) < dabs(work(ix,iy,iz) ) ) then
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
        work(ix,iy,iz) = steps( work(ix,iy,iz)-0.0*dx, 0.d0 )
      enddo
    enddo
  enddo

end subroutine signed_distance_from_pointcloud
