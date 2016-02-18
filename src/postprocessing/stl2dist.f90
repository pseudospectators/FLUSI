subroutine stl2dist(help)
  use vars
  use helpers
  use stl_file_reader
  use p3dfft_wrapper
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: stlfile, outfile, mode, dummy, flag
  real(kind=pr),dimension(:,:),allocatable :: triangles,normals
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, scale
  integer :: ntri


  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --stl2dist infile.stl outfile_000.h5 NORMALIZATION x0 y0 z0 xl yl zl nx ny nz FLAG safety"
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
  call read_stl_file(stlfile, ntri, triangles, normals)
  call normalize_stl_file(ntri, triangles, mode, scale)

  time=0.d0
  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  call fft_initialize() ! also initializes the domain decomp

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  if (trim(adjustl(flag))=="--surface") then
    if (root) write(*,*) "...computing distance only near fluid-solid interface"
    ! compute signed distance function on all grid points
    call signed_distance_from_triangles_surface(ntri, triangles, normals, work)

  elseif (trim(adjustl(flag))=="--everywhere") then
    if (root) write(*,*) "...computing distance everywhere in eulerian grid"
    ! compute signed distance function on all grid points
    call signed_distance_from_triangles_everywhere(ntri, triangles, normals, work)

  endif

  ! save result to file
  call save_field_hdf5 ( time, outfile, work )

  ! finalize
  deallocate (work)
  call fft_free()
end subroutine


! this the routine that does the actual job of computing the signed distance
! for the eulerian grid. it does so everywhere without any consideration
subroutine signed_distance_from_triangles_everywhere(ntri, triangles, normals, work)
  use fsi_vars
  use stl_file_reader ! for the function distance to triangle
  implicit none
  integer, intent(in) :: ntri
  real(kind=pr),dimension(1:3,3*ntri),intent(in):: triangles
  real(kind=pr),dimension(1:3,ntri),intent(in):: normals
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix,iy,iz,i,ivertex
  real(kind=pr) :: x,y,z,dists,tmp

  ! initialize
  work = 10.0d10

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x = dx*dble(ix)
        y = dy*dble(iy)
        z = dz*dble(iz)
        dists = 9.0d9

        ! for each eulerian grid point, check all triangles. this is horribly
        ! inefficient, but any more clever method requires a lot of thinking
        ! to be generally applicable
        do i = 1,ntri
          ivertex = 3*i - 2
          ! the distance to the current triangle:
          tmp = pointTriangleDistance(triangles(:,ivertex),triangles(:,ivertex+1),&
          triangles(:,ivertex+2),(/x,y,z/),normals(:,i) )
          ! if closer (in abs value!) then use this now
          if ( dabs(tmp) < dabs(dists)) then
            dists = tmp
          endif
        enddo
        ! store value in array
        work(ix,iy,iz) = dists
      enddo
    enddo
  enddo

end subroutine signed_distance_from_triangles_everywhere




! this the routine that does the actual job of computing the signed distance
! for the eulerian grid. it does so everywhere without any consideration
subroutine signed_distance_from_triangles_surface(ntri, triangles, normals, work)
  use fsi_vars
  use helpers
  use stl_file_reader ! for the function distance to triangle
  implicit none
  integer, intent(in) :: ntri
  real(kind=pr),dimension(1:3,3*ntri),intent(in):: triangles
  real(kind=pr),dimension(1:3,ntri),intent(in):: normals
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix,iy,iz,i,ivertex,xmin,xmax,ymin,ymax,zmin,zmax,safety,mpicode
  real(kind=pr) :: x,y,z,dists,tmp
  character(len=80) :: flag

  call get_command_argument(16,flag)
  ! initialize
  work = 0.d0
  read(flag,*) safety

  if (root) write(*,*) "safety distance is ",safety

  ivertex = 1

  ! stage I: mark "interesting" points. we mark the volume (cube) spanned
  ! by each triangle with a negative number
  do i = 1,ntri
    xmin = floor( minval( triangles(1,ivertex:ivertex+2) ) /dx) -safety
    ymin = floor( minval( triangles(2,ivertex:ivertex+2) ) /dy) -safety
    zmin = floor( minval( triangles(3,ivertex:ivertex+2) ) /dz) -safety

    xmax = ceiling( maxval( triangles(1,ivertex:ivertex+2) ) /dx)+safety
    ymax = ceiling( maxval( triangles(2,ivertex:ivertex+2) ) /dy)+safety
    zmax = ceiling( maxval( triangles(3,ivertex:ivertex+2) ) /dz)+safety

    xmin = max(xmin,ra(1))
    ymin = max(ymin,ra(2))
    zmin = max(zmin,ra(3))

    xmax = min(xmax,rb(1))
    ymax = min(ymax,rb(2))
    zmax = min(zmax,rb(3))
    ivertex = ivertex + 3
    work ( xmin:xmax,ymin:ymax,zmin:zmax ) = 1.d0
  end do

  ! inform about load balancing. it is a big pity that, since only a subvolume
  ! will be concerned, load balancing can be very poor. the indication is that
  ! many cpu do not have any points to treat, while few of them have the entire
  ! load to carry. overloading (more proc than CPU) might help
  if (root) write(*,*) "note load balancing problems might occur!"
  call MPI_barrier (MPI_COMM_world, mpicode)
  ! show how many points this rank will have to treat
  write(*,'("rank=",i4," npoints=", i8," (",i3,"%)")') mpirank, nint(sum(work)),&
   100*nint(sum(work))/( (rb(1)-ra(1)+1)*(rb(2)-ra(2)+1)*(rb(3)-ra(3)+1)  )
  call MPI_barrier (MPI_COMM_world, mpicode)


  ! inform how many points are globally concerned with computing the signed distance
  tmp = mpisum( sum(work) )
  if (root) then
    write(*,'("dist will be generated only near interface, concerning ",i3"% of &
    &total points")') nint(100.d0*tmp/dble(nx*ny*nz))
  endif

  ! compute signed distance for marked points
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x = dx*dble(ix)
        y = dy*dble(iy)
        z = dz*dble(iz)
        dists = 9.0d9

        ! is the point previously marked as interesting?
        if (work(ix,iy,iz) > 0.d0) then
          ! for each eulerian grid point, check all triangles. this is horribly
          ! inefficient, but any more clever method requires a lot of thinking
          ! to be generally applicable
          do i = 1,ntri
            ivertex = 3*i - 2
            ! the distance to the current triangle:
            tmp = pointTriangleDistance(triangles(:,ivertex),triangles(:,ivertex+1),&
            triangles(:,ivertex+2),(/x,y,z/),normals(:,i) )
            ! if closer (in abs value!) then use this now
            if ( dabs(tmp) < dabs(dists)) then
              dists = tmp
            endif
          enddo
          ! store value in array
          work(ix,iy,iz) = dists
        endif
      enddo
    enddo
  enddo

end subroutine signed_distance_from_triangles_surface
