subroutine stl2dist(help)
    use vars
    use module_helpers
    use stl_file_reader
    use p3dfft_wrapper
    implicit none
    logical, intent(in) :: help
    character(len=strlen) :: stlfile, outfile, mode, dummy, flag
    real(kind=pr),dimension(:,:),allocatable :: triangles,normals
    real(kind=pr),dimension(:,:,:),allocatable :: work
    real(kind=pr) :: time, scale
    integer :: ntri
    integer :: mpierror


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

    time = 0.d0
    ! initialize code and domain decomposition, but do not use FFTs
    call decomposition_initialize(.true.)

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

    call MPI_barrier(MPI_COMM_WORLD, mpierror)

    ! save result to file
    call save_field_hdf5 ( time, outfile, work )

    ! finalize
    deallocate (work)
    call fft_free()
end subroutine


! this the routine that does the actual job of computing the signed distance
! for the eulerian grid. it does so everywhere without any consideration
subroutine signed_distance_from_triangles_everywhere(ntri, triangles, normals, work)
    use vars
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
    use vars
    use module_helpers
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
    read(flag,*) safety
    write(*,*) "safety distance is ", safety
    write(*,*) "number of triangles is ", ntri

    ! initialize
    work = 7.0d9

    do i = 1, ntri
        ivertex = 3*i - 2

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

        do iz = zmin, zmax
            z = dz*dble(iz)
            do iy = ymin, ymax
                y = dy*dble(iy)
                do ix = xmin, xmax
                    x = dx*dble(ix)

                    ! the distance to the current triangle:
                    tmp = pointTriangleDistance( triangles(:,ivertex),triangles(:,ivertex+1),&
                    triangles(:,ivertex+2), (/x,y,z/), normals(:,i) )

                    ! write(*,*) tmp

                    if (abs(tmp)<dble(safety)*maxval((/dx,dy,dz/))) then
                        ! if closer (in abs value!) then use this now
                        if ( dabs(tmp) < dabs(work(ix,iy,iz))) then
                            work(ix,iy,iz) = tmp
                        endif
                    endif
                enddo
            enddo
        enddo

    end do
end subroutine signed_distance_from_triangles_surface
