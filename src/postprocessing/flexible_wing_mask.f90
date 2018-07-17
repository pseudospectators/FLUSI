subroutine flexible_wing_mask(help)
    use vars
    use stl_file_reader
    use module_helpers
    use module_ini_files_parser_mpi
    use p3dfft_wrapper
    use penalization, only : mask_color
    use module_insects
    implicit none
    logical, intent(in) :: help
    character(len=strlen) :: file
    integer :: num_lines, num_cols, n_header=0, npoints, ntri, ix, iy, iz, itri
    real(kind=pr) :: x,y,z, d
    real(kind=pr), allocatable :: points_coordinates(:,:)
    real(kind=pr), allocatable :: triangle_indices(:,:)
    real(kind=pr), allocatable :: work(:,:,:)


    if (help.and.root) then
      write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      write(*,*) "./flusi -p --flexible-wing-mask "
      write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      write(*,*) "! "
      write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      write(*,*) "Parallel: yes"
      return
    endif

    ! first thing to do would be load ascii DAT files
    file = "points_coor.dat"
    call count_lines_in_ascii_file_mpi(file, npoints, n_header)
    call count_cols_in_ascii_file_mpi(file, num_cols, n_header)
    allocate( points_coordinates(1:npoints, 1:num_cols) )
    call read_array_from_ascii_file_mpi(file, points_coordinates, n_header)


    file = "mesh_triangle_elements.dat"
    call count_lines_in_ascii_file_mpi(file, ntri, n_header)
    call count_cols_in_ascii_file_mpi(file, num_cols, n_header)
    allocate( triangle_indices(1:ntri, 1:num_cols) )
    call read_array_from_ascii_file_mpi(file, triangle_indices, n_header)

    xl = 0.2_pr
    yl = 0.2_pr
    zl = 0.2_pr
    nx = 64
    ny = 64
    nz = 64

    ! set up MPI and arrays bounds (ra, rb)
    call decomposition_initialize()
    ! this will be our mask array
    allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    ! initialize the array as very large distance
    ! in the first step, we use this array to create the DISTANCE function, later the MASK
    work = 9.9e9

    ! outer loop over triangles. in every triangle we loop over the union of its
    ! bounding box with the local CPUS oart of the mask array
    do itri = 1, ntri

        ! HACK NOTE: Modify bounding box
        do iz = ra(3),rb(3)
            do iy = ra(2),rb(2)
                do ix = ra(1),rb(1)
                    x = dx*dble(ix)-xl/2.0
                    y = dy*dble(iy)-xl/2.0
                    z = dz*dble(iz)-xl/2.0

                    d = pointTriangleDistance( points_coordinates( triangle_indices(itri,3),2:4), &
                    points_coordinates( triangle_indices(itri,4),2:4), &
                    points_coordinates( triangle_indices(itri,5),2:4), &
                    (/x,y,z/), (/0.0_pr,0.0_pr,-1.0_pr/) )

                    work(ix,iy,iz) = min( work(ix,iy,iz), d)
                end do
            end do
        end do
    end do


    call save_field_hdf5 ( 0.0d0, 'mask_000.h5', work )

    ! create the mask function

    deallocate( points_coordinates, triangle_indices, work)
end subroutine
