subroutine create_mask_from_triangular_mesh(wings)

    implicit none
    type(wing), intent(inout) :: wings
    integer :: num_lines, num_cols, n_header=0, npoints, ntri, ix, iy, iz, itri, i
    real(kind=pr) :: x,y,z, d
    real(kind=pr), allocatable :: points_coordinates(:,:)
    real(kind=pr), allocatable :: triangle_indices(:,:)
    real(kind=pr), allocatable :: work(:,:,:),us(:,:,:)

    work = 9.9e9

    ! outer loop over triangles. in every triangle we loop over the union of its
    ! bounding box with the local CPUS part of the mask array
    do i = 1, nWings

    do itri = 1, ntri

      ! determine bounding box for one triangle
       do i = 2,4
         do j = 3,5
           if ( points_coordinates(triangle_indices(itri,j),i) == &
            minval(points_coordinates(triangle_indices(itri,3:5),i))) then

                 imin(i-1) = triangle_indices(itri,j)
           endif

           if ( points_coordinates(triangle_indices(itri,j),i) == &
            maxval(points_coordinates(triangle_indices(itri,3:5),i))) then

                 imax(i-1) = triangle_indices(itri,j)
           end if
         end do
       end do


       xmin = floor((points_coordinates(imin(1),2) - tw - smoothing_wing)/dx) - safety
       xmax = ceiling((points_coordinates(imax(1),2) + tw + smoothing_wing)/dx) + safety

       ymin = floor((points_coordinates(imin(2),3) - tw - smoothing_wing)/dy) - safety
       ymax = ceiling((points_coordinates(imax(2),3) + tw + smoothing_wing)/dy) + safety

       zmin = floor((points_coordinates(imin(2),4) - tw - smoothing_wing)/dz) - safety
       zmax = ceiling((points_coordinates(imax(2),4) + tw + smoothing_wing)/dz) + safety

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

  end do
