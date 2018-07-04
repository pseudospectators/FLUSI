subroutine create_mask_from_triangular_mesh(wings)

    implicit none
    type(wing), intent(inout) :: wings
    integer :: ix, iy, iz, itri, i
    integer :: ixmin, ixmax, iymin, iymax, izmin, izmax
    real(kind=pr) :: x,y,z, d

    wings(i)%mask = 9.9e9

    ! outer loop over triangles. in every triangle we loop over the union of its
    ! bounding box with the local CPUS part of the mask array
    do i = 1, nWings

      do itri = 1, maxval(wings(i)%tri_elements(:,1))

        ! determine bounding box for one triangle
        do j = 3,5
            if ( wings(i)%x(wings(i)%tri_elements(itri,j)) == &
            minval(wings(i)%x(wings(i)%tri_elements(itri,2:4)))) then
            ixmin = wings(i)%tri_elements(itri,j)
            end if

            if ( wings(i)%x(wings(i)%tri_elements(itri,j)) == &
            maxval(wings(i)%x(wings(i)%tri_elements(itri,2:4)))) then
            ixmax = wings(i)%tri_elements(itri,j)
            end if

            if ( wings(i)%y(wings(i)%tri_elements(itri,j)) == &
            minval(wings(i)%y(wings(i)%tri_elements(itri,2:4)))) then
            iymin = wings(i)%tri_elements(itri,j)
            end if

            if ( wings(i)%y(wings(i)%tri_elements(itri,j)) == &
            maxval(wings(i)%y(wings(i)%tri_elements(itri,2:4)))) then
            iymax = wings(i)%tri_elements(itri,j)
            end if

            if ( wings(i)%z(wings(i)%tri_elements(itri,j)) == &
            minval(wings(i)%z(wings(i)%tri_elements(itri,2:4)))) then
            izmin = wings(i)%tri_elements(itri,j)
            end if

            if ( wings(i)%z(wings(i)%tri_elements(itri,j)) == &
            maxval(wings(i)%z(wings(i)%tri_elements(itri,2:4)))) then
            izmax = wings(i)%tri_elements(itri,j)
            end if
        end do

    xmin = floor((wings(i)%x(ixmin) - tw - smoothing_wing)/dx) - safety
    xmax = ceiling((wings(i)%(ixmax) + tw + smoothing_wing)/dx) + safety

    ymin = floor((wings(i)%y(iymin) - tw - smoothing_wing)/dy) - safety
    ymax = ceiling((wings(i)%y(iymax) + tw + smoothing_wing)/dy) + safety

    zmin = floor((wings(i)%z(izmin) - tw - smoothing_wing)/dz) - safety
    zmax = ceiling((wings(i)%z(izmax) + tw + smoothing_wing)/dz) + safety

    do iz = zmin,zmax
      do iy = ymin,ymax
        do ix = xmin,xmax
          x = dx*dble(ix)-xl/2.0
          y = dy*dble(iy)-xl/2.0
          z = dz*dble(iz)-xl/2.0

          d = pointTriangleDistance(&
              (/wings(i)%x(wings(i)%tri_elements(itri,2)),   &
                wings(i)%y(wings(i)%tri_elements(itri,2)),   &
                wings(i)%x(wings(i)%tri_elements(itri,2))/), &
              (/wings(i)%x(wings(i)%tri_elements(itri,3)),   &
                wings(i)%y(wings(i)%tri_elements(itri,3)),   &
                wings(i)%x(wings(i)%tri_elements(itri,3))/), &
              (/wings(i)%x(wings(i)%tri_elements(itri,4)),   &
                wings(i)%y(wings(i)%tri_elements(itri,4)),   &
                wings(i)%x(wings(i)%tri_elements(itri,4))/), &
              (/x,y,z/), (/0.0_pr,0.0_pr,-1.0_pr/) )

          wings(i)%mask(ix,iy,iz) = min(wings(i)%mask(ix,iy,iz), d)
        end do
      end do
    end do
  end do

  end do
