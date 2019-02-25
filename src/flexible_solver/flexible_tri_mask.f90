! TO DO migrate this one into create_mask_fsi subroutine
subroutine Draw_flexible_wing(time, wings, mask, mask_color, us)!, unsigned_distance)

  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  type(flexible_wing),dimension(1:nWings), intent(inout) :: wings

  ! initialize everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  ! Create mask function and us field from triangular mesh
  call create_mask_from_triangular_mesh(wings,mask,us,mask_color)

  !deallocate(unsigned_distance(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

end subroutine Draw_flexible_wing

subroutine calculate_normal_vectors_of_wing(wings)

 implicit none

 type(flexible_wing),dimension(1:nWings), intent(inout) :: wings
 real(kind=pr),allocatable :: normal(:,:)
 integer :: i, itri

do i = 1, nWings

allocate(normal(1:wings(i)%ntri,1:3))

 do itri = 1, wings(i)%ntri



   ! Calculate the normal vector of one triangle
     normal(itri,1:3) = cross((/wings(i)%x(wings(i)%tri_elements(itri,2)) - &
                              wings(i)%x(wings(i)%tri_elements(itri,3)),  &
                              wings(i)%y(wings(i)%tri_elements(itri,2)) - &
                              wings(i)%y(wings(i)%tri_elements(itri,3)),  &
                              wings(i)%z(wings(i)%tri_elements(itri,2)) - &
                              wings(i)%z(wings(i)%tri_elements(itri,3))/),&
                            (/wings(i)%x(wings(i)%tri_elements(itri,3)) - &
                              wings(i)%x(wings(i)%tri_elements(itri,4)),  &
                              wings(i)%y(wings(i)%tri_elements(itri,3)) - &
                              wings(i)%y(wings(i)%tri_elements(itri,4)),  &
                              wings(i)%z(wings(i)%tri_elements(itri,3)) - &
                              wings(i)%z(wings(i)%tri_elements(itri,4))/))

     ! dimentionalized to get a unit vector
     wings(i)%tri_element_normals(itri,1) = normal(itri,1)/norm2(normal(itri,1:3))
     wings(i)%tri_element_normals(itri,2) = normal(itri,2)/norm2(normal(itri,1:3))
     wings(i)%tri_element_normals(itri,3) = normal(itri,3)/norm2(normal(itri,1:3))

  enddo

deallocate(normal)

enddo


end subroutine


subroutine create_mask_from_triangular_mesh(wings,mask,us,mask_color)

    implicit none
    real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    !real(kind=pr),intent(inout)::unsigned_distance(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
    integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    type(flexible_wing),dimension(1:nWings), intent(inout) :: wings
    integer :: ix, iy, iz, itri, i, j
    integer :: ixmin, ixmax, iymin, iymax, izmin, izmax
    integer :: xmin, xmax, ymin, ymax, zmin, zmax
    integer, parameter :: safety = 2
    real(kind=pr) :: x,y,z, distance
    real(kind=pr),dimension(1:3) :: velocity

    real(kind=pr),allocatable :: unsigned_distance(:,:,:)


    allocate(unsigned_distance(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    unsigned_distance = 1.d10 !assign the distance to be really far away


    do i = 1, nWings


      ! outer loop over triangles. in every triangle we loop over the union of its
      ! bounding box with the local CPUS part of the mask array
      do itri = 1, wings(i)%ntri


        ! Determine bounding box for one triangle
        ixmin = wings(i)%tri_elements(itri,&
                minloc(wings(i)%x(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)

        ixmax = wings(i)%tri_elements(itri,&
                maxloc(wings(i)%x(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)

        iymin = wings(i)%tri_elements(itri,&
                minloc(wings(i)%y(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)

        iymax = wings(i)%tri_elements(itri,&
                maxloc(wings(i)%y(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)

        izmin = wings(i)%tri_elements(itri,&
                minloc(wings(i)%z(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)

        izmax = wings(i)%tri_elements(itri,&
                maxloc(wings(i)%z(wings(i)%tri_elements(itri,2:4)),DIM=1) + 1)


        xmin = floor((wings(i)%x(ixmin) - wings(i)%t_wing &
                    - wings(i)%wing_smoothing)/dx) - safety
        xmax = ceiling((wings(i)%x(ixmax) + wings(i)%t_wing &
                    + wings(i)%wing_smoothing)/dx) + safety

        ymin = floor((wings(i)%y(iymin) - wings(i)%t_wing &
                    - wings(i)%wing_smoothing)/dy) - safety
        ymax = ceiling((wings(i)%y(iymax) + wings(i)%t_wing &
                      + wings(i)%wing_smoothing)/dy) + safety

        zmin = floor((wings(i)%z(izmin) - wings(i)%t_wing &
                    - wings(i)%wing_smoothing)/dz) - safety
        zmax = ceiling((wings(i)%z(izmax) + wings(i)%t_wing &
                      + wings(i)%wing_smoothing)/dz) + safety

        ! inner loops over all Eulerian grid nodes inside the bounding box
        ! containing the Lagrangian triangular element to calculate the unsigned
        ! distance between the fluid grid nodes and the solid grid nodes
        do iz = max(ra(3),zmin),min(rb(3),zmax)
          do iy = max(ra(2),ymin),min(rb(2),ymax)
            do ix = max(ra(1),xmin),min(rb(1),xmax)

              x = dx*dble(ix)!-xl/2
              y = dy*dble(iy)!-yl/2
              z = dz*dble(iz)!-zl/2


              call calculate_unsigned_distance_and_us(distance, velocity, &
                  (/wings(i)%x(wings(i)%tri_elements(itri,2)),   &
                    wings(i)%y(wings(i)%tri_elements(itri,2)),   &
                    wings(i)%z(wings(i)%tri_elements(itri,2))/), &
                  (/wings(i)%x(wings(i)%tri_elements(itri,3)),   &
                    wings(i)%y(wings(i)%tri_elements(itri,3)),   &
                    wings(i)%z(wings(i)%tri_elements(itri,3))/), &
                  (/wings(i)%x(wings(i)%tri_elements(itri,4)),   &
                    wings(i)%y(wings(i)%tri_elements(itri,4)),   &
                    wings(i)%z(wings(i)%tri_elements(itri,4))/), &
                  (/wings(i)%vx(wings(i)%tri_elements(itri,2)),   &
                    wings(i)%vy(wings(i)%tri_elements(itri,2)),   &
                    wings(i)%vz(wings(i)%tri_elements(itri,2))/), &
                  (/wings(i)%vx(wings(i)%tri_elements(itri,3)),   &
                    wings(i)%vy(wings(i)%tri_elements(itri,3)),   &
                    wings(i)%vz(wings(i)%tri_elements(itri,3))/), &
                  (/wings(i)%vx(wings(i)%tri_elements(itri,4)),   &
                    wings(i)%vy(wings(i)%tri_elements(itri,4)),   &
                    wings(i)%vz(wings(i)%tri_elements(itri,4))/), &
                  (/x,y,z/), wings(i)%tri_element_normals(itri,1:3))


              if (distance < unsigned_distance(ix,iy,iz)) then
                  unsigned_distance(ix,iy,iz) = distance
                  us(ix,iy,iz,1:3) = velocity
              endif

            enddo
          enddo
        enddo

        ! Then we calculate mask function from unsigned distance
        do iz = max(ra(3),zmin),min(rb(3),zmax)
          do iy = max(ra(2),ymin),min(rb(2),ymax)
            do ix = max(ra(1),xmin),min(rb(1),xmax)

                call smoothstep(mask(ix,iy,iz),unsigned_distance(ix,iy,iz),&
                                wings(i)%t_wing,wings(i)%wing_smoothing)

                !-- assign mask color
                if (mask(ix,iy,iz) > 0.d0) mask_color(ix,iy,iz)=1

            enddo
          enddo
        enddo

      enddo

  enddo

  deallocate(unsigned_distance)

end subroutine create_mask_from_triangular_mesh

subroutine calculate_unsigned_distance_and_us(distance,us,tri1,tri2,tri3,utri1,utri2,utri3,point,normal)
  ! calculate distance between a point and a triangle in 3D
  ! SYNTAX
  !   dist = pointTriangleDistance(TRI,P)
  !   [dist,PP0] = pointTriangleDistance(TRI,P)
  !
  ! DESCRIPTION
  !   Calculate the distance of a given point P from a triangle TRI.
  !   Point P is a row vector of the form 1x3. The triangle is a matrix
  !   formed by three rows of points TRI = [P1P2P3] each of size 1x3.
  !   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
  !   to the triangle TRI.
  !   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
  !   closest point PP0 to P on the triangle TRI.
  !
  ! Author: Gwendolyn Fischer
  ! Release: 1.0
  ! Release date: 09/02/02
  ! Release: 1.1 Fixed Bug because of normalization
  ! Release: 1.2 Fixed Bug because of typo in region 5 20101013
  ! Release: 1.3 Fixed Bug because of typo in region 2 20101014

  ! Possible extention could be a version tailored not to return the distance
  ! and additionally the closest point, but instead return only the closest
  ! point. Could lead to a small speed gain.

  ! Example:
  ! !! The Problem
  ! P0 = [0.5 -0.3 0.5]
  !
  ! P1 = [0 -1 0]
  ! P2 = [1  0 0]
  ! P3 = [0  0 0]
  !
  ! vertices = [P1 P2 P3]
  ! faces = [1 2 3]
  !
  ! !! The Engine
  ! [dist,PP0] = pointTriangleDistance([P1P2P3],P0)
  !
  ! !! Visualization
  ! [x,y,z] = sphere(20)
  ! x = dist*x+P0(1)
  ! y = dist*y+P0(2)
  ! z = dist*z+P0(3)
  !
  ! figure
  ! hold all
  ! patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
  ! plot3(P0(1),P0(2),P0(3),'b*')
  ! plot3(PP0(1),PP0(2),PP0(3),'*g')
  ! surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
  ! view(3)

  ! The algorithm is based on
  ! "David Eberly, 'Distance Between Point and Triangle in 3D',
  ! Geometric Tools, LLC, (1999)"
  ! http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
  !
  !        ^t
  !  \     |
  !   \reg2|
  !    \   |
  !     \  |
  !      \ |
  !       \|
  !        *P2
  !        |\
  !        | \
  !  reg3  |  \ reg1
  !        |   \
  !        |reg0\
  !        |     \
  !        |      \ P1
  ! -------*-------*------->s
  !        |P0      \
  !  reg4  | reg5    \ reg6
  !
  !Modification by Hung 31/05/2018:
  !     Changing from function to subroutine to get out both distance and us
  !     The calculation of us is added for fluid structure interaction
  !
  implicit none
  real(kind=pr), intent(out) :: distance
  real(kind=pr), dimension(1:3), intent(out) :: us
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,tri3,utri1,utri2,utri3
  real(kind=pr), dimension(1:3), intent(in) :: point,normal
  real(kind=pr), dimension(1:3) :: BB,EE0,EE1,DD
  real(kind=pr) :: a,b,c,d,e,f,det,s,t,sqrDistance,tmp0,tmp1,numer,denom,invDet

  ! rewrite triangle in normal form
  BB = tri1
  EE0 = tri2-BB
  EE1 = tri3-BB


  DD = BB - point
  a = dot_product(EE0,EE0)
  b = dot_product(EE0,EE1)
  c = dot_product(EE1,EE1)
  d = dot_product(EE0,DD)
  e = dot_product(EE1,DD)
  f = dot_product(DD,DD)



  det = a*c - b*b ! do we have to use abs here?
  s   = b*e - c*d
  t   = b*d - a*e

  if (det < 1.0d-15) then
    distance = 9.0d9
    return
  endif


  ! write(*,'(12(es12.4,1x))') tri1,tri2,tri3,point
  ! write(*,'(12(es12.4,1x))') a,b,c,d,e,f,det,s,t

  ! Terible tree of conditionals to determine in which region of the diagram
  ! shown above the projection of the point into the triangle-plane lies.
  if ((s+t) <= det) then
    if (s < 0.d0) then
      if (t < 0.d0) then
        !region4
        if (d < 0.d0) then
          t = 0.d0
          if (-d >= a) then
            s = 1.d0
            sqrDistance = a + 2.d0*d + f

            us = utri2

          else
            s = -d/a
            sqrDistance = d*s + f

            us = interpolationPointLine(tri1,tri2,utri1,utri2,point)

          endif
        else
          s = 0.d0
          if (e >= 0.d0) then
            t = 0.d0
            sqrDistance = f

            us = utri1

          else
            if (-e >= c) then
              t = 1.d0
              sqrDistance = c + 2.d0*e + f

              us = utri3

            else
              t = -e/c
              sqrDistance = e*t + f

              us = interpolationPointLine(tri3,tri1,utri3,utri1,point)

            endif
          endif
        endif !of region 4
      else
        ! region 3
        s = 0.d0
        if (e >= 0.d0) then
          t = 0.d0
          sqrDistance = f

          us = utri1

        else
          if (-e >= c) then
            t = 1.d0
            sqrDistance = c + 2.d0*e +f

            us = utri3

          else
            t = -e/c
            sqrDistance = e*t + f

            us = interpolationPointLine(tri3,tri1,utri3,utri1,point)

          endif
        endif
      endif !of region 3
    else
      if (t < 0.d0) then
        ! region 5
        t = 0.d0
        if (d >= 0.d0) then
          s = 0.d0
          sqrDistance = f

          us = utri1

        else
          if (-d >= a) then
            s = 1.d0
            sqrDistance = a + 2.d0*d + f! GF 20101013 fixed typo d*s ->2*d

            us = utri2

          else
            s = -d/a
            sqrDistance = d*s + f

            us = interpolationPointLine(tri1,tri2,utri1,utri2,point)

          endif
        endif
      else
        ! region 0
        invDet = 1.d0/det
        s = s*invDet
        t = t*invDet
        sqrDistance = s*(a*s + b*t + 2.d0*d) &
                    + t*(b*s + c*t + 2.d0*e) + f

        us = interpolationPointTriangle(tri1,tri2,tri3,utri1,utri2,utri3,point,normal)

      endif
    endif
  else
    if (s < 0.d0) then
      ! region 2
      tmp0 = b + d
      tmp1 = c + e
      if (tmp1 > tmp0) then ! minimum on edge s+t=1
        numer = tmp1 - tmp0
        denom = a - 2.d0*b + c
        if (numer >= denom) then
          s = 1.d0
          t = 0.d0
          sqrDistance = a + 2.d0*d + f ! GF 20101014 fixed typo 2*b -> 2*d

          us = utri2

        else
          s = numer/denom
          t = 1.d0-s
          sqrDistance = s*(a*s + b*t + 2.d0*d) &
                      + t*(b*s + c*t + 2.d0*e) + f

          us = interpolationPointLine(tri2,tri3,utri2,utri3,point)

        endif
      else          ! minimum on edge s=0
        s = 0.d0
        if (tmp1 <= 0.d0) then
          t = 1.d0
          sqrDistance = c + 2.d0*e + f

          us = utri3

        else
          if (e >= 0.d0) then
            t = 0.d0
            sqrDistance = f

            us = utri1
          else
            t = -e/c
            sqrDistance = e*t + f

            us = interpolationPointLine(tri3,tri1,utri3,utri1,point)

          endif
        endif
      endif !of region 2
    else
      if (t < 0.d0) then
        !region6
        tmp0 = b + e
        tmp1 = a + d
        if (tmp1 > tmp0) then ! minimum on  edge  s + t = 1  with  t > 0
          numer = tmp1 - tmp0
          denom = a-2.d0*b+c
          if (numer >= denom) then
            t = 1.d0
            s = 0.d0
            sqrDistance = c + 2.d0*e + f

            us = utri3

          else
            t = numer/denom
            s = 1.d0 - t
            sqrDistance = s*(a*s + b*t + 2.d0*d) &
                        + t*(b*s + c*t + 2.d0*e) + f

            us = interpolationPointLine(tri2,tri3,utri2,utri3,point)

          endif
        else !minimum on  edge  t = 0  with  s <= 1
          t = 0.d0
          if (tmp1 <= 0) then
              s = 1.d0
              sqrDistance = a + 2.d0*d + f

              us = utri2

          else
            if (d >= 0.d0) then
                s = 0.d0
                sqrDistance = f

                us = utri1

            else
                s = -d/a
                sqrDistance = d*s + f

                us = interpolationPointLine(tri1,tri2,utri1,utri2,point)

            endif
          endif
        endif
        !end region 6
      else
        ! region 1
        numer = c + e - b - d
        if (numer <= 0.d0) then
          s = 0.d0
          t = 1.d0
          sqrDistance = c + 2.d0*e + f

          us = utri3

        else
          denom = a - 2.d0*b + c
          if (numer >= denom) then
            s = 1.d0
            t = 0.d0
            sqrDistance = a + 2.d0*d + f

            us = utri2

          else
            s = numer/denom
            t = 1-s
            sqrDistance = s*(a*s + b*t + 2.d0*d) &
                        + t*(b*s + c*t + 2.d0*e) + f

            us = interpolationPointLine(tri2,tri3,utri2,utri3,point)

          endif
        endif !of region 1
      endif
    endif
  endif

  ! account for numerical round-off error
  if (sqrDistance < 0.d0) then
    sqrDistance = 0.d0
  endif



  ! closest point on triangle
  DD = BB + s*EE0 + t*EE1;
  ! vector from target point to closest point on surface
  DD = point-DD
  t = dot_product(DD,normal)

  !write(*,*) "WARNING CHANGED TO UNSIGNED DISTANCE!!"
  ! if (t >= 0.d0) then
    distance = dsqrt(sqrDistance)
  ! else
  !   pointTriangleDistance = -dsqrt(sqrDistance)
  ! endif

end subroutine calculate_unsigned_distance_and_us

subroutine project_point_onto_triangle(point_projected,tri1,point,normal)

  !Give the projection of a point onto a plan constructed from a normal vector and
  !a point belonging to the plan given by the equation
  !                 ax + by + cz + d = 0

  implicit none
  real(kind=pr), dimension(1:3), intent(in) :: tri1,point,normal
  real(kind=pr), dimension(1:3), intent(out) :: point_projected
  real(kind=pr) :: a,b,c,d

  a = normal(1)
  b = normal(2)
  c = normal(3)
  d = - (a*tri1(1) + b*tri1(2) + c*tri1(3))

  point_projected(1) = ((b**2 + c**2)*point(1) - a*b*point(2) - a*c*point(3) - d*a)/ &
                       (a**2 + b**2 + c**2)

  point_projected(2) = (-a*b*point(1) + (a**2 + c**2)*point(2) - b*c*point(3) - d*b)/ &
                       (a**2 + b**2 + c**2)

  point_projected(3) = (-a*c*point(1) - b*c*point(2) + (a**2 + b**2)*point(3) - d*c)/ &
                       (a**2 + b**2 + c**2)

end subroutine

function interpolationPointTriangle(tri1,tri2,tri3,utri1,utri2,utri3,point,normal)

  !Interpolate value of three vertices of a triangle onto a point inside that
  !triangle using barycentric interpolation

  implicit none
  real(kind=pr), dimension(1:3) :: interpolationPointTriangle
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,tri3,utri1,utri2,utri3,normal
  real(kind=pr), dimension(1:3), intent(in) :: point
  real(kind=pr), dimension(1:3) :: point_projected, tmp
  real(kind=pr) :: A1, A2, A3, A

  !First we need to project the point onto the plan of the triangle
  call project_point_onto_triangle(point_projected,tri1,point,normal)

  !
  tmp = cross(tri2-point_projected,tri3-point_projected)
  A1 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1-point_projected,tri3-point_projected)
  A2 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1-point_projected,tri2-point_projected)
  A3 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1-tri3,tri2-tri3)
  A = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  !Check if the projection of the point is inside the triangle

  if ((A1 + A2 + A3) > (A + 1.0d-15)) then
    write(*,*) "WARNING: The projection of the point may be outside of the triangle"
    write(*,*) A1, A2, A3
    write(*,*) A
    write(*,*) tri1, tri2, tri3
    write(*,*) point
    write(*,*) point_projected
  endif

  interpolationPointTriangle(1:3) = (A1*utri1(1:3) + A2*utri2(1:3) + A3*utri3(1:3))/ &
                                    (A1 + A2 + A3)

end function

subroutine project_point_onto_line(point_projected,tri1,tri2,point)

  !Give the projection of a point onto a line constructed from two vertices tri1
  !and tri2 of a triangle

  implicit none
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,point
  real(kind=pr), dimension(1:3), intent(out) :: point_projected

  point_projected(1:3) = tri1 + (dot_product(point - tri1,tri2 - tri1)/ &
                         dot_product(tri2 - tri1,tri2 - tri1))*(tri2 - tri1)

end subroutine

function interpolationPointLine(tri1,tri2,utri1,utri2,point)

  !Interpolate value of three vertices of a triangle onto a point inside that
  !triangle using barycentric interpolation

  implicit none
  real(kind=pr), dimension(1:3) :: interpolationPointLine
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,utri1,utri2
  real(kind=pr), dimension(1:3), intent(in) :: point
  real(kind=pr), dimension(1:3) :: point_projected, tmp
  real(kind=pr) :: D1, D2, D

  !First we need to project the point onto the edge of the triangle
  call project_point_onto_line(point_projected,tri1,tri2,point)

  !
  tmp = point_projected - tri2
  D1 = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = point_projected - tri1
  D2 = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = tri1 - tri2
  D = sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  !Check if the projection of the point is inside the triangle
  if (root) then
    if ((D1 + D2) > D + 1.0d-15) then
    write(*,*) "WARNING: The projection of the point may be out of the line"
    write(*,*) D1, D2
    write(*,*) D
    write(*,*) tri1, tri2
    write(*,*) point
    write(*,*) point_projected
    endif
 endif

  interpolationPointLine(1:3) = (D1*utri1(1:3) + D2*utri2(1:3))/ &
                                (D1 + D2)

end function interpolationPointLine
