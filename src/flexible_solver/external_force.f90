!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct external force vector consists of gravity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine external_forces_construction(time,dt0,dt1, it,Wings)
! This is actually just for 1 wing
use mpi
implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing),dimension(1:nWings), intent(inout)  :: Wings
integer :: i, j, np

do i=1,nWings

  np = wings(i)%np
  ! Initialize
  wings(i)%Fext = 0.d0

  ! Gravitational forces
  call gravitational_forces_on_wing (wings(i))

  ! Forces from the fluid pressure field
  if (activate_press_force=="yes") call pressure_forces_on_wing (wings(i))

  ! Fictitious forces appear when the reference frames of wings are non-inertial frame
  !call fictitious_forces_of_moving_reference_frame (time,dt0,dt1,it,wings(i))



enddo

end subroutine

subroutine gravitational_forces_on_wing (wing)
!Calculate the forces acting on ONE wing by the gravitational field

implicit none

type(flexible_wing), intent (inout) :: wing
integer :: j, np

  np = wing%np

  do j=1,np
    wing%Fext(1:np)        = grav(1)*wing%m(j) !forces on the x-direction
    wing%Fext(np+1:2*np)   = grav(2)*wing%m(j) !forces on the y-direction
    wing%Fext(2*np+1:3*np) = grav(3)*wing%m(j) !forces on the z-direction
  enddo

end subroutine

subroutine pressure_forces_on_wing (wing)
!Calculate the forces acting on ONE wing by the fluid pressure field

implicit none

type(flexible_wing), intent (inout) :: wing
integer :: j, np

  call transform_pressure_into_point_forces_per_node(wing)

end subroutine

subroutine get_pressure_on_wing_surfaces(wings,pressure_field)
use mpi
implicit none

type(flexible_wing),dimension(1:nWings), intent(inout) :: Wings
real(kind=pr), intent(inout) :: pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
integer :: i, itri, ntri, np
real(kind=pr), allocatable :: upside(:,:,:), downside(:,:,:)
logical :: periodic

periodic =.false.

do i = 1, nWings

  np = wings(i)%np
  allocate(upside(1:wings(i)%ntri,1:3,1:3),downside(1:wings(i)%ntri,1:3,1:3))


    !Getting surface position from the centerline of the wing
    call calculate_wing_surfaces(upside(1:wings(i)%ntri,1:3,1:3), &
                                 wings(i)%x(1:np), &
                                 wings(i)%y(1:np), &
                                 wings(i)%z(1:np), &
                                 wings(i)%t_wing, &
                                 wings(i)%tri_elements(1:wings(i)%ntri,1:4), &
                                wings(i)%tri_element_normals(1:wings(i)%ntri,1:4),"upside")

    call calculate_wing_surfaces(downside(1:wings(i)%ntri,1:3,1:3), &
                                 wings(i)%x(1:np), &
                                 wings(i)%y(1:np), &
                                 wings(i)%z(1:np), &
                                 wings(i)%t_wing, &
                                 wings(i)%tri_elements(1:wings(i)%ntri,1:4), &
                                 wings(i)%tri_element_normals(1:wings(i)%ntri,1:4),"downside")
      !write(*,*) "downside surface z"
      !write(*,*) downside(1:wings(i)%ntri,1:3,3)

    !Interpolate pressure at wing surface from the pressure field of fluid
    call interpolate_pressure_on_wing_surfaces(wings(i)%press_upside(1:wings(i)%ntri,1:3),&
                                              (/ga(1)*dx,ga(2)*dy,ga(3)*dz/),(/dx,dy,dz/),&
                                               pressure_field,upside,wings(i)%ntri,periodic)

    call interpolate_pressure_on_wing_surfaces(wings(i)%press_downside(1:wings(i)%ntri,1:3),&
                                              (/ga(1)*dx,ga(2)*dy,ga(3)*dz/),(/dx,dy,dz/),&
                                               pressure_field,downside,wings(i)%ntri,periodic)


      !write(*,'("Im CPU ",i5," gives pressure",f7.3," and ",f7.3," &
      !          while upside values are ",f7.3," at ",i5," ",i5," ",i5," and ",f7.3," at ",i5," ",i5," ",i5,".")') mpirank, &
      !          maxval(pressure_field), minval(pressure_field), &
      !          minval(wings(i)%press_upside(1:wings(i)%ntri,1:3)), minloc(wings(i)%press_upside(1:wings(i)%ntri,1)),&
      !          minloc(wings(i)%press_upside(1:wings(i)%ntri,2)), minloc(wings(i)%press_upside(1:wings(i)%ntri,3)),&
      !          maxval(wings(i)%press_upside(1:wings(i)%ntri,1:3)), maxloc(wings(i)%press_upside(1:wings(i)%ntri,1)),&
      !          maxloc(wings(i)%press_upside(1:wings(i)%ntri,2)), maxloc(wings(i)%press_upside(1:wings(i)%ntri,3))
      !write(*,'("Im CPU ",i5," gives pressure",f7.3," and ",f7.3," &
      !          while downside values are ",f7.3," at ",i5," ",i5," ",i5," and ",f7.3," at ",i5," ",i5," ",i5,".")') mpirank, &
      !          maxval(pressure_field), minval(pressure_field), &
      !          minval(wings(i)%press_downside(1:wings(i)%ntri,1:3)), minloc(wings(i)%press_downside(1:wings(i)%ntri,1)),&
      !          minloc(wings(i)%press_downside(1:wings(i)%ntri,2)), minloc(wings(i)%press_downside(1:wings(i)%ntri,3)),&
      !          maxval(wings(i)%press_downside(1:wings(i)%ntri,1:3)), maxloc(wings(i)%press_downside(1:wings(i)%ntri,1)),&
      !          maxloc(wings(i)%press_downside(1:wings(i)%ntri,2)), maxloc(wings(i)%press_downside(1:wings(i)%ntri,3))



  deallocate(upside,downside)

enddo

end subroutine

subroutine calculate_wing_surfaces(surface,x,y,z,thickness,tri_elements,tri_element_normals,position)
use mpi
  implicit none

  real(kind=pr), intent (inout) :: surface(:,:,:)
  real(kind=pr), intent(in) :: x(:),y(:),z(:)
  real(kind=pr), intent(in) :: thickness
  real(kind=pr), intent(in) :: tri_element_normals(:,:)
  integer, intent(in) :: tri_elements(:,:)
  character(len=*), intent(in) :: position
  integer :: itri, ntri

  ntri = size(tri_elements,DIM=1)

  select case(position)
    case ("upside")
      do itri=1,ntri
           surface(itri,1,1:3) = (/x(tri_elements(itri,2)), &
                                  y(tri_elements(itri,2)), &
                                  z(tri_elements(itri,2))/) &
                                + thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                                !+ thickness*(/0.d0,0.d0,1.d0/)
           surface(itri,2,1:3) = (/x(tri_elements(itri,3)), &
                                 y(tri_elements(itri,3)), &
                                 z(tri_elements(itri,3))/) &
                               + thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                               !+ thickness*(/0.d0,0.d0,1.d0/)
           surface(itri,3,1:3) = (/x(tri_elements(itri,4)), &
                                 y(tri_elements(itri,4)), &
                                 z(tri_elements(itri,4))/) &
                               + thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                              !+ thickness*(/0.d0,0.d0,1.d0/)
      enddo
    case ("downside")
      do itri=1,ntri
          surface(itri,1,1:3) = (/x(tri_elements(itri,2)), &
                                 y(tri_elements(itri,2)), &
                                 z(tri_elements(itri,2))/) &
                               - thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                               !+ thickness*(/0.d0,0.d0,1.d0/)
          surface(itri,2,1:3) = (/x(tri_elements(itri,3)), &
                                y(tri_elements(itri,3)), &
                                z(tri_elements(itri,3))/) &
                              - thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                              !+ thickness*(/0.d0,0.d0,1.d0/)
          surface(itri,3,1:3) = (/x(tri_elements(itri,4)), &
                                y(tri_elements(itri,4)), &
                                z(tri_elements(itri,4))/) &
                             - thickness*tri_element_normals(itri,4)*tri_element_normals(itri,1:3)
                            !  + thickness*(/0.d0,0.d0,1.d0/)
      enddo
  end select

end subroutine

subroutine interpolate_pressure_on_wing_surfaces(press_surface,origin,grid_space,pressure_field,surface,ntri,periodic)

  use interpolation
  use mpi

  implicit none


  real(kind=pr), intent(inout) :: press_surface(:,:)
  real(kind=pr), intent(inout) :: pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr), intent(in) :: origin(1:3), grid_space(1:3)
  real(kind=pr), intent(in) :: surface(:,:,:)
  real(kind=pr), allocatable :: press_surface_local(:,:),p_surface(:,:)
  real(kind=pr) :: min_press, max_press
  integer, intent(in) :: ntri
  logical, intent(in) :: periodic
  integer :: itri, mpicode

  allocate(press_surface_local(1:ntri,1:3),p_surface(1:ntri,1:3))

  press_surface_local = 17.d0
  p_surface = 17.d0

  min_press = minval(pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))-1.d-6
  max_press = maxval(pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))+1.d-6

      do itri=1,ntri


          !press_surface_local(itri,1) =  pressure_trilinear_interp(origin, grid_space, pressure_field,&
          !                                                         surface(itri,1,1:3), periodic)
          !press_surface_local(itri,2) =  pressure_trilinear_interp(origin, grid_space, pressure_field,&
          !                                                         surface(itri,2,1:3), periodic)
          !press_surface_local(itri,3) =  pressure_trilinear_interp(origin, grid_space, pressure_field,&
          !                                                         surface(itri,3,1:3), periodic)

          call delta_interpolation( surface(itri,1,1:3), &
                                    pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)), press_surface_local(itri,1))
          call delta_interpolation( surface(itri,2,1:3), &
                                    pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)), press_surface_local(itri,2))
          call delta_interpolation( surface(itri,3,1:3), &
                                    pressure_field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)), press_surface_local(itri,3))


      enddo

      call MPI_ALLREDUCE (press_surface_local,p_surface,size(press_surface_local),&
                       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)



      press_surface = p_surface!press_surface_local!p_surface


end subroutine

subroutine transform_pressure_into_point_forces_per_triangle(wing)

  implicit none

  type(flexible_wing), intent (inout) :: wing
  real(kind=pr), allocatable :: force_upside(:), force_downside(:)
  real(kind=pr), allocatable :: centroid_upside_projected(:,:), centroid_downside_projected(:,:)
  real(kind=pr), allocatable :: upside(:,:,:), downside(:,:,:)
  real(kind=pr), allocatable :: upside_distributed_forces(:,:,:), downside_distributed_forces(:,:,:)
  integer :: itri,np

  allocate(upside(1:wing%ntri,1:3,1:3),downside(1:wing%ntri,1:3,1:3))
  allocate(force_upside(1:wing%ntri),force_downside(1:wing%ntri))
  allocate(centroid_upside_projected(1:wing%ntri,1:3), centroid_downside_projected(1:wing%ntri,1:3))
  allocate(upside_distributed_forces(1:wing%ntri,1:3,1:3), downside_distributed_forces(1:wing%ntri,1:3,1:3))

  np=wing%np

  !Getting surface position from the centerline of the wing
  upside = 0.d0
  downside = 0.d0

  call calculate_wing_surfaces(upside(1:wing%ntri,1:3,1:3), &
                               wing%u_new(1:np), &
                               wing%u_new(np+1:2*np), &
                               wing%u_new(2*np+1:3*np),   &
                               wing%t_wing, &
                               wing%tri_elements(1:wing%ntri,1:4), &
                               wing%tri_element_normals(1:wing%ntri,1:4),"upside")

  call calculate_wing_surfaces(downside(1:wing%ntri,1:3,1:3), &
                               wing%u_new(1:np), &
                               wing%u_new(np+1:2*np), &
                               wing%u_new(2*np+1:3*np),   &
                               wing%t_wing, &
                               wing%tri_elements(1:wing%ntri,1:4), &
                               wing%tri_element_normals(1:wing%ntri,1:4),"downside")


  force_upside=0.d0
  force_downside=0.d0
  centroid_upside_projected =0.d0
  centroid_downside_projected = 0.d0
  upside_distributed_forces = 0.d0
  downside_distributed_forces = 0.d0
  do itri=1,wing%ntri

  !Calculate pressure forces acting at the barycenter of the triangle elements
  force_upside(itri) = wing%tri_element_areas(itri)*sum(wing%press_upside(itri,1:3))/3
  force_downside(itri) = wing%tri_element_areas(itri)*sum(wing%press_downside(itri,1:3))/3

  !Calculate the centroid of the pressure truncated (since the pressure at the three
  !vertices are not the same) triangular prism element
  ! NOTE at first, for simplicity, the centroid here is calculated based directly on
  ! the three vertices of the triangle and the corresponding pressure values
  call truncated_triangular_prism_centroid(centroid_upside_projected(itri,1:3),&
                                            upside(itri,1,1:3), &
                                            upside(itri,2,1:3), &
                                            upside(itri,3,1:3), &
                                                                  wing%press_upside(itri,1), &
                                                                  wing%press_upside(itri,2), &
                                                                  wing%press_upside(itri,3))

  call truncated_triangular_prism_centroid(centroid_downside_projected(itri,1:3),downside(itri,1,1:3),   &
                                                                    downside(itri,2,1:3),   &
                                                                    downside(itri,3,1:3),   &
                                                                    wing%press_downside(itri,1), &
                                                                    wing%press_downside(itri,2), &
                                                                    wing%press_downside(itri,3))

  !
  !call project_point_onto_triangle(centroid_upside_projected,upside(itri,1,1:3),&
  !                                 centroid_upside(itri,1:3),wing%tri_element_normals(itri,1:3))
  !call project_point_onto_triangle(centroid_downside_projected,downside(itri,1,1:3),&
  !                                 centroid_downside(itri,1:3),wing%tri_element_normals(itri,1:3))

  call distribute_concentrated_force_into_three_vertices(upside_distributed_forces(itri,1:3,1:3),force_upside(itri),&
                                                          centroid_upside_projected(itri,1:3), &
                                                         (/wing%x(wing%tri_elements(itri,2)), &
                                                           wing%y(wing%tri_elements(itri,2)), &
                                                           wing%z(wing%tri_elements(itri,2))/), &
                                                         (/wing%x(wing%tri_elements(itri,3)), &
                                                           wing%y(wing%tri_elements(itri,3)), &
                                                           wing%z(wing%tri_elements(itri,3))/), &
                                                         (/wing%x(wing%tri_elements(itri,4)), &
                                                           wing%y(wing%tri_elements(itri,4)), &
                                                           wing%z(wing%tri_elements(itri,4))/), &
                                                           wing%tri_element_normals(itri,1:3), &
							                                             wing%tri_element_normals(itri,4))

  call distribute_concentrated_force_into_three_vertices(downside_distributed_forces(itri,1:3,1:3),force_downside(itri),&
                                                          centroid_downside_projected(itri,1:3), &
                                                         (/wing%x(wing%tri_elements(itri,2)), &
                                                           wing%y(wing%tri_elements(itri,2)), &
                                                           wing%z(wing%tri_elements(itri,2))/), &
                                                         (/wing%x(wing%tri_elements(itri,3)), &
                                                           wing%y(wing%tri_elements(itri,3)), &
                                                           wing%z(wing%tri_elements(itri,3))/), &
                                                         (/wing%x(wing%tri_elements(itri,4)), &
                                                           wing%y(wing%tri_elements(itri,4)), &
                                                           wing%z(wing%tri_elements(itri,4))/), &
                                                           wing%tri_element_normals(itri,1:3),&
							                                             wing%tri_element_normals(itri,4))

         wing%Fext(wing%tri_elements(itri,2)) = wing%Fext(wing%tri_elements(itri,2)) + &
                                                (downside_distributed_forces(itri,1,1) - upside_distributed_forces(itri,1,1))
    wing%Fext(wing%tri_elements(itri,2) + np) = wing%Fext(wing%tri_elements(itri,2) + np) + &
                                                (downside_distributed_forces(itri,1,2) - upside_distributed_forces(itri,1,2))
  wing%Fext(wing%tri_elements(itri,2) + 2*np) = wing%Fext(wing%tri_elements(itri,2) + 2*np) + &
                                                (downside_distributed_forces(itri,1,3) - upside_distributed_forces(itri,1,3))

         wing%Fext(wing%tri_elements(itri,3)) = wing%Fext(wing%tri_elements(itri,3)) + &
                                                (downside_distributed_forces(itri,2,1) - upside_distributed_forces(itri,2,1))
    wing%Fext(wing%tri_elements(itri,3) + np) = wing%Fext(wing%tri_elements(itri,3) + np) + &
                                                (downside_distributed_forces(itri,2,2) - upside_distributed_forces(itri,2,2))
  wing%Fext(wing%tri_elements(itri,3) + 2*np) = wing%Fext(wing%tri_elements(itri,3) + 2*np) + &
                                                (downside_distributed_forces(itri,2,3) - upside_distributed_forces(itri,2,3))

         wing%Fext(wing%tri_elements(itri,4)) = wing%Fext(wing%tri_elements(itri,4)) + &
                                                (downside_distributed_forces(itri,3,1) - upside_distributed_forces(itri,3,1))
    wing%Fext(wing%tri_elements(itri,4) + np) = wing%Fext(wing%tri_elements(itri,4) + np) + &
                                                (downside_distributed_forces(itri,3,2) - upside_distributed_forces(itri,3,2))
  wing%Fext(wing%tri_elements(itri,4) + 2*np) = wing%Fext(wing%tri_elements(itri,4) + 2*np) + &
                                                (downside_distributed_forces(itri,3,3) - upside_distributed_forces(itri,3,3))

  enddo

  deallocate(upside,downside)
  deallocate(force_upside,force_downside)
  deallocate(centroid_upside_projected, centroid_downside_projected)
  deallocate(upside_distributed_forces, downside_distributed_forces)

end subroutine

subroutine transform_pressure_into_point_forces_per_node(wing)

  implicit none

  type(flexible_wing), intent (inout) :: wing
  real(kind=pr), allocatable :: upside(:,:,:), downside(:,:,:)
  real(kind=pr), allocatable :: point_forces(:,:,:)
  integer :: itri,np

  allocate(upside(1:wing%ntri,1:3,1:3),downside(1:wing%ntri,1:3,1:3))
  allocate(point_forces(1:wing%ntri,1:3,1:3))

  np=wing%np

  !Getting surface position from the centerline of the wing
  upside = 0.d0
  downside = 0.d0

  call calculate_wing_surfaces(upside(1:wing%ntri,1:3,1:3), &
                               wing%u_new(1:np), &
                               wing%u_new(np+1:2*np), &
                               wing%u_new(2*np+1:3*np),   &
                               wing%t_wing, &
                               wing%tri_elements(1:wing%ntri,1:4), &
                               wing%tri_element_normals(1:wing%ntri,1:4),"upside")

  call calculate_wing_surfaces(downside(1:wing%ntri,1:3,1:3), &
                               wing%u_new(1:np), &
                               wing%u_new(np+1:2*np), &
                               wing%u_new(2*np+1:3*np),   &
                               wing%t_wing, &
                               wing%tri_elements(1:wing%ntri,1:4), &
                               wing%tri_element_normals(1:wing%ntri,1:4),"downside")


  point_forces = 0.d0

  do itri=1,wing%ntri

  !Calculate pressure forces acting at the barycenter of the triangle elements
  point_forces(itri,1,1:3) = wing%tri_element_areas(itri)/3* &
                            (wing%press_downside(itri,1)-wing%press_upside(itri,1))* &
                             wing%tri_element_normals(itri,4)*wing%tri_element_normals(itri,1:3)

  point_forces(itri,2,1:3) = wing%tri_element_areas(itri)/3* &
                            (wing%press_downside(itri,2)-wing%press_upside(itri,2))* &
                             wing%tri_element_normals(itri,4)*wing%tri_element_normals(itri,1:3)

  point_forces(itri,3,1:3) = wing%tri_element_areas(itri)/3* &
                            (wing%press_downside(itri,3)-wing%press_upside(itri,3))* &
                             wing%tri_element_normals(itri,4)*wing%tri_element_normals(itri,1:3)


         wing%Fext(wing%tri_elements(itri,2)) = wing%Fext(wing%tri_elements(itri,2)) + point_forces(itri,1,1)
    wing%Fext(wing%tri_elements(itri,2) + np) = wing%Fext(wing%tri_elements(itri,2) + np) + point_forces(itri,1,2)
  wing%Fext(wing%tri_elements(itri,2) + 2*np) = wing%Fext(wing%tri_elements(itri,2) + 2*np) + point_forces(itri,1,3)

         wing%Fext(wing%tri_elements(itri,3)) = wing%Fext(wing%tri_elements(itri,3)) + point_forces(itri,2,1)
    wing%Fext(wing%tri_elements(itri,3) + np) = wing%Fext(wing%tri_elements(itri,3) + np) + point_forces(itri,2,2)
  wing%Fext(wing%tri_elements(itri,3) + 2*np) = wing%Fext(wing%tri_elements(itri,3) + 2*np) + point_forces(itri,2,3)

         wing%Fext(wing%tri_elements(itri,4)) = wing%Fext(wing%tri_elements(itri,4)) + point_forces(itri,3,1)
    wing%Fext(wing%tri_elements(itri,4) + np) = wing%Fext(wing%tri_elements(itri,4) + np) + point_forces(itri,3,2)
  wing%Fext(wing%tri_elements(itri,4) + 2*np) = wing%Fext(wing%tri_elements(itri,4) + 2*np) + point_forces(itri,3,3)

  enddo

  deallocate(upside,downside)
  deallocate(point_forces)

end subroutine

subroutine distribute_concentrated_force_into_three_vertices(distributed_force,concentrated_force,&
  concentrated_force_position,tri1,tri2,tri3,normal,direction)
! Replace the concentrated force at the projection of the truncated-triangular-prism's centroid by a system
! of forces applied at the three vertices of the triangle. This is based on the principle of two equivalent
! force systems, these new three forces must apply the same force and moment on the triangle as the old
! concentrated force

  implicit none

  ! force vectors distributed at three vertices using the principle of two equivalent force systems
  real(kind=pr),dimension(1:3,1:3),intent(inout) :: distributed_force
  ! magnitude of the concentrated force
  real(kind=pr), intent(in) :: concentrated_force
  ! the application point of the concentrated_force (projection of the centroid on the triangle)
  real(kind=pr), dimension(1:3), intent(in) :: concentrated_force_position
  ! coordinates of the three vertices of the triangle
  real(kind=pr),dimension(1:3),intent(in) :: tri1, tri2, tri3
  ! the normal vector of the triangle, which will also be the direction of the forces
  real(kind=pr),dimension(1:3),intent(in) :: normal
  ! determine whether the force acting on the upside or the downside of the wing
  real(kind=pr),dimension(1:3) :: tmp
  real(kind=pr),intent(in) :: direction
  real(kind=pr) :: A1, A2, A3, A

  ! areas formed by three vertices and the application point of the concentrated_force for barycentric interpolation
  tmp = cross(tri2-concentrated_force_position,tri3-concentrated_force_position)
  A1 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1-concentrated_force_position,tri3-concentrated_force_position)
  A2 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1-concentrated_force_position,tri2-concentrated_force_position)
  A3 = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  tmp = cross(tri1,tri2)
  A = 0.5*sqrt(tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

  !Check if the projection of the centroid is inside the triangle
  !if (root) then
  !  if ((A1 + A2 + A3) > (A + 1.0d-15)) then
  !    write(*,*) "WARNING: The projection of the centroid may be outside of the triangle!"
  !    write(*,*) A1, A2, A3
  !    write(*,*) A
  !    write(*,*) tri1, tri2, tri3
  !    write(*,*) concentrated_force_position
  !  endif
  !endif

  distributed_force(1,1:3) = - A1/A*normal(1:3)*concentrated_force*direction
  distributed_force(2,1:3) = - A2/A*normal(1:3)*concentrated_force*direction
  distributed_force(3,1:3) = - A3/A*normal(1:3)*concentrated_force*direction

end subroutine


subroutine fictitious_forces_of_moving_reference_frame (time,dt0,dt1,it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), intent (inout) :: wings
integer :: i,j, np

  np = wings%np

  do j=1,np
    wings%Fext(j)        = wings%Fext(j) - wings%at_inertia(1)*wings%m(j)  !forces on the x-direction
    wings%Fext(j + np)   = wings%Fext(j + np) - wings%at_inertia(2)*wings%m(j) !forces on the y-direction
    wings%Fext(j + 2*np) = wings%Fext(j + 2*np) - wings%at_inertia(3)*wings%m(j) !forces on the z-direction
  enddo

  !do i=1,nVeins_BC
  !  j = nint(wings%veins_bending_BC(1,2,i))
  !  if (it==0) write(*,*) j
  !    wings%Fext(j)        = wings%Fext(j) - wings%at_inertia(1)*wings%m(j)  !forces on the x-direction
  !    wings%Fext(j + np)   = wings%Fext(j + np) - wings%at_inertia(2)*wings%m(j) !forces on the y-direction
  !    wings%Fext(j + 2*np) = wings%Fext(j + 2*np) - wings%at_inertia(3)*wings%m(j) !forces on the z-direction

  !enddo


end subroutine
