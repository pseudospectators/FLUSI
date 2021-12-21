subroutine length_calculation_wrapper(x,y,z,springs)
!This subroutine calculates actually for the whole vein or membrane

implicit none

real(kind=pr), intent(in) :: x(:),y(:),z(:)
real(kind=pr), intent(inout) :: springs(:,:)
integer :: ind

do ind=1,nint(maxval(springs(:,1)))

    call length_calculation(x(nint(springs(ind,2))), &
                            x(nint(springs(ind,3))), &
                            y(nint(springs(ind,2))), &
                            y(nint(springs(ind,3))), &
                            z(nint(springs(ind,2))), &
                            z(nint(springs(ind,3))), &
                            springs(ind,5))

enddo

end subroutine

subroutine length_calculation(x1,x2,y1,y2,z1,z2,length)

  real(kind=pr), intent(in) :: x1,x2,y1,y2,z1,z2
  real(kind=pr), intent(inout) :: length

  length = sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

end subroutine

subroutine angle_calculation_wrapper(x,y,z,vein)
!Calculate relative angles between two consecutive segments belonging to a vein
!   x(vein(ind,2))          x(vein(ind,3))                x(vein(ind,4))
!        O-----------------------O----------------------------O
!                  X1                          X2
!This subroutine calculates actually for the whole vein

implicit none

real(kind=pr), intent(in) :: x(:),y(:),z(:)
real(kind=pr), intent(inout) :: vein(:,:)
!real(kind=pr) :: x1,x2,y1,y2,z1,z2
integer :: ind

do ind=1,nint(maxval(vein(:,1)))

    call angle_calculation(x(nint(vein(ind,2))),x(nint(vein(ind,3))), &
                           x(nint(vein(ind,4))),y(nint(vein(ind,2))), &
                           y(nint(vein(ind,3))),y(nint(vein(ind,4))), &
                           z(nint(vein(ind,2))),z(nint(vein(ind,3))), &
                           z(nint(vein(ind,4))),vein(ind,7),vein(ind,8))

enddo

end subroutine

subroutine angle_calculation(x1,x2,x3,y1,y2,y3,z1,z2,z3,theta,phi)

implicit none

  real(kind=pr), intent(in) :: x1,x2,x3,y1,y2,y3,z1,z2,z3
  real(kind=pr), intent(inout) :: theta, phi
  real(kind=pr) :: dx1, dy1, dz1, dx2, dy2, dz2

  dx1 = x2 - x1
  dy1 = y2 - y1
  dz1 = z2 - z1
  dx2 = x3 - x2
  dy2 = y3 - y2
  dz2 = z3 - z2

 ! theta = atan2(dx1*dy2-dx2*dy1,dx1*dx2+dy1*dy2)

  !  phi = atan2(dx1*dz2-dx2*dz1,dx1*dx2+dz1*dz2)

!  if ((dx1<1.0d-16) .or. (dx2<1.0d-16)) then
 !   theta = 0.0d0
 ! else
    theta = atan2(dx1*dy2-dx2*dy1,dx1*dx2+dy1*dy2)
 ! endif
      phi = atan2(dx1*dz2-dx2*dz1,dx1*dx2+dz1*dz2)


end subroutine

subroutine set_diameters_for_veins(d_veins, d_veins_BC, &
  middle_point_indices, middle_point_indices_BC)

  implicit none
  real(kind=pr), intent(inout) :: d_veins(1:,1:), d_veins_BC(1:,1:)
  real(kind=pr), intent(inout) :: middle_point_indices(1:), middle_point_indices_BC(1:)

  d_veins_BC(1:3,1) = 200
  d_veins_BC(1:3,2) = 350
  d_veins_BC(1:3,3) = 250

  !vein 1
  d_veins(1:3,1) = 155
  d_veins(1,2)   = 150
  d_veins(2,2)   = 100
  d_veins(3,2)   = 145
  d_veins(1,3)   = 150
  d_veins(2,3)   = 115
  d_veins(3,3)   = 100
  d_veins(1,4)   = 100
  d_veins(2,4)   = 0
  d_veins(3,4)   = 80
  d_veins(1:3,5) = 80

  !vein 2
  d_veins(1:3,6) = 200
  d_veins(1:3,7) = 100

  !vein3
  d_veins(1:3,8) = 150
  d_veins(1:3,9) = 100

  !vein 4 middle_point_indices 89
  d_veins(1,10)  = 125
  d_veins(2,10)  = 0
  d_veins(3,10)  = 95
  d_veins(1,11)  = 85
  d_veins(2,11)  = 75
  d_veins(3,11)  = 60

  !vein 5 middle_point_indices 142 126
  d_veins(1,12)  = 120
  d_veins(2,12)  = 0
  d_veins(3,12)  = 150
  d_veins(1,13)  = 80
  d_veins(2,13)  = 0
  d_veins(3,13)  = 70
  d_veins(1,14)  = 70
  d_veins(2,14)  = 65
  d_veins(3,14)  = 60
  d_veins(1,15)  = 60
  d_veins(2,15)  = 55
  d_veins(3,15)  = 50

  !vein 5 middle_point_indices 204
  d_veins(1:3,16) = 120
  d_veins(1,17)   = 120
  d_veins(2,17)   = 0
  d_veins(3,17)   = 100
  d_veins(1,18)   = 120
  d_veins(2,18)   = 80
  d_veins(3,18)   = 55

  !vein 7
  d_veins(1:3,19) = 60

  !vein 8 middle_point_indices 351
  d_veins(1:3,20) = 120
  d_veins(1:3,21) = 120
  d_veins(1,22)   = 115
  d_veins(2,22)   = 50
  d_veins(3,22)   = 5

  !vein 9 middle_point_indices 192
  d_veins(1,23)   = 250
  d_veins(2,23)   = 200
  d_veins(3,23)   = 80

  !vein 10
  d_veins(1,24)   = 40
  d_veins(2,24)   = 0
  d_veins(3,24)   = 35

  !vein 11
  d_veins(1:3,25) = 130

  !vein 12
  d_veins(1:3,26) = 120

  !vein 13
  d_veins(1:3,27) = 150

  !vein 14, 15, 16 and 17
  d_veins(1:3,28) = 10
  d_veins(1:3,29) = 10
  d_veins(1:3,30) = 90
  d_veins(1:3,31) = 110

  middle_point_indices = (/0, 288, 295, 0, 0, 0, 0, 0, 0, 0, 89, 0, 0, 142, 126, &
                           195, 205, 0, 0, 0, 0, 351, 329, 0, 0, 0, 0, 0, 0, 0, 0/)

  middle_point_indices_BC = (/0, 0, 0/)

  d_veins = d_veins*4.55373406193078d-5
  d_veins_BC = d_veins_BC*4.55373406193078d-5
end subroutine

subroutine calculate_flexural_rigidity_from_Young_modulus(j,kby, kbz, E, d_veins, &
  middle_point_indices, veins_extension)

  implicit none
  integer, intent(in) :: j
  real(kind=pr), intent(in) :: E
  real(kind=pr), intent(in) :: d_veins(1:3)
  real(kind=pr), intent(in) :: veins_extension(1:,1:)
  real(kind=pr), intent(in) :: middle_point_indices
  real(kind=pr), intent(inout) :: kby(1:), kbz(1:)
  real(kind=pr) :: Lv, kb_head, kb_mid, kb_tail, l_head, l_mid, l_tail
  real(kind=pr), dimension(1:3) :: Iy, Iz, EIy, EIz
  integer :: i, nv, ind_mid
  logical :: stop_condition

  stop_condition = .true.

  !vein length
  Lv = sum(veins_extension(1:size(veins_extension,1),4))

  !second moment of inertia
  Iy=pi/2*(d_veins/2)**4
  Iz=pi/2*(d_veins/2)**4

  !flexural rigidity
  EIy=E*Iy
  EIz=E*Iz


    if (middle_point_indices == 0) then
        !vein with constant diameter
        if (d_veins(2)/=0) then
          do i=1,nint(maxval(veins_extension(:,1)))
            kby(i) = EIz(1)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                       (Lv*(2*maxval(veins_extension(:,1))+4))
            kbz(i) = EIy(1)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                       (Lv*(2*maxval(veins_extension(:,1))+4))
          enddo
        !vein whose diameter decreases or increases linearly
        elseif (d_veins(2)==0) then
            if ((EIz(2)/=0) .or. (EIz(1)==EIz(3))) then
                write(*,*) d_veins(1), d_veins(2), d_veins(3)
                write(*,*) j
                call abort(032020202,"Diameter data is wrong!")
            else
                kb_head = EIy(1)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                          (Lv*(2*maxval(veins_extension(:,1))+4))
                kb_tail = EIy(3)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                          (Lv*(2*maxval(veins_extension(:,1))+4))
                do i=1,nint(maxval(veins_extension(:,1)))
                    l_head = sum(veins_extension(1:i,4)) !distance from the current point to the head of the vein
                    l_tail = Lv - l_head; !distance from the current point to the tail of the vein
                    kby(i) = (l_tail*kb_head + l_head*kb_tail)/Lv
                    kbz(i) = kby(i)
                enddo
            endif
        endif
    elseif (middle_point_indices /= 0) then
                    kb_head = EIy(1)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                              (Lv*(2*maxval(veins_extension(:,1))+4))
                    kb_mid  = EIy(2)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                              (Lv*(2*maxval(veins_extension(:,1))+4))
                    kb_tail = EIy(3)*(maxval(veins_extension(:,1))+1)*(2*maxval(veins_extension(:,1))+3)/&
                              (Lv*(2*maxval(veins_extension(:,1))+4))
            i = 0
            do while (stop_condition)
              i = i + 1
              if (nint(veins_extension(i,3))==middle_point_indices) then
                ind_mid = i!index of the element corresponding to the measured middle point
                stop_condition = .false.
              endif
            enddo
            do i=1,ind_mid
                    l_head = sum(veins_extension(1:i,4)); !distance from the current point to the head of the vein
                    l_mid  = sum(veins_extension(1:ind_mid,4)) - l_head
                    kby(i) = (l_mid*kb_head + l_head*kb_mid)/sum(veins_extension(1:ind_mid,4))
                    kbz(i) = kby(i)
            enddo
            do i=ind_mid+1,nint(maxval(veins_extension(:,1)))
                    l_mid  = sum(veins_extension(ind_mid+1:i,4))
                    l_tail = Lv - sum(veins_extension(1:ind_mid,4)) - l_mid !distance from the current point to the tail of the vein
                    kby(i) = (l_mid*kb_tail + l_tail*kb_mid)/sum(veins_extension(ind_mid+1:nint(maxval(veins_extension(:,1))),4))
                    kbz(i) = kby(i)
            enddo
    endif


end subroutine

subroutine convert_flexural_rigidity_into_spring_stiffness(EIy, EIz, kby0, kbz0, extension_springs)

  implicit none
  real(kind=pr), intent(in) :: extension_springs(1:,:)
  real(kind=pr), intent(in) :: EIy, EIz
  real(kind=pr), intent(inout) :: kby0, kbz0
  real(kind=pr) :: L, ns

  ns = maxval(extension_springs(:,1)) + 1
  L = sum(extension_springs(1:(nint(ns)-1),4))

  kby0 = EIz/L*ns*(2*ns+1)/(2*ns+2)
  kbz0 = EIy/L*ns*(2*ns+1)/(2*ns+2)

end subroutine

subroutine solve_linear_system_wing ( A, b, x )
  implicit none
  real(kind=pr),dimension(1:,1:), intent(in) :: A
  real(kind=pr),dimension(1:), intent(inout) :: x
  real(kind=pr),dimension(1:), intent(in) :: b
  real(kind=pr) :: t0
  integer :: error, nn
  integer, allocatable, dimension(:) :: ipiv

  nn = size(b)

  if (size(b) /= size(x) ) call abort(91,'Solve_linear_system: size(b)/=size(x)')

  allocate(ipiv(1:nn))

  t0 = MPI_wtime()

  call dgetrf( nn, nn, A , nn, ipiv, error )
  if (error .ne. 0) then
    write(*,*) "!!! mmCrutial: dgetrf error.", error ,nn
    call abort(92, "Error in solve liner system (dgetrf)")
  endif

  x = b
  call dgetrs( 'N', nn, 1, A, nn, ipiv, x, nn, error )
  if (error .ne. 0) then
    write(*,*) "!!! mmCrutial: dgetrs error.", error ,nn
    call abort(93, "Error in solve liner system (dgetrs)")
  endif

  call toc("Flexible wing (solve_linear_system_wing)", MPI_wtime() - t0)

end subroutine

subroutine Moving_boundary_point(wings)


    implicit none
    type(flexible_wing), dimension(1:nWings), intent(inout) :: wings
    integer :: i,j

    do i=1,nWings
        do j=1,nVeins_BC

          wings(i)%x_BC(-1,j) = wings(i)%x_BC(-1,j) + 1.0d-2
          wings(i)%y_BC(-1,j) = wings(i)%y_BC(-1,j) + 1.0d-2
          wings(i)%z_BC(-1,j) = wings(i)%z_BC(-1,j) + 1.0d-2

          wings(i)%x_BC(0,j) = wings(i)%x_BC(-1,j)
          wings(i)%y_BC(0,j) = wings(i)%y_BC(-1,j)
          wings(i)%z_BC(0,j) = wings(i)%z_BC(-1,j)

        enddo
      enddo
end subroutine

subroutine rotate_wing(wing)

  implicit none
  !type(diptera), intent(inout) :: Insect
  type(flexible_wing), intent(inout) :: wing
  integer :: i
  real(kind=pr), dimension(1:3) :: u


  ! Rotate wing into the global coordinate system
  do i=1,wing%np
    u = matmul(wing%M_wing_inv,matmul(wing%M_solver_inv, &
                                  (/wing%u_old(i), wing%u_old(i+wing%np), wing%u_old(i+2*wing%np)/)))
    wing%x(i) = u(1)
    wing%y(i) = u(2)
    wing%z(i) = u(3)
  enddo


end subroutine

subroutine rotate_and_translate_wing_into_global_system(wing, Insect)

  implicit none
  type(diptera), intent(inout) :: Insect
  type(flexible_wing), intent(inout) :: wing
  integer :: i
  real(kind=pr), dimension(1:3) :: u


  ! Rotate wing into the global coordinate system
  do i=1,wing%np
    ! wing position in the body coordinate system
      ! Rotate the wing into the body coordinate system
      u = matmul(wing%M_wing_inv,matmul(wing%M_solver_inv, &
                                  (/wing%u_old(i), wing%u_old(i+wing%np), wing%u_old(i+2*wing%np)/)))

      ! Translate the wing to the pivot point
      wing%x(i) = u(1) + wing%x_pivot_b(1)
      wing%y(i) = u(2) + wing%x_pivot_b(2)
      wing%z(i) = u(3) + wing%x_pivot_b(3)

    ! wing position in the global coordinate system
      ! Rotate the wing into the global coordinate system
      u = matmul(Insect%M_body_inv,(/wing%x(i), wing%y(i), wing%z(i)/))

      ! Translate the wing to the body center
      wing%x(i) = u(1) + wing%x0(1)
      wing%y(i) = u(2) + wing%x0(2)
      wing%z(i) = u(3) + wing%x0(3)
  enddo

end subroutine


subroutine construct_total_velocity(wing,M_body,M_body_inv)

  implicit none
  real(kind=pr), intent(in) :: M_body(1:3,1:3), M_body_inv(1:3,1:3)
  type(flexible_wing), intent(inout) :: wing
  integer :: j, np
  real(kind=pr), dimension(1:3) :: u,v,r

  np = wing%np

  !if (root) write(*,*) wing%rot_rel_wing_g(1), wing%rot_rel_wing_g(2), wing%rot_rel_wing_g(3)
  !if (root) write(*,*) wing%rot_rel_wing_w(1), wing%rot_rel_wing_w(2), wing%rot_rel_wing_w(3)
  !if (root) write(*,*) wing%vr0(1), wing%vr0(2), wing%vr0(3)

  do j=1,np

    r(1) = wing%u_old(0*wing%np+j)
    r(2) = wing%u_old(1*wing%np+j)
    r(3) = wing%u_old(2*wing%np+j)

    call rotate_vector_into_global_system(wing,M_body,r)

    u = cross(wing%rot_rel_wing_g(1:3),r)

    v(1) = wing%u_old(3*wing%np+j)
    v(2) = wing%u_old(4*wing%np+j)
    v(3) = wing%u_old(5*wing%np+j)

    call rotate_vector_into_global_system(wing,M_body_inv,v)

    wing%vx(j) = wing%vt0(1) + v(1) + u(1)
    wing%vy(j) = wing%vt0(2) + v(2) + u(2)
    wing%vz(j) = wing%vt0(3) + v(3) + u(3)

  enddo

end subroutine

subroutine translate_wing(wing,M_body)

  implicit none
  real(kind=pr), intent(in) :: M_body(1:3,1:3)
  type(flexible_wing), intent(inout) :: wing
  integer :: j, np
  real(kind=pr), dimension(1:3) :: u,v

  np = wing%np

  !if (root) write(*,*) wing%rot_rel_wing_g(1), wing%rot_rel_wing_g(2), wing%rot_rel_wing_g(3)
  !if (root) write(*,*) wing%rot_rel_wing_w(1), wing%rot_rel_wing_w(2), wing%rot_rel_wing_w(3)
  !if (root) write(*,*) wing%vr0(1), wing%vr0(2), wing%vr0(3)

  do j=1,np

    u = cross(wing%rot_rel_wing_g(1:3),(/wing%x(j),wing%y(j),wing%z(j)/))


    v(1) = wing%u_old(3*wing%np+j)
    v(2) = wing%u_old(4*wing%np+j)
    v(3) = wing%u_old(5*wing%np+j)

    call rotate_vector_into_global_system(wing,M_body,v)

    wing%vx(j) = wing%vt0(1) + v(1) + u(1)
    wing%vy(j) = wing%vt0(2) + v(2) + u(2)
    wing%vz(j) = wing%vt0(3) + v(3) + u(3)

    wing%x(j) = wing%x(j) + wing%x_pivot_b(1)
    wing%y(j) = wing%y(j) + wing%x_pivot_b(2)
    wing%z(j) = wing%z(j) + wing%x_pivot_b(3)

  enddo

end subroutine

subroutine rotate_vector_into_wing_system(wing,M_body,Vector)

implicit none

!real(kind=pr),intent(in) :: time
real(kind=pr), intent(in) :: M_body(1:3,1:3)
type(flexible_wing), intent (in) :: wing
real(kind=pr), intent(inout) :: Vector(1:3)
real(kind=pr), dimension(1:3,1:3) :: Mx, My, Mz
real(kind=pr), dimension(1:3) :: u


!The wing is rotated based on conventional Euler angles. It is then rotated in
!order: around z axis (yaw) first, then y axis (pitch) and finally x axis (roll)

u = matmul(wing%M_solver,matmul(wing%M_wing,matmul(M_body,(/Vector(1), Vector(2), Vector(3)/))))
Vector(1) = u(1)
Vector(2) = u(2)
Vector(3) = u(3)


end subroutine

subroutine rotate_vector_into_global_system(wing,M_body_inv,Vector)

implicit none

!real(kind=pr),intent(in) :: time
real(kind=pr), intent(in) :: M_body_inv(1:3,1:3)
type(flexible_wing), intent (inout) :: wing
real(kind=pr), intent(inout) :: Vector(1:3)
real(kind=pr), dimension(1:3,1:3) :: Mx, My, Mz
real(kind=pr), dimension(1:3) :: u


!The wing is rotated based on conventional Euler angles. It is then rotated in
!order: around z axis (yaw) first, then y axis (pitch) and finally x axis (roll)

!  call Rx(Mx,wing%phi)
!  call Ry(My,wing%alpha)
!  call Rz(Mz,wing%theta)

! Rotate wing around x axis
  u = matmul(M_body_inv,matmul(wing%M_wing_inv,matmul(wing%M_solver_inv,(/Vector(1), Vector(2), Vector(3)/))))
  Vector(1) = u(1)
  Vector(2) = u(2)
  Vector(3) = u(3)

! Rotate wing around x axis
!  u = matmul(My,(/Vector(1), Vector(2), Vector(3)/))!
!  Vector(1) = u(1)
!  Vector(2) = u(2)
!  Vector(3) = u(3)

! Rotate wing around x axis
!  u = matmul(Mz,(/Vector(1), Vector(2), Vector(3)/))
!  Vector(1) = u(1)
!  Vector(2) = u(2)
!  Vector(3) = u(3)


end subroutine

subroutine flexible_wing_angular_velocities (time, Wing, Insect, M_body )

  implicit none

  real(kind=pr), intent(in) :: time
  !real(kind=pr), intent(in) :: M_wing(1:3,1:3)
  real(kind=pr), intent(in) :: M_body(1:3,1:3)
  type(diptera), intent(inout) :: Insect
  type(flexible_wing), intent (inout) :: wing

  real(kind=pr) :: eta_stroke
  real(kind=pr) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  !real(kind=rk) :: phi_l, alpha_l, theta_l, phi_dt_l, alpha_dt_l, theta_dt_l
  real(kind=pr), dimension(1:3) :: rot_alpha, rot_theta, rot_phi, &
  rot_r_alpha, rot_r_theta, rot_r_phi
  real(kind=pr), dimension(1:3,1:3) :: M_wing, M_wing_r, &
  M1_tmp, M2_tmp, M1, M2, M3, M1_r, M2_r, M3_r, &
  M_stroke, M_stroke_r

  phi      = Wing%phi
  alpha    = Wing%alpha
  theta    = Wing%theta
  phi_dt   = Wing%phi_dt
  alpha_dt = Wing%alpha_dt
  theta_dt = Wing%theta_dt

  !phi_l      = Insect%phi_l
  !alpha_l    = Insect%alpha_l
  !theta_l    = Insect%theta_l
  !phi_dt_l   = Insect%phi_dt_l
  !alpha_dt_l = Insect%alpha_dt_l
  !theta_dt_l = Insect%theta_dt_l

  eta_stroke = Insect%eta_stroke

  !-----------------------------------------------------------------------------
  ! define the rotation matrices to change between coordinate systems
  !-----------------------------------------------------------------------------
  if (wing%ID == "left") then
    call Ry(M1_tmp,eta_stroke)
    M_stroke = M1_tmp
  elseif (wing%ID == "right") then
    call Rx(M1_tmp,pi)
    call Ry(M2_tmp,eta_stroke)
    M_stroke = matmul(M1_tmp,M2_tmp)
  endif

  if (wing%ID == "left") then
    call Ry(M1,alpha)
    call Rz(M2,theta)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,phi)
    M_wing = matmul(M1,matmul(M2,matmul(M3,M_stroke)))
  elseif (wing%ID == "right") then
    call Ry(M1,-alpha)
    call Rz(M2,theta)   ! Order changed (Dmitry, 7 Nov 2013)
    call Rx(M3,-phi)
    M_wing = matmul(M1,matmul(M2,matmul(M3,M_stroke)))
  endif


  !-----------------------------------------------------------------------------
  ! angular velocity vectors (in wing system)
  !-----------------------------------------------------------------------------
  if (wing%ID == "left") then
    rot_alpha = (/ 0.0d0, alpha_dt, 0.0d0 /)
    rot_theta = (/ 0.0d0, 0.0d0, theta_dt /)
    rot_phi   = (/ phi_dt, 0.0d0, 0.0d0   /)
  elseif (wing%ID == "right") then
    rot_alpha = (/ 0.0d0, -alpha_dt, 0.0d0/)
    rot_theta = (/ 0.0d0, 0.0d0, theta_dt /)
    rot_phi   = (/ -phi_dt, 0.0d0, 0.0d0  /)
  endif

  ! in the wing coordinate system
  wing%rot_rel_wing_w = matmul(M_wing,matmul(transpose(M_stroke),matmul(transpose(M3), &
  rot_phi+matmul(transpose(M2),rot_theta+matmul(transpose(M1), &
  rot_alpha)))))
  !Insect%rot_rel_wing_r_w = matmul(M_wing_r,matmul(transpose(M_stroke_r),matmul(transpose(M3_r), &
  !rot_r_phi+matmul(transpose(M2_r),rot_r_theta+matmul(transpose(M1_r), &
  !rot_r_alpha)))))

  ! direct definition, equivalent to what is above.
  ! Insect%rot_rel_wing_l_w = (/phi_dt_l*cos(alpha_l)*cos(theta_l)-theta_dt_l*sin(alpha_l),&
  !   alpha_dt_l-phi_dt_l*sin(theta_l),&
  !   theta_dt_l*cos(alpha_l)+phi_dt_l*sin(alpha_l)*cos(theta_l)/)

  ! prior to the call of this routine, the routine body_angular_velocity has
  ! computed the body angular velocity (both g/b frames) so here now we can also
  ! compute global and absolute wing angular velocities.
  wing%rot_rel_wing_b = matmul( transpose(M_wing), wing%rot_rel_wing_w )
  !Insect%rot_rel_wing_r_b = matmul( transpose(M_wing_r), Insect%rot_rel_wing_r_w )

  Wing%rot_rel_wing_g = matmul( transpose(M_body), Wing%rot_rel_wing_b )
  !Insect%rot_rel_wing_r_g = matmul( transpose(M_body), Insect%rot_rel_wing_r_b )

  wing%rot_abs_wing_g = wing%rot_rel_wing_g !+ wing%rot_body_g
  !Insect%rot_abs_wing_r_g = Insect%rot_body_g + Insect%rot_rel_wing_r_g

  wing%vr0 = wing%rot_rel_wing_g

  call rotate_vector_into_wing_system(wing,Insect%M_body,wing%vr0)


  !if (Insect%wing_fsi == "yes") then
    !**********************************
    !** Wing fsi model               **
    !**********************************
    ! overwrite the left wing
  !  Insect%rot_rel_wing_l_w = Insect%STATE(18:20)
  !  Insect%rot_rel_wing_l_b = matmul( transpose(M_wing_l), Insect%rot_rel_wing_l_w )
  !  Insect%rot_rel_wing_l_g = matmul( transpose(M_body), Insect%rot_rel_wing_l_b )
    ! the last guy is actually unused, as we have non-rotating body
  !  Insect%rot_abs_wing_l_g = Insect%rot_body_g + Insect%rot_rel_wing_l_g
  !endif


end subroutine

subroutine flexible_wing_angular_accel( time, wing, Insect )
  implicit none
  real(kind=pr), intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  type(diptera), intent(inout) :: Insect

  real(kind=pr) :: M_body(1:3,1:3), rot_dt_wing_g(1:3), M_wing_r(1:3,1:3), M_wing_l(1:3,1:3)
  type(diptera) :: Insect2
  type(flexible_wing) :: wing2
  real(kind=pr) :: dt,t

  dt = 1.0d-8
  Insect2 = Insect
  Wing2 = Wing

  Wing%rot_dt_wing_w = 0.d0
  !Wing%rot_dt_wing_r_w = 0.d0

  ! fetch motion state at time+dt
  call BodyMotion (time+dt, Insect2)
  if (wing2%ID == "left") then
    call Flexible_wing_motions ( time+dt, wing2, Insect2%kine_wing_l )
  elseif (wing2%ID == "right") then
    call Flexible_wing_motions ( time+dt, wing2, Insect2%kine_wing_r )
  endif
  call StrokePlane (time+dt, Insect2)
  call body_rotation_matrix( Insect2, M_body )
  call flexible_wing_angular_velocities ( time+dt, Wing2, Insect2, M_body )
  ! the current state is already calculated and saved in Wing and Insect variables

  ! use one-sided finite differences to derive the absolute angular velocity with
  ! respect to time. NOte in older code versions, this was wrong, as we derived
  ! the ang. vel. in the wing coordinate system, which is a moving reference frame.

  ! now happily, the old results are still correct, as long as the body does not rotate
  ! see, e.g., this document https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjzl-XP6_LNAhWoC5oKHUdDCHwQFggoMAI&url=http%3A%2F%2Focw.mit.edu%2Fcourses%2Faeronautics-and-astronautics%2F16-07-dynamics-fall-2009%2Flecture-notes%2FMIT16_07F09_Lec08.pdf&usg=AFQjCNHzEB-n_NMm6K3J1eRpIaGnuKpW0Q&sig2=yEPNin3bL5DnWauNJk2hcw&bvm=bv.126993452,d.bGs&cad=rjt
  ! however, if the body moves, an additional term occurs, and this was indeed missing
  ! in previous results.
  Wing%rot_dt_wing_g = (Wing2%rot_rel_wing_g - Wing%rot_rel_wing_g) / dt
  Wing%rot_dt_wing_w = matmul(Wing%M_solver,matmul(Wing%M_wing,matmul(Insect%M_body, Wing%rot_dt_wing_g)))

  Wing%ar0 = Wing%rot_dt_wing_w

  !rot_dt_wing_g = (Insect2%rot_rel_wing_r_g - Insect%rot_rel_wing_r_g) / dt
  !Insect%rot_dt_wing_r_w = matmul(M_wing_r,matmul(M_body, rot_dt_wing_g))


  ! if (root) then
  !   write(*,*) "L new code", Insect%rot_dt_wing_l_w
  !   write(*,*) "L old code", (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt
  !   write(*,*) "rot_rel_wing_l_g", Insect%rot_rel_wing_l_g
  !   write(*,*) "rot_body_g", Insect%rot_body_g
  !   write(*,*) "rot_body_b", Insect%rot_body_b
  !   write(*,*) "rot_dt_body_g", (Insect2%rot_body_g-Insect%rot_body_g)/dt
  ! endif
  !
  ! if (root) then
  !   t = 0.d0
  !   open  (17,file='test.t',status='replace')
  !   do while (t<=6.d0)
  !     call FlappingMotion_left ( t, Insect)
  !     write (17,'(7(es15.8,1x))') t,  &
  !     Insect%alpha_l, Insect%phi_l, Insect%theta_l, &
  !     Insect%alpha_dt_l, Insect%phi_dt_l, Insect%theta_dt_l
  !
  !     t = t + 1.0d-5
  !   end do
  !   close (17)
  !   call abort(7)
  ! endif

  ! ! this is the OLD CODE. It is actually wrong, since it computes the time derivative
  ! ! in the moving refrence frame of the wing, which is not the actual angular acceleration
  ! Insect%rot_dt_wing_r_w = (Insect2%rot_rel_wing_r_w - Insect%rot_rel_wing_r_w)/dt
  ! Insect%rot_dt_wing_l_w = (Insect2%rot_rel_wing_l_w - Insect%rot_rel_wing_l_w)/dt


  !if (Insect%wing_fsi == "yes") then
    !**********************************
    !** Wing fsi model               **
    !**********************************
    ! overwrite the left wings acceleration. the right hand side 18,19,20 is the
    ! time derivative of the angular velocity, so the acceleration is readily available
  !  Insect%rot_dt_wing_l_w = Insect%RHS_THIS(18:20)
  !endif

end subroutine flexible_wing_angular_accel

!-----------------------------------------------------------------------------
! Return the rotation matrix for the mass-spring solver
!-----------------------------------------------------------------------------
subroutine MSM_solver_rotation_matrix( Wing, M_solver )
! Because the MSM solver solves the wing deformation in a coordinate system where
! the wing is aligned along the x-axis and the upper surface is upward, oriented
! by the z-axis

  implicit none

  type(flexible_wing), intent (inout) :: wing
  real(kind=pr),intent(out) :: M_solver(1:3,1:3)
  real(kind=pr),dimension(1:3,1:3) :: M1, M2

    call Rx(M1,pi/2)
    call Rz(M2,pi/2)
    M_solver = matmul(M1,M2)

end subroutine MSM_solver_rotation_matrix

!-----------------------------------------------------------------------------
! return the rotation matrix for the flexible wing
!-----------------------------------------------------------------------------
subroutine flexible_wing_rotation_matrix( Wing, Insect, M_wing )
  implicit none

  type(diptera),intent(inout) :: Insect
  type(flexible_wing), intent (inout) :: wing
  real(kind=pr),intent(out) :: M_wing(1:3,1:3)
  real(kind=pr),dimension(1:3,1:3) :: M1, M2, M3, M_stroke, M1_tmp, M2_tmp

  !if ( Insect%wing_fsi /= "yes" ) then
    ! we're not using the wing fsi solver, so the wings follow a prescribed
    ! motion and we can compute the rotation matrix from the angles
    if (wing%ID == "left") then
      call Ry(M1_tmp,Insect%eta_stroke)
      M_stroke = M1_tmp
    elseif (wing%ID == "right") then
      call Rx(M1_tmp,pi)
      call Ry(M2_tmp,Insect%eta_stroke)
      M_stroke = matmul(M1_tmp,M2_tmp)
    endif

    if (wing%ID == "left") then
      call Ry(M1,Wing%alpha)
      call Rz(M2,Wing%theta)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3,Wing%phi)
      M_wing = matmul(M1,matmul(M2,matmul(M3,M_stroke)))
    elseif (wing%ID == "right") then
      call Ry(M1,-Wing%alpha)
      call Rz(M2,Wing%theta)   ! Order changed (Dmitry, 7 Nov 2013)
      call Rx(M3,-Wing%phi)
      M_wing = matmul(M1,matmul(M2,matmul(M3,M_stroke)))
    endif

  !else
    !**********************************
    !** Wing fsi model               **
    !**********************************
    ! in the wing FSI case, a quaternion-based formulation is used to get the
    ! rotation matrix from the wing quaterion. note the wing quaternion are the
    ! entries 14,15,16,17 of the Insect%STATE vector
    ! entries 18,19,20 are the angular VELOCITY of the wing
  !  call rotation_matrix_from_quaternion( Insect%STATE(14:17), M_wing_l)
  !endif
end subroutine flexible_wing_rotation_matrix



logical function Vector_isNAN(f)
  real(kind=pr), intent(in) :: f(:)
  integer :: a, i
  a = size(f)
  Vector_isNAN = .false.
  do i=1,a
    if (.not.(f(i).eq.f(i))) then
      Vector_isNAN = .true.
      return
    endif
  enddo
end function

logical function Matrix_isNAN(f,i,j)
  real(kind=pr), intent(in) :: f(:,:)
  integer, intent(out) :: i, j

  Matrix_isNAN = .false.
  do i=1,size(f,DIM=1)
    do j=1,size(f,DIM=2)
      if (.not.(f(i,j).eq.f(i,j))) then
        Matrix_isNAN = .true.
        if (root) then
          write(*,*) i, j
        endif
        !return
      endif
    enddo
  enddo
end function

subroutine truncated_triangular_prism_centroid(centroid,tri1,tri2,tri3,press_tri1,press_tri2,press_tri3)
!This function calculates the centroid of a truncated triangular prism by simple
!barycentric interpolation
! NOTICE: this can also be done by calculating directly the centroid of the truncated
! formed by the triangle and the pressure at the three vertices by deviding it
! into an ordinary triangular prism "minus" a pyramid

  implicit none
  real(kind=pr), dimension(1:3), intent(inout) :: centroid
  real(kind=pr), dimension(1:3), intent(in) :: tri1,tri2,tri3
  !array contains the pressure at the three vertices of the triangle tri1, tri2, tri3
  real(kind=pr), intent(in) :: press_tri1,press_tri2,press_tri3
  !the normal vector of the triangle
!  real(kind=pr), dimension(1:3), intent(in) :: normal

  centroid(1:3) =(press_tri1*tri1(1:3) + &
                                             press_tri2*tri2(1:3) + &
                                             press_tri3*tri3(1:3))/ &
                                            (press_tri1 + press_tri2 + press_tri3)

end subroutine

real(kind=pr) function pressure_trilinear_interp(x0, dx, field, x_target, periodic)
  use vars, only: pr, abort, mpirank, ga, gb
  use mpi
  implicit none
  ! origin of array grid. Note this is a compatibility issue with WABBIT. In FLUSI
  ! we usually have only one grid start starts at 0,0,0 and then the variables ra(:) and rb(:)
  ! inidicate which part is on the mpirank. In WABBIT we have bocks, always 1:n,1:n,1:n but
  ! each with its own origin and spacing. For transferability, this routine is written in the
  ! general style. In FLUSI:
  !     x0=(/dble(ra(1))*dx, dble(ra(2))*dy, dble(ra(3))*dz/) and dx=(/dx,dy,dz/)
  ! or, with ghost nodes:
  !     x0=(/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/) and dx=(/dx,dy,dz/)
  real(kind=pr),dimension(1:3),intent(in) :: x0
  ! spacing of array grid
  real(kind=pr),dimension(1:3),intent(in) :: dx
  ! the point at which we wish to interpolate
  real(kind=pr),dimension(1:3),intent(in) :: x_target
  ! actual field. zero-based indexing.
  real(kind=pr),intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  ! assume periodicity of field or not?
  ! ATTENTION: this means we suppose the array FIELDS to be PERIODIC ON ITS OWN
  ! the global field may well be PERIODIC, but if you pass mpi-chunks to this routine
  ! you MUST set periodic=.false. even if the global field is periodic.
  logical, intent(in) :: periodic

  ! array bounds of the field array
  integer,dimension(1:3) :: lbounds, ubounds
  real(kind=pr) :: xd, yd, zd
  real(kind=pr) :: c00, c10, c01, c11, c0, c1
  integer :: ix, iy, iz, ix1, iy1, iz1

  if (ga(1)==0) then
    lbounds(1) = ga(1)
  else
    lbounds(1) = ga(1)! - 1
  endif

  if (ga(2)==0) then
    lbounds(2) = ga(2)
  else
    lbounds(2) = ga(2) !- 1
  endif

  if (ga(3)==0) then
    lbounds(3) = ga(3)
  else
    lbounds(3) = ga(3) !- 1
  endif

  ubounds = gb(1:3)!(/size(field,1), size(field,2), size(field,3)/) - 1

  !write(*,'("Im CPU ",i5," gives pressure in x-direction from",i5," to ",i5," &
  !          in y-direction from",i5," to ",i5," &
  !          and in z-direction from",i5," to ",i5,".")') mpirank, ga(1), gb(1), &
  !          ga(2), gb(2), ga(3), gb(3)

  ! indices of cube containing the target point, lower end
  ix = floor( (x_target(1)-x0(1))/dx(1) ) + ga(1)
  iy = floor( (x_target(2)-x0(2))/dx(2) ) + ga(2)
  iz = floor( (x_target(3)-x0(3))/dx(3) ) + ga(3)

  ! distance to lower point, normalized (0..1)
  xd = ( x_target(1)-(dble(ix)*dx(1)+x0(1)) ) / dx(1)
  yd = ( x_target(2)-(dble(iy)*dx(2)+x0(2)) ) / dx(2)
  zd = ( x_target(3)-(dble(iz)*dx(3)+x0(3)) ) / dx(3)

  ! if the point is not on the grid, return a very large, negative value
  pressure_trilinear_interp = -9.9d10

  if (periodic) then
    ! *** periodic case ***
    if ( (ix>=lbounds(1)).and.(ix<=ubounds(1)) ) then
      if ( (iy>=lbounds(2)).and.(iy<=ubounds(2)) ) then
        if ( (iz>=lbounds(3)).and.(iz<=ubounds(3)) ) then
          ix1 = ix+1
          iy1 = iy+1
          iz1 = iz+1

          ! periodization. note ix,iy,iz can be ubounds at most, so the next point
          ! ix+1 would be the first point again.
          if (ix1>ubounds(1)) ix1 = lbounds(1)
          if (iy1>ubounds(2)) iy1 = lbounds(2)
          if (iz1>ubounds(3)) iz1 = lbounds(3)

          c00 = field(ix,iy  ,iz )*(1.d0-xd)+field(ix1 ,iy  ,iz )*xd
          c10 = field(ix,iy1 ,iz )*(1.d0-xd)+field(ix1 ,iy1 ,iz )*xd
          c01 = field(ix,iy  ,iz1)*(1.d0-xd)+field(ix1 ,iy  ,iz1)*xd
          c11 = field(ix,iy1 ,iz1)*(1.d0-xd)+field(ix1 ,iy1 ,iz1)*xd

          c0 = c00*(1.d0-yd) + c10*yd
          c1 = c01*(1.d0-yd) + c11*yd

          pressure_trilinear_interp = c0*(1.d0-zd)+c1*zd
        endif
      endif
    endif
  else

    ! *** non-periodic case ***
    if ( (ix>=lbounds(1)).and.(ix<ubounds(1)) ) then
      if ( (iy>=lbounds(2)).and.(iy<ubounds(2)) ) then
        if ( (iz>=lbounds(3)).and.(iz<ubounds(3)) ) then
          ix1 = ix+1
          iy1 = iy+1
          iz1 = iz+1

          c00 = field(ix,iy  ,iz )*(1.d0-xd)+field(ix1 ,iy  ,iz )*xd
          c10 = field(ix,iy1 ,iz )*(1.d0-xd)+field(ix1 ,iy1 ,iz )*xd
          c01 = field(ix,iy  ,iz1)*(1.d0-xd)+field(ix1 ,iy  ,iz1)*xd
          c11 = field(ix,iy1 ,iz1)*(1.d0-xd)+field(ix1 ,iy1 ,iz1)*xd

          c0 = c00*(1.d0-yd) + c10*yd
          c1 = c01*(1.d0-yd) + c11*yd

          pressure_trilinear_interp = c0*(1.d0-zd)+c1*zd
        endif
      endif
    endif

  endif

  ! with periodic it can work, but non-periodic not
  if (size(field,1)<=1) then
    call abort(9997111,"linear interpolation for 2d simulations currently not implemented...have fun")
  end if
end function pressure_trilinear_interp

subroutine dense_to_HB_matrix(A, n, nnz, colptr, rowind, values)

  implicit none

  integer, intent(in) :: n
  real(kind=pr), intent(in) :: A(n,n)
  integer, intent(out) :: nnz
  integer, allocatable, intent(out) :: colptr(:)
  integer, allocatable, intent(out) :: rowind(:)
  real(kind=pr), allocatable, intent(out) :: values(:)

  integer :: i, j, p

  if (allocated(rowind) .or. allocated(values)) then
     print *, "rowind and values must be unallocated on entro to dense_to_sparse()"
     stop 1
  endif

  nnz = count((A > 1.0d-15) .or. (A < -1.0d-15))
!   nnz = count(A /= 0.0d0)
 ! write(*,*) nnz
  allocate(colptr(1:n+1))
  allocate(rowind(1:nnz))
  allocate(values(1:nnz))

  p = 1 ! Index into the rowind and values arrays (1 based)
  do j=1,n ! Iterate over each column
     colptr(j) = p ! Save the starting index for this column in 0 based format
     do i=1,n ! Iterate over the rows, looking for non zero entries
        if ((A(i,j) > 1.0d-15) .or. (A(i,j) < -1.0d-15)) then
           rowind(p) = i  ! Zero based
           values(p) = A(i,j)
           p = p + 1
        endif
     end do
  end do
  colptr(n+1) = p ! Save the final index in 0 based format
!  write(*,*) colptr(n+1), nnz
  ! Check that the final index points to the nz'th entry
  if (colptr(n+1) /= (nnz + 1)) then
     write(*,*) "Something went wrong"
     write(*,*) colptr(n+1), nnz
  endif

end subroutine
