!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct internal force vector consists of all elastic forces coming from
! extension and bending springs exerted to mass points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine internal_forces_construction(Wings)
! This is actually just for 1 wing

implicit none

type(flexible_wing), intent(inout)  :: Wings
integer :: j, np, ind

! Get the number of mass points for the sake of simplicity in coding
np = wings%np

! Initialize
wings%Fint = 0.d0

  !Calculate current lenghts and angles of extension and bending springs repsectively
  do j=1,nMembranes
      call length_calculation_wrapper(wings%u_new(1:np), &
                              wings%u_new(np+1:2*np), &
                              wings%u_new(2*np+1:3*np),    &
                              wings%membranes_extension(:,:,j))

  enddo

  do j=1,nMembrane_edges
      call length_calculation_wrapper(wings%u_new(1:np), &
                              wings%u_new(np+1:2*np), &
                              wings%u_new(2*np+1:3*np), &
                              wings%membrane_edge(:,:,j))
  enddo

  do j=1,nVeins

      call length_calculation_wrapper(wings%u_new(1:np), &
                              wings%u_new(np+1:2*np), &
                              wings%u_new(2*np+1:3*np), &
                              wings%veins_extension(:,:,j))

      call angle_calculation_wrapper(wings%u_new(1:np), &
                             wings%u_new(np+1:2*np), &
                             wings%u_new(2*np+1:3*np), &
                             wings%veins_bending(:,:,j))
  enddo

  do j=1,nVeins_BC
      call length_calculation_wrapper(wings%u_new(1:np), &
                              wings%u_new(np+1:2*np), &
                              wings%u_new(2*np+1:3*np),   &
                              wings%veins_extension_BC(1:,:,j))
      call angle_calculation_wrapper(wings%u_new(1:np), &
                             wings%u_new(np+1:2*np), &
                             wings%u_new(2*np+1:3*np), &
                             wings%veins_bending_BC(1:,:,j))
  end do

  call angle_calculation_wrapper(wings%u_new(1:np), &
                          wings%u_new(np+1:2*np), &
                          wings%u_new(2*np+1:3*np), &
                          wings%vein_connectors(:,:))

  !Construct internal force vector
  do j=1,nMembranes
       call internal_extension_force(wings%Fint, &
                                     wings%u_new(1:np), &
                                     wings%u_new(np+1:2*np), &
                                     wings%u_new(2*np+1:3*np), &
                                     wings%membranes_extension(:,:,j), &
                                     np, &
                                     wings%ke_m(:,j))
  enddo

  do j=1,nMembrane_edges
  call internal_extension_force(wings%Fint, &
                                wings%u_new(1:np), &
                                wings%u_new(np+1:2*np), &
                                wings%u_new(2*np+1:3*np), &
                                wings%membrane_edge(:,:,j), &
                                np, &
                                wings%ke_me(:,j))
  enddo

  do j=1,nVeins
      call internal_extension_force(wings%Fint, &
                                    wings%u_new(1:np), &
                                    wings%u_new(np+1:2*np), &
                                    wings%u_new(2*np+1:3*np), &
                                    wings%veins_extension(:,:,j), np, &
                                    wings%ke_v(:,j))
      call internal_bending_force(wings%Fint, &
                                  wings%u_new(1:np),wings%u_new(np+1:2*np),wings%u_new(2*np+1:3*np), &
                                  wings%veins_bending(:,:,j), np, &
                                  wings%kby(:,j),wings%kbz(:,j))
  enddo

  do j=1,nVeins_BC

      call internal_extension_force_BC(wings%Fint, &
                                       wings%u_new(1:np), &
                                       wings%u_new(np+1:2*np),&
                                       wings%u_new(2*np+1:3*np), &
                                       wings%x_BC(0,j), wings%y_BC(0,j), wings%z_BC(0,j), &
                                       wings%veins_extension_BC(0:,:,j), np, &
                                       wings%ke_vBC(0:,j))

      call internal_bending_force_BC(wings%Fint, &
                                     wings%u_new(1:np),&
                                     wings%u_new(np+1:2*np),&
                                     wings%u_new(2*np+1:3*np), &
                                     wings%x_BC(-1:0,j), wings%y_BC(-1:0,j), wings%z_BC(-1:0,j), &
                                     wings%veins_bending_BC(-1:,:,j), np, &
                                     wings%kby_BC(-1:,j),wings%kbz_BC(-1:,j))

  end do

  call internal_bending_force(wings%Fint, &
                              wings%u_new(1:np),wings%u_new(np+1:2*np),wings%u_new(2*np+1:3*np), &
                              wings%vein_connectors(1:,1:), np, &
                              wings%kby_c,wings%kbz_c)

end subroutine

subroutine internal_extension_force(Fint,x,y,z,extension_springs,np,ke)

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: ke(:)
  real(kind=pr), intent(in)     :: extension_springs(:,:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: Fint(:)
  real(kind=pr) :: Fx, Fy, Fz
  integer :: ind


  do ind=1,nint(maxval(extension_springs(:,1)))

      call calculate_extension_spring_force( x(nint(extension_springs(ind,2))), &
                                             x(nint(extension_springs(ind,3))), &
                                             y(nint(extension_springs(ind,2))), &
                                             y(nint(extension_springs(ind,3))), &
                                             z(nint(extension_springs(ind,2))), &
                                             z(nint(extension_springs(ind,3))), &
                                             extension_springs(ind,4), &
                                             extension_springs(ind,5), &
                                             ke(ind), Fx, Fy, Fz)

             Fint(nint(extension_springs(ind,2))) = Fint(nint(extension_springs(ind,2))) + Fx

        Fint(nint(extension_springs(ind,2)) + np) = Fint(nint(extension_springs(ind,2)) + np) + Fy

      Fint(nint(extension_springs(ind,2)) + 2*np) = Fint(nint(extension_springs(ind,2)) + 2*np) + Fz

             Fint(nint(extension_springs(ind,3))) = Fint(nint(extension_springs(ind,3))) - Fx

        Fint(nint(extension_springs(ind,3)) + np) = Fint(nint(extension_springs(ind,3)) + np) - Fy

      Fint(nint(extension_springs(ind,3)) + 2*np) = Fint(nint(extension_springs(ind,3)) + 2*np) - Fz

  enddo

end subroutine

subroutine internal_extension_force_BC(Fint,x,y,z,x_BC,y_BC,z_BC,extension_springs_BC,np,ke_BC)

  real(kind=pr), intent(in)     :: x(1:),y(1:),z(1:)
  real(kind=pr), intent(in)     :: x_BC, y_BC, z_BC
  real(kind=pr), intent(in)     :: ke_BC(0:)
  real(kind=pr), intent(inout)     :: extension_springs_BC(0:,1:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: Fint(1:)
  real(kind=pr) :: Fx, Fy, Fz
  integer :: ind

  !Initialize
  Fx=0
  Fy=0
  Fz=0

  ! Calculate forces caused by springs inside the vein
  do ind=1,nint(maxval(extension_springs_BC(:,1)))

      call calculate_extension_spring_force( x(nint(extension_springs_BC(ind,2))), &
                                             x(nint(extension_springs_BC(ind,3))), &
                                             y(nint(extension_springs_BC(ind,2))), &
                                             y(nint(extension_springs_BC(ind,3))), &
                                             z(nint(extension_springs_BC(ind,2))), &
                                             z(nint(extension_springs_BC(ind,3))), &
                                             extension_springs_BC(ind,4), &
                                             extension_springs_BC(ind,5), &
                                             ke_BC(ind), Fx, Fy, Fz)

             Fint(nint(extension_springs_BC(ind,2))) = Fint(nint(extension_springs_BC(ind,2))) + Fx

        Fint(nint(extension_springs_BC(ind,2)) + np) = Fint(nint(extension_springs_BC(ind,2)) + np) + Fy

      Fint(nint(extension_springs_BC(ind,2)) + 2*np) = Fint(nint(extension_springs_BC(ind,2)) + 2*np) + Fz

             Fint(nint(extension_springs_BC(ind,3))) = Fint(nint(extension_springs_BC(ind,3))) - Fx

        Fint(nint(extension_springs_BC(ind,3)) + np) = Fint(nint(extension_springs_BC(ind,3)) + np) - Fy

      Fint(nint(extension_springs_BC(ind,3)) + 2*np) = Fint(nint(extension_springs_BC(ind,3)) + 2*np) - Fz

  enddo

  !---------------------------------------------------------------------------------
  ! Calculate forces caused by extension spring connecting the vein with the boundary
  !---------------------------------------------------------------------------------
  ! Vein schematic: point -1 and 0 are boundary points and they are clamped or
  ! rotated. Extension spring (numbered -1 and 0) connects points 0 and 1
  !
  !  -1       0         1             2       3
  !             ke_BC(0)
  !   X-------X---------O-------------O-------O .......
  !             l0(0)
  !

  ! For extension spring number 0 connecting two points 0 and 1
      call length_calculation(x_BC,x(nint(extension_springs_BC(1,2))), &
                              y_BC,y(nint(extension_springs_BC(1,2))), &
                              z_BC,z(nint(extension_springs_BC(1,2))), &
                              extension_springs_BC(0,5))


      ! Calculate force of the extension spring 0 acting on point 1
      call calculate_extension_spring_force(x_BC,x(nint(extension_springs_BC(1,2))), &
                                            y_BC,y(nint(extension_springs_BC(1,2))), &
                                            z_BC,z(nint(extension_springs_BC(1,2))), &
                                            extension_springs_BC(0,4),extension_springs_BC(0,5), &
                                            ke_BC(0), Fx, Fy, Fz)

       Fint(nint(extension_springs_BC(1,2))) = Fint(nint(extension_springs_BC(1,2))) - Fx

  Fint(nint(extension_springs_BC(1,2)) + np) = Fint(nint(extension_springs_BC(1,2)) + np) - Fy

Fint(nint(extension_springs_BC(1,2)) + 2*np) = Fint(nint(extension_springs_BC(1,2)) + 2*np) - Fz

end subroutine

subroutine internal_bending_force(Fint,x,y,z,bending_springs,np,kby,kbz)

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: kby(:), kbz(:)
  real(kind=pr), intent(in)     :: bending_springs(:,:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: Fint(:)
  real(kind=pr) :: Fx, Fy, Fz
  integer :: ind

  do ind=1,nint(maxval(bending_springs(:,1)))

    !Construct force vector for the left part of the bending spring
    call calculate_bending_spring_force(x(nint(bending_springs(ind,2))), &
                                        x(nint(bending_springs(ind,3))), &
                                        y(nint(bending_springs(ind,2))), &
                                        y(nint(bending_springs(ind,3))), &
                                        z(nint(bending_springs(ind,2))), &
                                        z(nint(bending_springs(ind,3))), &
                                        bending_springs(ind,5), &
                                        bending_springs(ind,6), &
                                        bending_springs(ind,7), &
                                        bending_springs(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           Fint(nint(bending_springs(ind,2))) = Fint(nint(bending_springs(ind,2))) + Fx

      Fint(nint(bending_springs(ind,2)) + np) = Fint(nint(bending_springs(ind,2)) + np) + Fy

    Fint(nint(bending_springs(ind,2)) + 2*np) = Fint(nint(bending_springs(ind,2)) + 2*np) + Fz

           Fint(nint(bending_springs(ind,3))) = Fint(nint(bending_springs(ind,3))) - Fx

      Fint(nint(bending_springs(ind,3)) + np) = Fint(nint(bending_springs(ind,3)) + np) - Fy

    Fint(nint(bending_springs(ind,3)) + 2*np) = Fint(nint(bending_springs(ind,3)) + 2*np) - Fz

    !Construct force vector for the right part of the bending spring
    call calculate_bending_spring_force(x(nint(bending_springs(ind,3))), &
                                        x(nint(bending_springs(ind,4))), &
                                        y(nint(bending_springs(ind,3))), &
                                        y(nint(bending_springs(ind,4))), &
                                        z(nint(bending_springs(ind,3))), &
                                        z(nint(bending_springs(ind,4))), &
                                        bending_springs(ind,5), &
                                        bending_springs(ind,6), &
                                        bending_springs(ind,7), &
                                        bending_springs(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           Fint(nint(bending_springs(ind,3))) = Fint(nint(bending_springs(ind,3))) - Fx

      Fint(nint(bending_springs(ind,3)) + np) = Fint(nint(bending_springs(ind,3)) + np) - Fy

    Fint(nint(bending_springs(ind,3)) + 2*np) = Fint(nint(bending_springs(ind,3)) + 2*np) - Fz

           Fint(nint(bending_springs(ind,4))) = Fint(nint(bending_springs(ind,4))) + Fx

      Fint(nint(bending_springs(ind,4)) + np) = Fint(nint(bending_springs(ind,4)) + np) + Fy

    Fint(nint(bending_springs(ind,4)) + 2*np) = Fint(nint(bending_springs(ind,4)) + 2*np) + Fz

  enddo

end subroutine

subroutine internal_bending_force_BC(Fint,x,y,z,x_BC,y_BC,z_BC,&
                                     bending_springs_BC,np,kby_BC,kbz_BC)

implicit none

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: x_BC(-1:),y_BC(-1:),z_BC(-1:)
  real(kind=pr), intent(in)     :: kby_BC(-1:), kbz_BC(-1:)
  real(kind=pr), intent(inout)  :: bending_springs_BC(-1:,:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: Fint(:)
  real(kind=pr) :: Fx, Fy, Fz
  integer :: ind

  do ind=1,nint(maxval(bending_springs_BC(:,1)))

    !Construct force vector for the left part of the bending spring
    call calculate_bending_spring_force(x(nint(bending_springs_BC(ind,2))), &
                                        x(nint(bending_springs_BC(ind,3))), &
                                        y(nint(bending_springs_BC(ind,2))), &
                                        y(nint(bending_springs_BC(ind,3))), &
                                        z(nint(bending_springs_BC(ind,2))), &
                                        z(nint(bending_springs_BC(ind,3))), &
                                        bending_springs_BC(ind,5), &
                                        bending_springs_BC(ind,6), &
                                        bending_springs_BC(ind,7), &
                                        bending_springs_BC(ind,8), &
                                        kby_BC(ind), kbz_BC(ind), Fx, Fy, Fz)

           Fint(nint(bending_springs_BC(ind,2))) = Fint(nint(bending_springs_BC(ind,2))) + Fx

      Fint(nint(bending_springs_BC(ind,2)) + np) = Fint(nint(bending_springs_BC(ind,2)) + np) + Fy

    Fint(nint(bending_springs_BC(ind,2)) + 2*np) = Fint(nint(bending_springs_BC(ind,2)) + 2*np) + Fz

           Fint(nint(bending_springs_BC(ind,3))) = Fint(nint(bending_springs_BC(ind,3))) - Fx

      Fint(nint(bending_springs_BC(ind,3)) + np) = Fint(nint(bending_springs_BC(ind,3)) + np) - Fy

    Fint(nint(bending_springs_BC(ind,3)) + 2*np) = Fint(nint(bending_springs_BC(ind,3)) + 2*np) - Fz

    !Construct force vector for the right part of the bending spring
    call calculate_bending_spring_force(x(nint(bending_springs_BC(ind,3))), &
                                        x(nint(bending_springs_BC(ind,4))), &
                                        y(nint(bending_springs_BC(ind,3))), &
                                        y(nint(bending_springs_BC(ind,4))), &
                                        z(nint(bending_springs_BC(ind,3))), &
                                        z(nint(bending_springs_BC(ind,4))), &
                                        bending_springs_BC(ind,5), &
                                        bending_springs_BC(ind,6), &
                                        bending_springs_BC(ind,7), &
                                        bending_springs_BC(ind,8), &
                                        kby_BC(ind), kbz_BC(ind), Fx, Fy, Fz)

           Fint(nint(bending_springs_BC(ind,3))) = Fint(nint(bending_springs_BC(ind,3))) - Fx

      Fint(nint(bending_springs_BC(ind,3)) + np) = Fint(nint(bending_springs_BC(ind,3)) + np) - Fy

    Fint(nint(bending_springs_BC(ind,3)) + 2*np) = Fint(nint(bending_springs_BC(ind,3)) + 2*np) - Fz

           Fint(nint(bending_springs_BC(ind,4))) = Fint(nint(bending_springs_BC(ind,4))) + Fx

      Fint(nint(bending_springs_BC(ind,4)) + np) = Fint(nint(bending_springs_BC(ind,4)) + np) + Fy

    Fint(nint(bending_springs_BC(ind,4)) + 2*np) = Fint(nint(bending_springs_BC(ind,4)) + 2*np) + Fz

  enddo
  !---------------------------------------------------------------------------------
  ! Calculate forces caused by bending springs connecting the vein with the boundary
  !---------------------------------------------------------------------------------
  ! Vein schematic: point -1 and 0 are boundary points and they are clamped or
  ! rotated. Two bending springs (numbered -1 and 0) are placed at points 0 and 1
  !
  !  -1       0         1             2       3 .......
  !        kb_BC(-1)  kb_BC(0)
  !   X-------X---------O-------------O-------O .......
  !
  !       theta(-1)   theta(0)
  !        phi(-1)     phi(0)
  !

  ! For bending spring number -1 connecting three points -1,0,1
      call angle_calculation(x_BC(-1),x_BC(0),x(nint(bending_springs_BC(1,2))), &
                             y_BC(-1),y_BC(0),y(nint(bending_springs_BC(1,2))), &
                             z_BC(-1),z_BC(0),z(nint(bending_springs_BC(1,2))), &
                             bending_springs_BC(-1,7),bending_springs_BC(-1,8))


      ! Calculate force of right part of the bending spring -1 acting on point 1
      call calculate_bending_spring_force(x_BC(0),x(nint(bending_springs_BC(1,2))), &
                                          y_BC(0),y(nint(bending_springs_BC(1,2))), &
                                          z_BC(0),z(nint(bending_springs_BC(1,2))), &
                                          bending_springs_BC(-1,5),bending_springs_BC(-1,6), &
                                          bending_springs_BC(-1,7),bending_springs_BC(-1,8), &
                                          kby_BC(-1), kbz_BC(-1), Fx, Fy, Fz)


          Fint(nint(bending_springs_BC(1,2))) = Fint(nint(bending_springs_BC(1,2))) + Fx

     Fint(nint(bending_springs_BC(1,2)) + np) = Fint(nint(bending_springs_BC(1,2)) + np) + Fy

   Fint(nint(bending_springs_BC(1,2)) + 2*np) = Fint(nint(bending_springs_BC(1,2)) + 2*np) + Fz

  ! For bending spring number 0 connecting three points 0,1,2
      call angle_calculation(x_BC(0),x(nint(bending_springs_BC(1,2))),x(nint(bending_springs_BC(1,3))), &
                            y_BC(0),y(nint(bending_springs_BC(1,2))),y(nint(bending_springs_BC(1,3))), &
                            z_BC(0),z(nint(bending_springs_BC(1,2))),z(nint(bending_springs_BC(1,3))), &
                            bending_springs_BC(0,7),bending_springs_BC(0,8))


      ! Calculate force of left part of the bending spring 0 acting on point 1
      call calculate_bending_spring_force(x_BC(0),x(nint(bending_springs_BC(1,2))), &
                                          y_BC(0),y(nint(bending_springs_BC(1,2))), &
                                          z_BC(0),z(nint(bending_springs_BC(1,2))), &
                                          bending_springs_BC(0,5),bending_springs_BC(0,6), &
                                          bending_springs_BC(0,7),bending_springs_BC(0,8), &
                                          kby_BC(0), kbz_BC(0), Fx, Fy, Fz)

          Fint(nint(bending_springs_BC(1,2))) = Fint(nint(bending_springs_BC(1,2))) - Fx

     Fint(nint(bending_springs_BC(1,2)) + np) = Fint(nint(bending_springs_BC(1,2)) + np) - Fy

   Fint(nint(bending_springs_BC(1,2)) + 2*np) = Fint(nint(bending_springs_BC(1,2)) + 2*np) - Fz

      ! Calculate force of right part of the bending spring 0 acting on point 1
      call calculate_bending_spring_force(x(nint(bending_springs_BC(1,2))), x(nint(bending_springs_BC(1,3))), &
                                          y(nint(bending_springs_BC(1,2))), y(nint(bending_springs_BC(1,3))), &
                                          z(nint(bending_springs_BC(1,2))), z(nint(bending_springs_BC(1,3))), &
                                          bending_springs_BC(0,5),bending_springs_BC(0,6), &
                                          bending_springs_BC(0,7),bending_springs_BC(0,8), &
                                          kby_BC(0), kbz_BC(0), Fx, Fy, Fz)


          Fint(nint(bending_springs_BC(1,2))) = Fint(nint(bending_springs_BC(1,2))) - Fx

     Fint(nint(bending_springs_BC(1,2)) + np) = Fint(nint(bending_springs_BC(1,2)) + np) - Fy

   Fint(nint(bending_springs_BC(1,2)) + 2*np) = Fint(nint(bending_springs_BC(1,2)) + 2*np) - Fz

      ! Calculate force of right part of the bending spring 0 acting on point 2

          Fint(nint(bending_springs_BC(1,3))) = Fint(nint(bending_springs_BC(1,3))) + Fx

     Fint(nint(bending_springs_BC(1,3)) + np) = Fint(nint(bending_springs_BC(1,3)) + np) + Fy

   Fint(nint(bending_springs_BC(1,3)) + 2*np) = Fint(nint(bending_springs_BC(1,3)) + 2*np) + Fz

end subroutine

subroutine calculate_extension_spring_force( x1, x2, y1, y2, z1, z2, &
  l0, l, ke, Fx, Fy, Fz )

real(kind=pr), intent(in)  :: x1, x2, y1, y2, z1, z2
real(kind=pr), intent(in)  :: l0, l
real(kind=pr), intent(in)  :: ke
real(kind=pr), intent(out) :: Fx, Fy, Fz


Fx = ke*(l-l0)*(x2-x1)/l

Fy = ke*(l-l0)*(y2-y1)/l

Fz = ke*(l-l0)*(z2-z1)/l

end subroutine

subroutine calculate_bending_spring_force( x1, x2, y1, y2, z1, z2, &
  theta0, phi0, theta, phi, kby, kbz, Fx, Fy, Fz )

real(kind=pr), intent(in)  :: x1, x2, y1, y2, z1, z2
real(kind=pr), intent(in)  :: theta0, theta, phi0, phi
real(kind=pr), intent(in)  :: kby, kbz
real(kind=pr), intent(out) :: Fx, Fy, Fz
real(kind=pr) :: idxy, idxz

if (((x2-x1)**2+(y2-y1)**2)<1.0d-12) then
    idxy=0
else
    idxy = 1/((x2-x1)**2+(y2-y1)**2)
endif

if (((x2-x1)**2+(z2-z1)**2)<1.0d-12) then
    idxz=0
else
    idxz = 1/((x2-x1)**2+(z2-z1)**2)
endif

Fx = kby*(theta - theta0)*(y2-y1)*idxy + kbz*(phi - phi0)*(z2-z1)*idxz

Fy = - kby*(theta - theta0)*(x2-x1)*idxy

Fz = - kbz*(phi - phi0)*(x2-x1)*idxz

end subroutine
