!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct internal force vector consists of all elastic forces coming from
! extension and bending springs exerted to mass points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine internal_forces_construction(Wings)

type(Wing), dimension(1:nWings), intent(inout)  :: Wings
integer :: i,j

!Calculate current lenghts and angles
do j=1:nMembranes
    call length_calculation(wings(i)%x,wings(i)%y,wings(i)%z,    &
                            wings(i)%membranes_extension(:,:,j), &
                            wings(i)%membranes_extension(:,5,j))
enddo

call length_calculation(wings(i)%x,wings(i)%y,wings(i)%z, &
                        wings(i)%membrane_edge,           &
                        wings(i)%membrane_edge(:,5))

do j=1:nVeins

    call length_calculation(wings(i)%x,wings(i)%y,wings(i)%z, &
                            wings(i)%veins_extension(:,:,j),  &
                            wings(i)%veins_extension(:,5,j))

    call angle_calculation(wings(i)%x,wings(i)%y,wings(i)%z, &
                           wings(i)%veins_bending(:,:,j),    &
                           wings(i)%veins_bending(:,7,j),    &
                           wings(i)%veins_bending(:,8,j))
enddo

do j=1:nVeins_BC
    call length_calculation(wings(i)%x,wings(i)%y,wings(i)%z,   &
                            wings(i)%veins_extension_BC(:,:,j), &
                            wings(i)%veins_extension_BC(:,5,j))
    call angle_calculation(wings(i)%x,wings(i)%y,wings(i)%z, &
                           wings(i)%veins_bending_BC(:,:,j), &
                           wings(i)%veins_bending_BC(:,7,j), &
                           wings(i)%veins_bending_BC(:,8,j))
end do

!Construct internal force vector
do j=1:nMembranes
     call internal_extension_force(wings(i)%F_int,
                                   wings(i)%x,wings(i)%y,wings(i)%z, &
                                   wings(i)%membranes_extension(:,:,j), &
                                   wings(i)%ke_m(:,j))
enddo

call internal_extension_force(wings(i)%F_int, &
                              wings(i)%x,wings(i)%y,wings(i)%z, &
                              wings(i)%membrane_edge, &
                              wings(i)%ke_me)

do j=1:nVeins
    call internal_extension_force(wings(i)%F_int, &
                                  wings(i)%x,wings(i)%y,wings(i)%z, &
                                  wings(i)%veins_extension(:,:,j), &
                                  wings(i)%ke_v(:,j))
    call internal_bending_force(wings(i)%F_int, &
                                wings(i)%x,wings(i)%y,wings(i)%z, &
                                wings(i)%veins_bending(:,:,j), &
                                wings(i)%kby(:,j),wings(i)%kbz(:,j))
enddo


do j=1:nVeins_BC
    call internal_extension_force_BC(wings(i)%F_int, &
                                     wings(i)%x,wings(i)%y,wings(i)%z, &
                                     wings(i)%veins_extension_BC(:,:,j), &
                                     wings(i)%ke_vBC(:,j))
    call internal_bending_force_BC(wings(i)%F_int, &
                                   wings(i)%x,wings(i)%y,wings(i)%z, &
                                   wings(i)%veins_bending_BC(:,:,j), &
                                   wings(i)%kby_BC(:,j),wings(i)%kbz_BC(:,j))

end do

end subroutine

subroutine internal_extension_force(F_int,x,y,z,extension_springs,ke)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: ke
  real(kind=pr), intent(in)     :: extension_springs
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=1:nint(maxval(extension_springs(:,1)))

      call calculate_extension_spring_force( x(nint(extension_springs(ind,2))), &
                                             x(nint(extension_springs(ind,3))), &
                                             y(nint(extension_springs(ind,2))), &
                                             y(nint(extension_springs(ind,3))), &
                                             z(nint(extension_springs(ind,2))), &
                                             z(nint(extension_springs(ind,3))), &
                                             extension_springs(ind,4), &
                                             extension_springs(ind,5), &
                                             ke(ind), Fx, Fy, Fz)

             F_int(nint(extension_springs(ind,2))) = F_int(nint(extension_springs(ind,2))) + Fx

        F_int(nint(extension_springs(ind,2)) + np) = F_int(nint(extension_springs(ind,2)) + np) + Fy

      F_int(nint(extension_springs(ind,2)) + 2*np) = F_int(nint(extension_springs(ind,2)) + 2*np) + Fz

             F_int(nint(extension_springs(ind,3))) = F_int(nint(extension_springs(ind,3))) - Fx

        F_int(nint(extension_springs(ind,3)) + np) = F_int(nint(extension_springs(ind,3)) + np) - Fy

      F_int(nint(extension_springs(ind,3)) + 2*np) = F_int(nint(extension_springs(ind,3)) + 2*np) - Fz

  enddo

end subroutine

subroutine internal_extension_force_BC(F_int,x,y,z,extension_springs,ke)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: ke_BC
  real(kind=pr), intent(in)     :: extension_springs
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=-1:nint(maxval(extension_springs(:,1)))

      call calculate_extension_spring_force( x(nint(extension_springs(ind,2))), &
                                             x(nint(extension_springs(ind,3))), &
                                             y(nint(extension_springs(ind,2))), &
                                             y(nint(extension_springs(ind,3))), &
                                             z(nint(extension_springs(ind,2))), &
                                             z(nint(extension_springs(ind,3))), &
                                             extension_springs(ind,4), &
                                             extension_springs(ind,5), &
                                             ke_BC(ind), Fx, Fy, Fz)

             F_int(nint(extension_springs(ind,2))) = F_int(nint(extension_springs(ind,2))) + Fx

        F_int(nint(extension_springs(ind,2)) + np) = F_int(nint(extension_springs(ind,2)) + np) + Fy

      F_int(nint(extension_springs(ind,2)) + 2*np) = F_int(nint(extension_springs(ind,2)) + 2*np) + Fz

             F_int(nint(extension_springs(ind,3))) = F_int(nint(extension_springs(ind,3))) - Fx

        F_int(nint(extension_springs(ind,3)) + np) = F_int(nint(extension_springs(ind,3)) + np) - Fy

      F_int(nint(extension_springs(ind,3)) + 2*np) = F_int(nint(extension_springs(ind,3)) + 2*np) - Fz

  enddo

end subroutine

subroutine internal_bending_force(F_int,x,y,z,bending_springs,kby,kbz)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: kby, kbz
  real(kind=pr), intent(in)     :: bending_springs
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=1:nint(maxval(bending_springs(:,1)))

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

           F_int(nint(bending_springs(ind,2))) = F_int(nint(bending_springs(ind,2))) + Fx

      F_int(nint(bending_springs(ind,2)) + np) = F_int(nint(bending_springs(ind,2)) + np) + Fy

    F_int(nint(bending_springs(ind,2)) + 2*np) = F_int(nint(bending_springs(ind,2)) + 2*np) + Fz

           F_int(nint(bending_springs(ind,3))) = F_int(nint(bending_springs(ind,3))) - Fx

      F_int(nint(bending_springs(ind,3)) + np) = F_int(nint(bending_springs(ind,3)) + np) - Fy

    F_int(nint(bending_springs(ind,3)) + 2*np) = F_int(nint(bending_springs(ind,3)) + 2*np) - Fz

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

           F_int(nint(bending_springs(ind,3))) = F_int(nint(bending_springs(ind,3))) - Fx

      F_int(nint(bending_springs(ind,3)) + np) = F_int(nint(bending_springs(ind,3)) + np) - Fy

    F_int(nint(bending_springs(ind,3)) + 2*np) = F_int(nint(bending_springs(ind,3)) + 2*np) - Fz

           F_int(nint(bending_springs(ind,4))) = F_int(nint(bending_springs(ind,4))) + Fx

      F_int(nint(bending_springs(ind,4)) + np) = F_int(nint(bending_springs(ind,4)) + np) + Fy

    F_int(nint(bending_springs(ind,4)) + 2*np) = F_int(nint(bending_springs(ind,4)) + 2*np) + Fz

  enddo

end subroutine

subroutine internal_bending_force_BC(F_int,x,y,z,bending_springs_BC,kby_BC,kbz_BC)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: kby_BC, kbz_BC
  real(kind=pr), intent(in)     :: bending_springs_BC
  real(kind=pr), intent(inout)  :: F_int
  real(kind=pr) :: Fx, Fy, Fz
  integer :: ind

  do ind=-1:maxval(bending_springs_BC(:,1))

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

           F_int(nint(bending_springs_BC(ind,2))) = F_int(nint(bending_springs_BC(ind,2))) + Fx

      F_int(nint(bending_springs_BC(ind,2)) + np) = F_int(nint(bending_springs_BC(ind,2)) + np) + Fy

    F_int(nint(bending_springs_BC(ind,2)) + 2*np) = F_int(nint(bending_springs_BC(ind,2)) + 2*np) + Fz

           F_int(nint(bending_springs_BC(ind,3))) = F_int(nint(bending_springs_BC(ind,3))) - Fx

      F_int(nint(bending_springs_BC(ind,3)) + np) = F_int(nint(bending_springs_BC(ind,3)) + np) - Fy

    F_int(nint(bending_springs_BC(ind,3)) + 2*np) = F_int(nint(bending_springs_BC(ind,3)) + 2*np) - Fz

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

           F_int(nint(bending_springs_BC(ind,3))) = F_int(nint(bending_springs_BC(ind,3))) - Fx

      F_int(nint(bending_springs_BC(ind,3)) + np) = F_int(nint(bending_springs_BC(ind,3)) + np) - Fy

    F_int(nint(bending_springs_BC(ind,3)) + 2*np) = F_int(nint(bending_springs_BC(ind,3)) + 2*np) - Fz

           F_int(nint(bending_springs_BC(ind,4))) = F_int(nint(bending_springs_BC(ind,4))) + Fx

      F_int(nint(bending_springs_BC(ind,4)) + np) = F_int(nint(bending_springs_BC(ind,4)) + np) + Fy

    F_int(nint(bending_springs_BC(ind,4)) + 2*np) = F_int(nint(bending_springs_BC(ind,4)) + 2*np) + Fz

  enddo

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
