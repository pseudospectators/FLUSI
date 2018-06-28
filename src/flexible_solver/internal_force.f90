!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct internal force vector consists of all elastic forces coming from
! extension and bending springs exerted to mass points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine internal_forces_construction(x,y,z,ke_m,ke_me,ke_v,kby,kbz,ke_vBC,kby_BC,kbz_BC, &
F_int,membranes_extension,membrane_edge,veins_bending,veins_extension, &
veins_bending_BC,veins_extension_BC)

real(kind=pr), intent(in)  :: x, y, z
real(kind=pr), intent(in)  :: ke_me,ke_v,kby,kbz,ke_vBC,kby_BC,kbz_BC
real(kind=pr), intent(inout)  :: membranes_extension
real(kind=pr), intent(inout)  :: membrane_edge
real(kind=pr), intent(inout)  :: veins_bending
real(kind=pr), intent(inout)  :: veins_extension
real(kind=pr), intent(inout)  :: veins_bending_BC
real(kind=pr), intent(inout)  :: veins_extension_BC
real(kind=pr), allocatable, intent(out)  :: F_int
integer :: j

allocate(F_int(1:3*np))

!Calculate current lenghts and angles
do j=1:nMembranes
    call length_calculation(x,y,z,membranes_extension(:,:,j),membranes_extension(:,5,j))
enddo

call length_calculation(x,y,z,membrane_edge,membrane_edge(:,5))

do j=1:nVeins
    call length_calculation(x,y,z,veins_extension(:,:,j),veins_extension(:,5,j))
    call angle_calculation( x,y,z,veins_bending(:,:,j),veins_bending(:,7,j), veins_bending(:,8,j))
enddo

do j=1:nVeins_BC
    call length_calculation(x,y,z,veins_extension_BC(:,:,j),veins_extension_BC(:,5,j))
    call angle_calculation( x,y,z,veins_bending_BC(:,:,j),&
    veins_bending_BC(:,7,j), veins_bending_BC(:,8,j))
end

!Construct internal force vector
do j=1:nMembranes
     call internal_extension_force(F_int,x,y,z,np,membranes_extension(:,:,j),ke_m(:,j))
enddo

call internal_extension_force(F_int,x,y,z,np,membrane_edge,ke_me)

do j=1:nVeins
    call internal_extension_force(F_int,x,y,z,np,veins_extension(:,:,j),ke_v(:,j))
    call internal_bending_force(F_int,x,y,z,np,veins_bending(:,:,j),kby(:,j),kbz(:,j))
enddo


do j=1:nVeins_BC
    call internal_extension_force_BC(F_int,x,y,z,np,veins_extension_BC(:,:,j),ke_vBC(:,j))
    call internal_bending_force_BC(F_int,x,y,z,np,veins_bending_BC(:,:,j),kby_BC(:,j),kbz_BC(:,j))
    ![alpha1(j),alpha2(j),beta1(j),beta2(j)] = angle_BC_calculation(...
    !                             x1(j),x2(j),x(veins_bending_BC(1,2,j)),x(veins_bending_BC(1,3,j)),...
  !                               y1(j),y2(j),y(veins_bending_BC(1,2,j)),y(veins_bending_BC(1,3,j)),...
  !                               z1(j),z2(j),z(veins_bending_BC(1,2,j)),z(veins_bending_BC(1,3,j)))
  !  l1(j) = sqrt((x(veins_bending_BC(1,2,j))-x2(j))^2+(y(veins_bending_BC(1,2,j))-y2(j))^2+(z(veins_bending_BC(1,2,j))-z2(j))^2)

    !F_int = Fint_BC(F_int,x,y,z,np,veins_bending_BC(:,:,j),x2(j),y2(j),z2(j),...
    !                    alpha2_0(j),beta2_0(j),alpha1(j),beta1(j),...
    !                    alpha2(j),beta2(j),kby1(j),kbz1(j),kby2(j),...
    !                    kbz2(j),l1_0(j),l1(j),ke1(j))
enddo


end subroutine

subroutine internal_extension_force(F_int,x,y,z,extension_springs,ke)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: ke
  real(kind=pr), intent(in)     :: extension_springs
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=1:maxval(extension_springs(:,1))

      call calculate_extension_spring_force( x(extension_springs(ind,2)), &
                                             x(extension_springs(ind,3)), &
                                             y(extension_springs(ind,2)), &
                                             y(extension_springs(ind,3)), &
                                             z(extension_springs(ind,2)), &
                                             z(extension_springs(ind,3)), &
                                             extension_springs(ind,4), &
                                             extension_springs(ind,5), &
                                             ke(ind), Fx, Fy, Fz)

             F_int(extension_springs(ind,2)) = F_int(extension_springs(ind,2)) + Fx

        F_int(extension_springs(ind,2) + np) = F_int(extension_springs(ind,2) + np) + Fy

      F_int(extension_springs(ind,2) + 2*np) = F_int(extension_springs(ind,2) + 2*np) + Fz

             F_int(extension_springs(ind,3)) = F_int(extension_springs(ind,3)) - Fx

        F_int(extension_springs(ind,3) + np) = F_int(extension_springs(ind,3) + np) - Fy

      F_int(extension_springs(ind,3) + 2*np) = F_int(extension_springs(ind,3) + 2*np) - Fz

  enddo

end subroutine

subroutine internal_bending_force(F_int,x,y,z,bending_springs,kby,kbz)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: kby, kbz
  real(kind=pr), intent(in)     :: bending_springs
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=1:maxval(bending_springs(:,1))

    !Construct force vector for the left part of the bending spring
    call calculate_bending_spring_force(x(bending_springs(ind,2)), &
                                        x(bending_springs(ind,3)), &
                                        y(bending_springs(ind,2)), &
                                        y(bending_springs(ind,3)), &
                                        z(bending_springs(ind,2)), &
                                        z(bending_springs(ind,3)), &
                                        bending_springs(ind,5), &
                                        bending_springs(ind,6), &
                                        bending_springs(ind,7), &
                                        bending_springs(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           F_int(bending_springs(ind,2)) = F_int(bending_springs(ind,2)) + Fx

      F_int(bending_springs(ind,2) + np) = F_int(bending_springs(ind,2) + np) + Fy

    F_int(bending_springs(ind,2) + 2*np) = F_int(bending_springs(ind,2) + 2*np) + Fz

           F_int(bending_springs(ind,3)) = F_int(bending_springs(ind,3)) - Fx

      F_int(bending_springs(ind,3) + np) = F_int(bending_springs(ind,3) + np) - Fy

    F_int(bending_springs(ind,3) + 2*np) = F_int(bending_springs(ind,3) + 2*np) - Fz

    !Construct force vector for the right part of the bending spring
    call calculate_bending_spring_force(x(bending_springs(ind,3)), &
                                        x(bending_springs(ind,4)), &
                                        y(bending_springs(ind,3)), &
                                        y(bending_springs(ind,4)), &
                                        z(bending_springs(ind,3)), &
                                        z(bending_springs(ind,4)), &
                                        bending_springs(ind,5), &
                                        bending_springs(ind,6), &
                                        bending_springs(ind,7), &
                                        bending_springs(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           F_int(bending_springs(ind,3)) = F_int(bending_springs(ind,3)) - Fx

      F_int(bending_springs(ind,3) + np) = F_int(bending_springs(ind,3) + np) - Fy

    F_int(bending_springs(ind,3) + 2*np) = F_int(bending_springs(ind,3) + 2*np) - Fz

           F_int(bending_springs(ind,4)) = F_int(bending_springs(ind,4)) + Fx

      F_int(bending_springs(ind,4) + np) = F_int(bending_springs(ind,4) + np) + Fy

    F_int(bending_springs(ind,4) + 2*np) = F_int(bending_springs(ind,4) + 2*np) + Fz

  enddo

end subroutine

subroutine internal_bending_force_BC(F_int,x,y,z,bending_springs_BC,kby_BC,kbz_BC)

  real(kind=pr), intent(in)     :: x, y, z
  real(kind=pr), intent(in)     :: kby_BC, kbz_BC
  real(kind=pr), intent(in)     :: bending_springs_BC
  real(kind=pr), intent(inout)  :: F_int
  integer :: ind

  do ind=1:maxval(bending_springs_BC(:,1))

    !Construct force vector for the left part of the bending spring
    call calculate_bending_spring_force(x(bending_springs_BC(ind,2)), &
                                        x(bending_springs_BC(ind,3)), &
                                        y(bending_springs_BC(ind,2)), &
                                        y(bending_springs_BC(ind,3)), &
                                        z(bending_springs_BC(ind,2)), &
                                        z(bending_springs_BC(ind,3)), &
                                        bending_springs_BC(ind,5), &
                                        bending_springs_BC(ind,6), &
                                        bending_springs_BC(ind,7), &
                                        bending_springs_BC(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           F_int(bending_springs_BC(ind,2)) = F_int(bending_springs_BC(ind,2)) + Fx

      F_int(bending_springs_BC(ind,2) + np) = F_int(bending_springs_BC(ind,2) + np) + Fy

    F_int(bending_springs_BC(ind,2) + 2*np) = F_int(bending_springs_BC(ind,2) + 2*np) + Fz

           F_int(bending_springs_BC(ind,3)) = F_int(bending_springs_BC(ind,3)) - Fx

      F_int(bending_springs_BC(ind,3) + np) = F_int(bending_springs_BC(ind,3) + np) - Fy

    F_int(bending_springs_BC(ind,3) + 2*np) = F_int(bending_springs_BC(ind,3) + 2*np) - Fz

    !Construct force vector for the right part of the bending spring
    call calculate_bending_spring_force(x(bending_springs_BC(ind,3)), &
                                        x(bending_springs_BC(ind,4)), &
                                        y(bending_springs_BC(ind,3)), &
                                        y(bending_springs_BC(ind,4)), &
                                        z(bending_springs_BC(ind,3)), &
                                        z(bending_springs_BC(ind,4)), &
                                        bending_springs_BC(ind,5), &
                                        bending_springs_BC(ind,6), &
                                        bending_springs_BC(ind,7), &
                                        bending_springs_BC(ind,8), &
                                        kby(ind), kbz(ind), Fx, Fy, Fz)

           F_int(bending_springs_BC(ind,3)) = F_int(bending_springs_BC(ind,3)) - Fx

      F_int(bending_springs_BC(ind,3) + np) = F_int(bending_springs_BC(ind,3) + np) - Fy

    F_int(bending_springs_BC(ind,3) + 2*np) = F_int(bending_springs_BC(ind,3) + 2*np) - Fz

           F_int(bending_springs_BC(ind,4)) = F_int(bending_springs_BC(ind,4)) + Fx

      F_int(bending_springs_BC(ind,4) + np) = F_int(bending_springs_BC(ind,4) + np) + Fy

    F_int(bending_springs_BC(ind,4) + 2*np) = F_int(bending_springs_BC(ind,4) + 2*np) + Fz

  enddo

  

end subroutine
