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

subroutine calculate_extension_spring_force_derivative(x1, x2, y1, y2, z1, z2, &
l0, l, ke, Fxx, Fyy, Fzz, Fxy, Fxz, Fyz )

real(kind=pr), intent(in)  :: x1, x2, y1, y2, z1, z2
real(kind=pr), intent(in)  :: l0, l
real(kind=pr), intent(in)  :: ke
real(kind=pr), intent(out) :: Fxx, Fyy, Fzz, Fxy, Fxz, Fyz
real(kind=pr) :: dl, xl, yl, zl

dl = (l-l0)/l
xl = (x2-x1)/l
yl = (y2-y1)/l
zl = (z2-z1)/l

Fxx = ke*(dl + xl^2 - dl*xl^2)

Fyy = ke*(dl + yl^2 - dl*yl^2)

Fzz = ke*(dl + zl^2 - dl*zl^2)

Fxy = ke*(xl*yl - dl*xl*yl)

Fxz = ke*(xl*zl - dl*xl*zl)

Fyz = ke*(yl*zl - dl*yl*zl)

end subroutine

subroutine calculate_bending_spring_force( x1, x2, y1, y2, z1, z2, &
  theta0, phi0, theta, phi, kby, kbz, Fx, Fy, Fz )

real(kind=pr), intent(in)  :: x1, x2, y1, y2, z1, z2
real(kind=pr), intent(in)  :: theta0, theta, phi0, phi
real(kind=pr), intent(in)  :: kby, kbz
real(kind=pr), intent(out) :: Fx, Fy, Fz
real(kind=pr) :: idxy, idxz

if (((x2-x1)^2+(y2-y1)^2)<1e-12) then
    idxy=0
else
    idxy = 1/((x2-x1)^2+(y2-y1)^2)
endif

if (((x2-x1)^2+(z2-z1)^2)<1e-12) then
    idxz=0
else
    idxz = 1/((x2-x1)^2+(z2-z1)^2)
endif

Fx = kby*(theta - theta0)*(y2-y1)*idxy + kbz*(phi - phi0)*(z2-z1)*idxz

Fy = - kby*(theta - theta0)*(x2-x1)*idxy

Fz = - kbz*(phi - phi0)*(x2-x1)*idxz

end subroutine

subroutine calculate_bending_spring_force_derivative(x1, x2, x3, y1, y2, y3,&
z1, z2, z3, theta0, phi0, theta, phi, kby, kbz, position, &
Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1,Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

!
!   X1        X2           X3
!   O---------O------------O
real(kind=pr), intent(in)  :: x1, x2, x3, y1, y2, y3, z1, z2, z3
real(kind=pr), intent(in)  :: theta0, theta, phi0, phi
real(kind=pr), intent(in)  :: kby, kbz
real(kind=pr), intent(out) :: Fx, Fy, Fz
real(kind=pr) :: xy1, xy2, xz1, xz2, xxy1, xxz1, yxy1, zxz1, xxy2, xxz2, yxy2, zxz2
real(kind=pr) :: Ay, By, Az, Bz, kbyT, kbzP

 xy1 = sqrt((x2-x1)^2+(y2-y1)^2)
 xy2 = sqrt((x3-x2)^2+(y3-y2)^2)
 xz1 = sqrt((x2-x1)^2+(z2-z1)^2)
 xz2 = sqrt((x3-x2)^2+(z3-z2)^2)

xxy1 = (x2-x1)/(xy1)^2
xxz1 = (x2-x1)/(xz1)^2
yxy1 = (y2-y1)/(xy1)^2
zxz1 = (z2-z1)/(xz1)^2

xxy2 = (x3-x2)/(xy2)^2
xxz2 = (x3-x2)/(xz2)^2
yxy2 = (y3-y2)/(xy2)^2
zxz2 = (z3-z2)/(xz2)^2

   Ay = ((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))/ &
       (((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))^2+((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))^2)
   By = ((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))/ &
       (((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))^2+((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))^2)
   Az = ((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))/ &
       (((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))^2+((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))^2)
   Bz = ((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))/ &
       (((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))^2+((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))^2)

kbyT = kby*(theta - theta0)
kbzP = kbz*(phi - phi0)

select case (position)
  case 2 !left point
      !Derivative with respect to the left point X1
        Fxx1 = kby*yxy1*((y3-y2)*Ay + (x3-x2)*By) - 2*kbyT*yxy1*xxy1 +...
               kbz*zxz1*((z3-z2)*Az + (x3-x2)*Bz) - 2*kbzP*zxz1*xxz1

        Fyy1 = kby*xxy1*((x3-x2)*Ay - (y3-y2)*By) + 2*kbyT*yxy1*xxy1

        Fzz1 = kbz*xxz1*((x3-x2)*Az - (z3-z2)*Bz) + 2*kbzP*zxz1*xxz1

        Fxy1 = - kby*yxy1*((x3-x2)*Ay - (y3-y2)*By) + kbyT*(1/(xy1)^2 - 2*(yxy1)^2)

        Fxz1 = - kbz*zxz1*((x3-x2)*Az - (z3-z2)*Bz) + kbzP*(1/(xz1)^2 - 2*(zxz1)^2)

        Fyx1 = - kby*xxy1*((x3-x2)*By + (y3-y2)*Ay) - kbyT*(1/(xy1)^2 - 2*(xxy1)^2)

        Fzx1 = - kbz*xxz1*((x3-x2)*Bz + (z3-z2)*Az) - kbzP*(1/(xz1)^2 - 2*(xxz1)^2)

        !Derivative with respect to the middle point X2
        Fxx2 = - kby*yxy1*((y3-y1)*Ay + (x3-2*x2+x1)*By) + 2*kbyT*yxy1*xxy1 + &
               - kbz*zxz1*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + 2*kbzP*zxz1*xxz1

        Fyy2 = - kby*xxy1*((x3-x1)*Ay - (y3-2*y2+y1)*By) - 2*kbyT*yxy1*xxy1

        Fzz2 = - kbz*xxz1*((x3-x1)*Az - (z3-2*z2+z1)*Bz) - 2*kbzP*zxz1*xxz1

        Fxy2 = kby*yxy1*((x3-x1)*Ay - (y3-2*y2+y1)*By) - kbyT*(1/(xy1)^2 - 2*(yxy1)^2)

        Fxz2 = kbz*zxz1*((x3-x1)*Az - (z3-2*z2+z1)*Bz) - kbzP*(1/(xz1)^2 - 2*(zxz1)^2)

        Fyx2 = kby*xxy1*((y3-y1)*Ay + (x3-2*x2+x1)*By) + kbyT*(1/(xy1)^2 - 2*(xxy1)^2)

        Fzx2 = kbz*xxz1*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + kbzP*(1/(xz1)^2 - 2*(xxz1)^2)

        !Derivative with respect to the third point X3
        Fxx3 = kby*yxy1*((y2-y1)*Ay - (x2-x1)*By) + kbz*zxz1*((z2-z1)*Az - (x2-x1)*Bz)

        Fyy3 = kby*xxy1*((x2-x1)*Ay + (y2-y1)*By)

        Fzz3 = kbz*xxz1*((x2-x1)*Az + (z2-z1)*Bz)

        Fxy3 = - kby*yxy1*((x2-x1)*Ay + (y2-y1)*By)

        Fxz3 = - kbz*zxz1*((x2-x1)*Az + (z2-z1)*Bz)

        Fyx3 = - kby*xxy1*((y2-y1)*Ay - (x2-x1)*By)

        Fzx3 = - kbz*xxz1*((z2-z1)*Az - (x2-x1)*Bz)

    case 3 %middle point
        %Derivative with respect to the left point X1
        %(atan2((x2-x1)*(y3-y2) - (x3-x2)*(y2-y1),(x2-x1)*(x3-x2) + (y3-y2)*(y2-y1))- a)*((y2-y1)/((x2-x1)^2+(y2-y1)^2)+(y3-y2)/((x3-x2)^2+(y3-y2)^2))
        Fxx1 = - kby*(yxy1+yxy2)*((y3-y2)*Ay + (x3-x2)*By) + 2*kbyT*yxy1*xxy1 +...
               - kbz*(zxz1+zxz2)*((z3-z2)*Az + (x3-x2)*Bz) + 2*kbzP*zxz1*xxz1

        Fyy1 = - kby*(xxy1+xxy2)*((x3-x2)*Ay - (y3-y2)*By) - 2*kbyT*yxy1*xxy1

        Fzz1 = - kbz*(xxz1+xxz2)*((x3-x2)*Az - (z3-z2)*Bz) - 2*kbzP*zxz1*xxz1

        Fxy1 = kby*(yxy1+yxy2)*((x3-x2)*Ay - (y3-y2)*By) + kbyT*(2*(yxy1)^2 - 1/(xy1)^2)

        Fxz1 = kbz*(zxz1+zxz2)*((x3-x2)*Az - (z3-z2)*Bz) + kbzP*(2*(zxz1)^2 - 1/(xz1)^2)

        Fyx1 = kby*(xxy2+xxy1)*((y3-y2)*Ay + (x3-x2)*By) - kbyT*(2*(xxy1)^2 - 1/(xy1)^2)

        Fzx1 = kbz*(xxz1+xxz2)*((z3-z2)*Az + (x3-x2)*Bz) - kbzP*(2*(xxz1)^2 - 1/(xz1)^2)

        %Derivative with respect to the middle point X2
        Fxx2 = kby*(yxy1+yxy2)*((y3-y1)*Ay + (x3-2*x2+x1)*By) + 2*kbyT*(yxy2*xxy2 - yxy1*xxy1) +...
               kbz*(zxz1+zxz2)*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + 2*kbzP*(zxz2*xxz2 - zxz1*xxz1)

        Fyy2 = kby*(xxy2+xxy1)*((x3-x1)*Ay - (y3-2*y2+y1)*By) + ...
               - 2*kbyT*(yxy2*xxy2 - yxy1*xxy1)

        Fzz2 = kbz*(xxz2+xxz1)*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + ...
               - 2*kbzP*(zxz2*xxz2 - zxz1*xxz1)

        Fxy2 = - kby*(yxy2+yxy1)*((x3-x1)*Ay - (y3-2*y2+y1)*By) + ...
               - kbyT*(1/(xy2)^2 - 1/(xy1)^2 - 2*(yxy2)^2 + 2*(yxy1)^2)

        Fxz2 = - kbz*(zxz2+zxz1)*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + ...
               - kbzP*(1/(xz2)^2 - 1/(xz1)^2 - 2*(zxz2)^2 + 2*(zxz1)^2)

        Fyx2 = - kby*(xxy1+xxy2)*((y3-y1)*Ay + (x3-2*x2+x1)*By) + ...
               kbyT*(1/(xy2)^2 - 1/(xy1)^2 + 2*(xxy1)^2 - 2*(xxy2)^2)

        Fzx2 = - kbz*(xxz1+xxz2)*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + ...
               kbzP*(1/(xz2)^2 - 1/(xz1)^2 + 2*(xxz1)^2 - 2*(xxz2)^2)

        %Derivative with respect to the third point X3
        Fxx3 = - kby*(yxy1+yxy2)*((y2-y1)*Ay - (x2-x1)*By) - 2*kbyT*yxy2*xxy2 +...
               - kbz*(zxz1+zxz2)*((z2-z1)*Az - (x2-x1)*Bz) - 2*kbzP*zxz2*xxz2

        Fyy3 = - kby*(xxy1+xxy2)*((x2-x1)*Ay + (y2-y1)*By) + 2*kbyT*yxy2*xxy2

        Fzz3 = - kbz*(xxz1+xxz2)*((x2-x1)*Az + (z2-z1)*Bz) + 2*kbzP*zxz2*xxz2

        Fxy3 = kby*(yxy1+yxy2)*((x2-x1)*Ay + (y2-y1)*By) - kbyT*(2*(yxy2)^2 - 1/(xy2)^2)

        Fxz3 = kbz*(zxz1+zxz2)*((x2-x1)*Az + (z2-z1)*Bz) - kbzP*(2*(zxz2)^2 - 1/(xz2)^2)

        Fyx3 = kby*(xxy2+xxy1)*((y2-y1)*Ay - (x2-x1)*By) + kbyT*(2*(xxy2)^2 - 1/(xy2)^2)

        Fzx3 = kbz*(xxz1+xxz2)*((z2-z1)*Az - (x2-x1)*Bz) + kbzP*(2*(xxz2)^2 - 1/(xz2)^2)

    case 4 !right point
       !Derivative with respect to the left point X1
        Fxx1 = kby*yxy2*((y3-y2)*Ay + (x3-x2)*By) +...
               kbz*zxz2*((z3-z2)*Az + (x3-x2)*Bz)

        Fyy1 = kby*xxy2*((x3-x2)*Ay - (y3-y2)*By)

        Fzz1 = kbz*xxz2*((x3-x2)*Az - (z3-z2)*Bz)

        Fxy1 = - kby*yxy2*((x3-x2)*Ay - (y3-y2)*By)

        Fxz1 = - kbz*zxz2*((x3-x2)*Az - (z3-z2)*Bz)

        Fyx1 = - kby*xxy2*((y3-y2)*Ay + (x3-x2)*By)

        Fzx1 = - kbz*xxz2*((z3-z2)*Az + (x3-x2)*Bz)

        !Derivative with respect to the middle point X2
        Fxx2 = - kby*yxy2*((y3-y1)*Ay + (x3-2*x2+x1)*By) - 2*kbyT*yxy2*xxy2 +...
               - kbz*zxz2*((z3-z1)*Az + (x3-2*x2+x1)*Bz) - 2*kbzP*zxz2*xxz2

        Fyy2 = - kby*xxy2*((x3-x1)*Ay - (y3-2*y2+y1)*By) + 2*kbyT*yxy2*xxy2

        Fzz2 = - kbz*xxz2*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + 2*kbzP*zxz2*xxz2

        Fxy2 = kby*yxy2*((x3-x1)*Ay - (y3-2*y2+y1)*By) + kbyT*(1/(xy2)^2 - 2*(yxy2)^2)

        Fxz2 = kbz*zxz2*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + kbzP*(1/(xz2)^2 - 2*(zxz2)^2)

        Fyx2 = kby*xxy2*((y3-y1)*Ay + (x3-2*x2+x1)*By) - kbyT*(1/(xy2)^2 - 2*(xxy2)^2)

        Fzx2 = kbz*xxz2*((z3-z1)*Az + (x3-2*x2+x1)*Bz) - kbzP*(1/(xz2)^2 - 2*(xxz2)^2)


        !Derivative with respect to the third point X3
        Fxx3 = kby*yxy2*((y2-y1)*Ay - (x2-x1)*By) + 2*kbyT*yxy2*xxy2 +...
               kbz*zxz2*((z2-z1)*Az - (x2-x1)*Bz) + 2*kbzP*zxz2*xxz2

        Fyy3 = kby*xxy2*((x2-x1)*Ay + (y2-y1)*By) - 2*kbyT*yxy2*xxy2

        Fzz3 = kbz*xxz2*((x2-x1)*Az + (z2-z1)*Bz) - 2*kbzP*zxz2*xxz2

        Fxy3 = - kby*yxy2*((x2-x1)*Ay + (y2-y1)*By) - kbyT*(1/(xy2)^2 - 2*(yxy2)^2)

        Fxz3 = - kbz*zxz2*((x2-x1)*Az + (z2-z1)*Bz) - kbzP*(1/(xz2)^2 - 2*(zxz2)^2)

        Fyx3 = - kby*xxy2*((y2-y1)*Ay - (x2-x1)*By) + kbyT*(1/(xy2)^2 - 2*(xxy2)^2)

        Fzx3 = - kbz*xxz2*((z2-z1)*Az - (x2-x1)*Bz) + kbzP*(1/(xz2)^2 - 2*(xxz2)^2)

end

end subroutine

subroutine length_calculation(x,y,z,springs,springs_length)

real(kind=pr), intent(in) :: x,y,z
real(kind=pr), intent(in) :: springs
real(kind=pr), allocatable, intent(out) :: springs_length
integer :: ind

allocate(springs_length(1:maxval(springs(:,1))))

do ind=1:maxval(springs(:,1))
    springs_length(ind,1) = sqrt((x(springs(ind,2))-x(springs(ind,3)))^2 + &
                                 (y(springs(ind,2))-y(springs(ind,3)))^2 + &
                                 (z(springs(ind,2))-z(springs(ind,3)))^2)
enddo

end subroutine

subroutine angle_calculation(x,y,z,vein,theta,phi)
!Calculate relative angles between two consecutive segments belonging to a vein
!   x(vein(ind,2))          x(vein(ind,3))                x(vein(ind,4))
!        O-----------------------O----------------------------O
!                  X1                          X2
!

real(kind=pr), intent(in) :: x,y,z
real(kind=pr), intent(in) :: vein
real(kind=pr), allocatable, intent(out) :: theta, phi
integer :: ind

allocate(theta(1:maxval(vein(:,1))))
allocate(phi(1:maxval(vein(:,1))))

do ind=1:maxval(vein(:,1))

    x1 = x(vein(ind,3)) - x(vein(ind,2))
    y1 = y(vein(ind,3)) - y(vein(ind,2))
    z1 = z(vein(ind,3)) - z(vein(ind,2))
    x2 = x(vein(ind,4)) - x(vein(ind,3))
    y2 = y(vein(ind,4)) - y(vein(ind,3))
    z2 = z(vein(ind,4)) - z(vein(ind,3))

    theta(ind) = atan2(x1*y2-x2*y1,x1*x2+y1*y2)

      phi(ind) = atan2(x1*z2-x2*z1,x1*x2+z1*z2)
enddo

end subroutine

subroutine convert_flexural_rigidity_into_spring_stiffnes(EIy, EIz, kby0, kbz0, L, ns)

  implicit none
  real(kind=pr), intent(in) :: EIy, EIz, L, ns
  real(kind=pr), intent(out) :: kby0, kbz0

  kby0 = EIy*(ns-1)/L
  kbz0 = EIz*(ns-1)/L

end subroutine
