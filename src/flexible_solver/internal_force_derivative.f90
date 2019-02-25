!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Construct Jacobian matrix of internal force vector for Newton-Raphson method
! since we use implicit scheme for time stepping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine internal_forces_derivatives_construction(Wings)

type(flexible_wing), intent(inout)  :: Wings
integer :: np,j, ij

! Get the number of mass points for the sake of simplicity in coding
np = wings%np

! Initialize
wings%FJ = 0.d0

  !Construct internal force derivative matrix aka Jacobian matrix
  do j=1,nMembranes
       call internal_extension_force_derivative(wings%FJ, &
                                     wings%u_new(1:np), &
                                     wings%u_new(np+1:2*np), &
                                     wings%u_new(2*np+1:3*np), &
                                     wings%membranes_extension(:,:,j), &
                                     wings%np, &
                                     wings%ke_m(:,j))
  enddo

  do j=1,nMembrane_edges
    call internal_extension_force_derivative(wings%FJ, &
                                  wings%u_new(1:np), &
                                  wings%u_new(np+1:2*np), &
                                  wings%u_new(2*np+1:3*np), &
                                  wings%membrane_edge(:,:,j), &
                                  wings%np, &
                                  wings%ke_me(:,j))
  enddo

  do j=1,nVeins


      call internal_extension_force_derivative(wings%FJ, &
                                    wings%u_new(1:np), &
                                    wings%u_new(np+1:2*np), &
                                    wings%u_new(2*np+1:3*np), &
                                    wings%veins_extension(:,:,j), &
                                    wings%np, &
                                    wings%ke_v(:,j))


      call internal_bending_force_derivative(wings%FJ, &
                                  wings%u_new(1:np), &
                                  wings%u_new(np+1:2*np), &
                                  wings%u_new(2*np+1:3*np), &
                                  wings%veins_bending(:,:,j), &
                                  wings%np, &
                                  wings%kby(:,j),wings%kbz(:,j))
  enddo

  do j=1,nVeins_BC
      call internal_extension_force_BC_derivative(wings%FJ, &
                                       wings%u_new(1:np), &
                                       wings%u_new(np+1:2*np), &
                                       wings%u_new(2*np+1:3*np), &
                                       wings%x_BC(0,j), wings%y_BC(0,j), wings%z_BC(0,j), &
                                       wings%veins_extension_BC(:,:,j), &
                                       wings%np, &
                                       wings%ke_vBC(:,j))


      call internal_bending_force_BC_derivative(wings%FJ, &
                                     wings%u_new(1:np), &
                                     wings%u_new(np+1:2*np), &
                                     wings%u_new(2*np+1:3*np), &
                                     wings%x_BC(-1:0,j), wings%y_BC(-1:0,j), wings%z_BC(-1:0,j), &
                                     wings%veins_bending_BC(:,:,j), &
                                     wings%np, &
                                     wings%kby_BC(:,j),wings%kbz_BC(:,j))

  end do

  call internal_bending_force_derivative(wings%FJ, &
                              wings%u_new(1:np), &
                              wings%u_new(np+1:2*np), &
                              wings%u_new(2*np+1:3*np), &
                              wings%vein_connectors, &
                              wings%np, &
                              wings%kby_c,wings%kbz_c)

end subroutine

subroutine internal_extension_force_derivative(FJa,x,y,z,extension_springs,np,ke)
! Calculate Jacobian matrix elements of forces caused by extension springs

real(kind=pr), intent(in)     :: x(:),y(:),z(:)
real(kind=pr), intent(in)     :: ke(:)
real(kind=pr), intent(in)     :: extension_springs(:,:)
integer, intent(in)           :: np
real(kind=pr), intent(inout)  :: FJa(:,:)
real(kind=pr) :: Fxx, Fyy, Fzz, Fxy, Fxz, Fyz
integer :: ind

do ind=1,nint(maxval(extension_springs(:,1)))

    call calculate_extension_spring_force_derivative(x(nint(extension_springs(ind,2))), x(nint(extension_springs(ind,3))), &
                                                     y(nint(extension_springs(ind,2))), y(nint(extension_springs(ind,3))), &
                                                     z(nint(extension_springs(ind,2))), z(nint(extension_springs(ind,3))), &
                                                     extension_springs(ind,4), extension_springs(ind,5), ke(ind), &
                                                     Fxx, Fyy, Fzz, Fxy, Fxz, Fyz)


    !First point with respect to first point
                FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2))) = &
                FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2))) + Fxx

          FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2)) + np) = &
          FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2)) + np) + Fxy

        FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2)) + 2*np) = &
        FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,2)) + 2*np) + Fxz

          FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2))) = &
          FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2))) + Fxy

    FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2)) + np) = &
    FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2)) + np) + Fyy

  FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2)) + 2*np) = &
  FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,2)) + 2*np) + Fyz

        FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2))) = &
        FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2))) + Fxz

  FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2)) + np) = &
  FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2)) + np) +  Fyz

FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2)) + 2*np) = &
FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,2)) + 2*np) + Fzz

    ! First point with respect to second point

                FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3))) = &
                FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3))) - Fxx

          FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3)) + np) = &
          FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3)) + np) - Fxy

        FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3)) + 2*np) = &
        FJa(nint(extension_springs(ind,2)),nint(extension_springs(ind,3)) + 2*np) - Fxz

          FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3))) = &
          FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3))) - Fxy

    FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3)) + np) = &
    FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3)) + np) - Fyy

  FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3)) + 2*np) = &
  FJa(nint(extension_springs(ind,2)) + np,nint(extension_springs(ind,3)) + 2*np) - Fyz

        FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3))) = &
        FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3))) - Fxz

  FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3)) + np) = &
  FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3)) + np) - Fyz

FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3)) + 2*np) = &
FJa(nint(extension_springs(ind,2)) + 2*np,nint(extension_springs(ind,3)) + 2*np) - Fzz

    ! Second point with respect to first point

                FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2))) = &
                FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2))) - Fxx

          FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2)) + np) = &
          FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2)) + np) - Fxy

        FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2)) + 2*np) = &
        FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,2)) + 2*np) - Fxz

          FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2))) = &
          FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2))) - Fxy

    FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2)) + np) = &
    FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2)) + np) - Fyy

  FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2)) + 2*np) = &
  FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,2)) + 2*np) - Fyz

        FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2))) = &
        FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2))) - Fxz

  FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2)) + np) = &
  FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2)) + np) - Fyz

FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2)) + 2*np) = &
FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,2)) + 2*np) - Fzz

    ! Second point with respect to second point

                FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3))) = &
                FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3))) + Fxx

          FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3)) + np) = &
          FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3)) + np) + Fxy

        FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3)) + 2*np) = &
        FJa(nint(extension_springs(ind,3)),nint(extension_springs(ind,3)) + 2*np) + Fxz

          FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3))) = &
          FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3))) + Fxy

    FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3)) + np) = &
    FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3)) + np) + Fyy

  FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3)) + 2*np) = &
  FJa(nint(extension_springs(ind,3)) + np,nint(extension_springs(ind,3)) + 2*np) + Fyz

        FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3))) = &
        FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3))) + Fxz

  FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3)) + np) = &
  FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3)) + np) + Fyz

FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3)) + 2*np) = &
FJa(nint(extension_springs(ind,3)) + 2*np,nint(extension_springs(ind,3)) + 2*np) + Fzz

enddo

end subroutine

subroutine internal_extension_force_BC_derivative(FJa,x,y,z,x_BC,y_BC,z_BC,extension_springs_BC,np,ke_BC)

! Calculate Jacobian matrix elements of forces caused by extension springs of veins
! with boundary conditions

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: x_BC,y_BC,z_BC
  real(kind=pr), intent(in)     :: ke_BC(0:)
  real(kind=pr), intent(in)     :: extension_springs_BC(0:,1:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: FJa(:,:)
  real(kind=pr) :: Fxx, Fyy, Fzz, Fxy, Fxz, Fyz
  integer :: ind

  do ind=1,nint(maxval(extension_springs_BC(:,1)))

      call calculate_extension_spring_force_derivative(x(nint(extension_springs_BC(ind,2))), x(nint(extension_springs_BC(ind,3))), &
                                                       y(nint(extension_springs_BC(ind,2))), y(nint(extension_springs_BC(ind,3))), &
                                                       z(nint(extension_springs_BC(ind,2))), z(nint(extension_springs_BC(ind,3))), &
                                                       extension_springs_BC(ind,4), extension_springs_BC(ind,5), ke_BC(ind), &
                                                       Fxx, Fyy, Fzz, Fxy, Fxz, Fyz)


      !First point with respect to first point
                  FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2))) = &
                  FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2))) + Fxx

            FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2)) + np) = &
            FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2)) + np) + Fxy

          FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2)) + 2*np) = &
          FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,2)) + 2*np) + Fxz

            FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2))) = &
            FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2))) + Fxy

      FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2)) + np) = &
      FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2)) + np) + Fyy

    FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2)) + 2*np) = &
    FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,2)) + 2*np) + Fyz

          FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2))) = &
          FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2))) + Fxz

    FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2)) + np) = &
    FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2)) + np) +  Fyz

  FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2)) + 2*np) = &
  FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,2)) + 2*np) + Fzz

      ! First point with respect to second point

                  FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3))) = &
                  FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3))) - Fxx

            FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3)) + np) = &
            FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3)) + np) - Fxy

          FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3)) + 2*np) = &
          FJa(nint(extension_springs_BC(ind,2)),nint(extension_springs_BC(ind,3)) + 2*np) - Fxz

            FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3))) = &
            FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3))) - Fxy

      FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3)) + np) = &
      FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3)) + np) - Fyy

    FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3)) + 2*np) = &
    FJa(nint(extension_springs_BC(ind,2)) + np,nint(extension_springs_BC(ind,3)) + 2*np) - Fyz

          FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3))) = &
          FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3))) - Fxz

    FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3)) + np) = &
    FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3)) + np) - Fyz

  FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3)) + 2*np) = &
  FJa(nint(extension_springs_BC(ind,2)) + 2*np,nint(extension_springs_BC(ind,3)) + 2*np) - Fzz

      ! Second point with respect to first point

                  FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2))) = &
                  FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2))) - Fxx

            FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2)) + np) = &
            FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2)) + np) - Fxy

          FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2)) + 2*np) = &
          FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,2)) + 2*np) - Fxz

            FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2))) = &
            FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2))) - Fxy

      FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2)) + np) = &
      FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2)) + np) - Fyy

    FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2)) + 2*np) = &
    FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,2)) + 2*np) - Fyz

          FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2))) = &
          FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2))) - Fxz

    FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2)) + np) = &
    FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2)) + np) - Fyz

  FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2)) + 2*np) = &
  FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,2)) + 2*np) - Fzz

      ! Second point with respect to second point

                  FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3))) = &
                  FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3))) + Fxx

            FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3)) + np) = &
            FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3)) + np) + Fxy

          FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3)) + 2*np) = &
          FJa(nint(extension_springs_BC(ind,3)),nint(extension_springs_BC(ind,3)) + 2*np) + Fxz

            FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3))) = &
            FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3))) + Fxy

      FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3)) + np) = &
      FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3)) + np) + Fyy

    FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3)) + 2*np) = &
    FJa(nint(extension_springs_BC(ind,3)) + np,nint(extension_springs_BC(ind,3)) + 2*np) + Fyz

          FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3))) = &
          FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3))) + Fxz

    FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3)) + np) = &
    FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3)) + np) + Fyz

  FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3)) + 2*np) = &
  FJa(nint(extension_springs_BC(ind,3)) + 2*np,nint(extension_springs_BC(ind,3)) + 2*np) + Fzz

  enddo

  !---------------------------------------------------------------------------------
  ! Calculate derivative of forces caused by extension spring connecting the vein with the boundary
  !---------------------------------------------------------------------------------
  ! Vein schematic: point -1 and 0 are boundary points and they are clamped or
  ! rotated. Extension spring (numbered 0) connects points 0 and 1
  !
  !  -1       0         1             2       3
  !             ke_BC(0)
  !   X-------X---------O-------------O-------O .......
  !             l0(0)
  !

  ! For extension spring number 0 connecting two points 0 and 1

      ! Calculate derivative of force of the extension spring 0 acting on point 1
        call calculate_extension_spring_force_derivative(x_BC,x(nint(extension_springs_BC(1,2))), &
                                                         y_BC,y(nint(extension_springs_BC(1,2))), &
                                                         z_BC,z(nint(extension_springs_BC(1,2))), &
                                                         extension_springs_BC(0,4),extension_springs_BC(0,5), &
                                                         ke_BC(0),Fxx, Fyy, Fzz, Fxy, Fxz, Fyz)

      ! Adding these terms into Jacobian matrix
        ! Second point with respect to second point
        FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2))) = &
        FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2))) + Fxx

  FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2)) + np) = &
  FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2)) + np) + Fxy

FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2)) + 2*np) = &
FJa(nint(extension_springs_BC(1,2)),nint(extension_springs_BC(1,2)) + 2*np) + Fxz

  FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2))) = &
  FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2))) + Fxy

FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2)) + np) = &
FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2)) + np) + Fyy

FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2)) + 2*np) = &
FJa(nint(extension_springs_BC(1,2)) + np,nint(extension_springs_BC(1,2)) + 2*np) + Fyz

FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2))) = &
FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2))) + Fxz

FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2)) + np) = &
FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2)) + np) +  Fyz

FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2)) + 2*np) = &
FJa(nint(extension_springs_BC(1,2)) + 2*np,nint(extension_springs_BC(1,2)) + 2*np) + Fzz

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

Fxx = ke*(dl + xl**2 - dl*xl**2)

Fyy = ke*(dl + yl**2 - dl*yl**2)

Fzz = ke*(dl + zl**2 - dl*zl**2)

Fxy = ke*(xl*yl - dl*xl*yl)

Fxz = ke*(xl*zl - dl*xl*zl)

Fyz = ke*(yl*zl - dl*yl*zl)

end subroutine calculate_extension_spring_force_derivative


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine internal_bending_force_derivative(FJa,x,y,z,bending_springs,np,kby,kbz)
! Calculate Jacobian matrix elements of forces caused by bending springs

  implicit none

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: kby(:), kbz(:)
  real(kind=pr), intent(in)     :: bending_springs(:,:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: FJa(:,:)
  real(kind=pr) :: Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1
  real(kind=pr) :: Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2
  real(kind=pr) :: Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3
  integer :: ind

  do ind=1,nint(maxval(bending_springs(:,1)))

      !left point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs(ind,2))),x(nint(bending_springs(ind,3))),x(nint(bending_springs(ind,4))),&
           y(nint(bending_springs(ind,2))),y(nint(bending_springs(ind,3))),y(nint(bending_springs(ind,4))),&
           z(nint(bending_springs(ind,2))),z(nint(bending_springs(ind,3))),z(nint(bending_springs(ind,4))),&
           bending_springs(ind,5),bending_springs(ind,6),bending_springs(ind,7),bending_springs(ind,8),&
           kby(ind), kbz(ind), 2, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))) = &
                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))) + Fxx1

          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,2))+np) = &
          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,2))+np) + Fyy1

      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,2))+2*np) = &
      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))+np) = &
              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))+np) + Fxy1

            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))+2*np) = &
            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,2))) = &
              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,2))) + Fyx1

            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,2))) = &
            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,2))) + Fzx1

                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))) = &
                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))) + Fxx2

          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,3))+np) = &
          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,3))+np) + Fyy2

      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,3))+2*np) = &
      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))+np) = &
              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))+np) + Fxy2

            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))+2*np) = &
            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,3))) = &
              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,3))) + Fyx2

            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,3))) = &
            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,3))) + Fzx2

                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))) = &
                  FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))) + Fxx3

          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,4))+np) = &
          FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,4))+np) + Fyy3

      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,4))+2*np) = &
      FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))+np) = &
              FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))+np) + Fxy3

            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))+2*np) = &
            FJa(nint(bending_springs(ind,2)),nint(bending_springs(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,4))) = &
              FJa(nint(bending_springs(ind,2))+np,nint(bending_springs(ind,4))) + Fyx3

            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,4))) = &
            FJa(nint(bending_springs(ind,2))+2*np,nint(bending_springs(ind,4))) + Fzx3

      !middle point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs(ind,2))),x(nint(bending_springs(ind,3))),x(nint(bending_springs(ind,4))),&
           y(nint(bending_springs(ind,2))),y(nint(bending_springs(ind,3))),y(nint(bending_springs(ind,4))),&
           z(nint(bending_springs(ind,2))),z(nint(bending_springs(ind,3))),z(nint(bending_springs(ind,4))),&
           bending_springs(ind,5),bending_springs(ind,6),bending_springs(ind,7),bending_springs(ind,8),&
           kby(ind), kbz(ind), 3, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))) = &
                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))) + Fxx1

          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,2))+np) = &
          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,2))+np) + Fyy1

      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,2))+2*np) = &
      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))+np) = &
              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))+np) + Fxy1

            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))+2*np) = &
            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,2))) = &
              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,2))) + Fyx1

            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,2))) = &
            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,2))) + Fzx1

                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))) = &
                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))) + Fxx2

          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,3))+np) = &
          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,3))+np) + Fyy2

      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,3))+2*np) = &
      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))+np) = &
              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))+np) + Fxy2

            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))+2*np) = &
            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,3))) = &
              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,3))) + Fyx2

            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,3))) = &
            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,3))) + Fzx2

                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))) = &
                  FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))) + Fxx3

          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,4))+np) = &
          FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,4))+np) + Fyy3

      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,4))+2*np) = &
      FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))+np) = &
              FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))+np) + Fxy3

            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))+2*np) = &
            FJa(nint(bending_springs(ind,3)),nint(bending_springs(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,4))) = &
              FJa(nint(bending_springs(ind,3))+np,nint(bending_springs(ind,4))) + Fyx3

            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,4))) = &
            FJa(nint(bending_springs(ind,3))+2*np,nint(bending_springs(ind,4))) + Fzx3

      !right point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs(ind,2))),x(nint(bending_springs(ind,3))),x(nint(bending_springs(ind,4))),&
           y(nint(bending_springs(ind,2))),y(nint(bending_springs(ind,3))),y(nint(bending_springs(ind,4))),&
           z(nint(bending_springs(ind,2))),z(nint(bending_springs(ind,3))),z(nint(bending_springs(ind,4))),&
           bending_springs(ind,5),bending_springs(ind,6),bending_springs(ind,7),bending_springs(ind,8),&
           kby(ind), kbz(ind), 4, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))) = &
                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))) + Fxx1

          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,2))+np) = &
          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,2))+np) + Fyy1

      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,2))+2*np) = &
      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))+np) = &
              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))+np) + Fxy1

            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))+2*np) = &
            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,2))) = &
              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,2))) + Fyx1

            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,2))) = &
            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,2))) + Fzx1

                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))) = &
                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))) + Fxx2

          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,3))+np) = &
          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,3))+np) + Fyy2

      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,3))+2*np) = &
      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))+np) = &
              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))+np) + Fxy2

            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))+2*np) = &
            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,3))) = &
              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,3))) + Fyx2

            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,3))) = &
            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,3))) + Fzx2

                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))) = &
                  FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))) + Fxx3

          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,4))+np) = &
          FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,4))+np) + Fyy3

      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,4))+2*np) = &
      FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))+np) = &
              FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))+np) + Fxy3

            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))+2*np) = &
            FJa(nint(bending_springs(ind,4)),nint(bending_springs(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,4))) = &
              FJa(nint(bending_springs(ind,4))+np,nint(bending_springs(ind,4))) + Fyx3

            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,4))) = &
            FJa(nint(bending_springs(ind,4))+2*np,nint(bending_springs(ind,4))) + Fzx3
  end do

end subroutine internal_bending_force_derivative

subroutine internal_bending_force_BC_derivative(FJa,x,y,z,x_BC,y_BC,z_BC,bending_springs_BC,np,kby_BC,kbz_BC)
! Calculate Jacobian matrix elements of forces caused by bending springs of veins
! with boundary conditions

  implicit none

  real(kind=pr), intent(in)     :: x(:),y(:),z(:)
  real(kind=pr), intent(in)     :: x_BC(-1:),y_BC(-1:),z_BC(-1:)
  real(kind=pr), intent(in)     :: kby_BC(-1:), kbz_BC(-1:)
  real(kind=pr), intent(in)     :: bending_springs_BC(-1:,1:)
  integer, intent(in)           :: np
  real(kind=pr), intent(inout)  :: FJa(:,:)
  real(kind=pr) :: Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1
  real(kind=pr) :: Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2
  real(kind=pr) :: Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3
  integer :: ind

  do ind=1,nint(maxval(bending_springs_BC(:,1)))

      !left point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs_BC(ind,2))),x(nint(bending_springs_BC(ind,3))),x(nint(bending_springs_BC(ind,4))),&
           y(nint(bending_springs_BC(ind,2))),y(nint(bending_springs_BC(ind,3))),y(nint(bending_springs_BC(ind,4))),&
           z(nint(bending_springs_BC(ind,2))),z(nint(bending_springs_BC(ind,3))),z(nint(bending_springs_BC(ind,4))),&
           bending_springs_BC(ind,5),bending_springs_BC(ind,6),bending_springs_BC(ind,7),bending_springs_BC(ind,8),&
           kby_BC(ind), kbz_BC(ind), 2, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))) = &
                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))) + Fxx1

          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,2))+np) = &
          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,2))+np) + Fyy1

      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,2))+2*np) = &
      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))+np) = &
              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))+np) + Fxy1

            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))+2*np) = &
            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,2))) = &
              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,2))) + Fyx1

            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,2))) = &
            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,2))) + Fzx1

                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))) = &
                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))) + Fxx2

          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,3))+np) = &
          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,3))+np) + Fyy2

      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,3))+2*np) = &
      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))+np) = &
              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))+np) + Fxy2

            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))+2*np) = &
            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,3))) = &
              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,3))) + Fyx2

            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,3))) = &
            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,3))) + Fzx2

                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))) = &
                  FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))) + Fxx3

          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,4))+np) = &
          FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,4))+np) + Fyy3

      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,4))+2*np) = &
      FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))+np) = &
              FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))+np) + Fxy3

            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))+2*np) = &
            FJa(nint(bending_springs_BC(ind,2)),nint(bending_springs_BC(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,4))) = &
              FJa(nint(bending_springs_BC(ind,2))+np,nint(bending_springs_BC(ind,4))) + Fyx3

            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,4))) = &
            FJa(nint(bending_springs_BC(ind,2))+2*np,nint(bending_springs_BC(ind,4))) + Fzx3

      !middle point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs_BC(ind,2))),x(nint(bending_springs_BC(ind,3))),x(nint(bending_springs_BC(ind,4))),&
           y(nint(bending_springs_BC(ind,2))),y(nint(bending_springs_BC(ind,3))),y(nint(bending_springs_BC(ind,4))),&
           z(nint(bending_springs_BC(ind,2))),z(nint(bending_springs_BC(ind,3))),z(nint(bending_springs_BC(ind,4))),&
           bending_springs_BC(ind,5),bending_springs_BC(ind,6),bending_springs_BC(ind,7),bending_springs_BC(ind,8),&
           kby_BC(ind), kbz_BC(ind), 3, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))) = &
                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))) + Fxx1

          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,2))+np) = &
          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,2))+np) + Fyy1

      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,2))+2*np) = &
      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))+np) = &
              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))+np) + Fxy1

            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))+2*np) = &
            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,2))) = &
              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,2))) + Fyx1

            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,2))) = &
            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,2))) + Fzx1

                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))) = &
                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))) + Fxx2

          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,3))+np) = &
          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,3))+np) + Fyy2

      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,3))+2*np) = &
      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))+np) = &
              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))+np) + Fxy2

            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))+2*np) = &
            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,3))) = &
              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,3))) + Fyx2

            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,3))) = &
            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,3))) + Fzx2

                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))) = &
                  FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))) + Fxx3

          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,4))+np) = &
          FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,4))+np) + Fyy3

      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,4))+2*np) = &
      FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))+np) = &
              FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))+np) + Fxy3

            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))+2*np) = &
            FJa(nint(bending_springs_BC(ind,3)),nint(bending_springs_BC(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,4))) = &
              FJa(nint(bending_springs_BC(ind,3))+np,nint(bending_springs_BC(ind,4))) + Fyx3

            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,4))) = &
            FJa(nint(bending_springs_BC(ind,3))+2*np,nint(bending_springs_BC(ind,4))) + Fzx3

      !right point
      call calculate_bending_spring_force_derivative(&
           x(nint(bending_springs_BC(ind,2))),x(nint(bending_springs_BC(ind,3))),x(nint(bending_springs_BC(ind,4))),&
           y(nint(bending_springs_BC(ind,2))),y(nint(bending_springs_BC(ind,3))),y(nint(bending_springs_BC(ind,4))),&
           z(nint(bending_springs_BC(ind,2))),z(nint(bending_springs_BC(ind,3))),z(nint(bending_springs_BC(ind,4))),&
           bending_springs_BC(ind,5),bending_springs_BC(ind,6),bending_springs_BC(ind,7),bending_springs_BC(ind,8),&
           kby_BC(ind), kbz_BC(ind), 4, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))) = &
                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))) + Fxx1

          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,2))+np) = &
          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,2))+np) + Fyy1

      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,2))+2*np) = &
      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,2))+2*np) + Fzz1

              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))+np) = &
              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))+np) + Fxy1

            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))+2*np) = &
            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,2))+2*np) + Fxz1

              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,2))) = &
              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,2))) + Fyx1

            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,2))) = &
            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,2))) + Fzx1

                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))) = &
                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))) + Fxx2

          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,3))+np) = &
          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,3))+np) + Fyy2

      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,3))+2*np) = &
      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,3))+2*np) + Fzz2

              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))+np) = &
              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))+np) + Fxy2

            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))+2*np) = &
            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,3))+2*np) + Fxz2

              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,3))) = &
              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,3))) + Fyx2

            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,3))) = &
            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,3))) + Fzx2

                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))) = &
                  FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))) + Fxx3

          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,4))+np) = &
          FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,4))+np) + Fyy3

      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,4))+2*np) = &
      FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,4))+2*np) + Fzz3

              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))+np) = &
              FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))+np) + Fxy3

            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))+2*np) = &
            FJa(nint(bending_springs_BC(ind,4)),nint(bending_springs_BC(ind,4))+2*np) + Fxz3

              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,4))) = &
              FJa(nint(bending_springs_BC(ind,4))+np,nint(bending_springs_BC(ind,4))) + Fyx3

            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,4))) = &
            FJa(nint(bending_springs_BC(ind,4))+2*np,nint(bending_springs_BC(ind,4))) + Fzx3
  end do

  !---------------------------------------------------------------------------------
  ! Calculate derivative of forces caused by bending springs connecting the vein with the boundary
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

      ! Calculate derivative of force of right part of the bending spring -1 acting
      ! on point 1 (right point correspongding to 4)
      call calculate_bending_spring_force_derivative(&
           x_BC(-1),x_BC(0),x(nint(bending_springs_BC(1,2))), &
           y_BC(-1),y_BC(0),y(nint(bending_springs_BC(1,2))), &
           z_BC(-1),z_BC(0),z(nint(bending_springs_BC(1,2))), &
           bending_springs_BC(-1,5),bending_springs_BC(-1,6), &
           bending_springs_BC(-1,7),bending_springs_BC(-1,8),&
           kby_BC(-1), kbz_BC(-1), 4, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

        !   write(*,*) 'first spring on point 1 BC bending Jacobian Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1'
        !   write(*,*) x_BC(-1),x_BC(0),x(nint(bending_springs_BC(1,2)))
      !     write(*,*) y_BC(-1),y_BC(0),y(nint(bending_springs_BC(1,2)))
        !   write(*,*) z_BC(-1),z_BC(0),z(nint(bending_springs_BC(1,2)))
        !   write(*,*) bending_springs_BC(-1,5),bending_springs_BC(-1,6), &
        !   bending_springs_BC(-1,7),bending_springs_BC(-1,8), kby_BC(-1), kbz_BC(-1)
        !   write(*,*) Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1
        !   write(*,*) Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2
        !   write(*,*) Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3

      ! Adding these terms to the Jacobian matrix
       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))) = &
       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))) + Fxx3

FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))+np) = &
FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))+np) + Fyy3

FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))+2*np) = &
FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))+2*np) + Fzz3

   FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+np) = &
   FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+np) + Fxy3

 FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+2*np) = &
 FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+2*np) + Fxz3

   FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))) = &
   FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))) + Fyx3

 FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))) = &
 FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))) + Fzx3

  ! For bending spring number 0 connecting three points 0,1,2

      ! Calculate derivative of force of left and right parts of the bending spring 0 acting
      ! on point 1 (middle point correspongding to 3)
      call calculate_bending_spring_force_derivative(&
           x_BC(0),x(nint(bending_springs_BC(1,2))),x(nint(bending_springs_BC(1,3))), &
           y_BC(0),y(nint(bending_springs_BC(1,2))),y(nint(bending_springs_BC(1,3))), &
           z_BC(0),z(nint(bending_springs_BC(1,2))),z(nint(bending_springs_BC(1,3))), &
           bending_springs_BC(0,5),bending_springs_BC(0,6), &
           bending_springs_BC(0,7),bending_springs_BC(0,8),&
           kby_BC(0), kbz_BC(0), 3, &
           Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
           Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
           Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

        !   write(*,*) 'second spring on point 1 BC bending Jacobian Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1'
        !   write(*,*) x_BC(0),x(nint(bending_springs_BC(1,2))),x(nint(bending_springs_BC(1,3)))
        !   write(*,*) y_BC(0),y(nint(bending_springs_BC(1,2))),y(nint(bending_springs_BC(1,3)))
        !   write(*,*)  z_BC(0),z(nint(bending_springs_BC(1,2))),z(nint(bending_springs_BC(1,3)))
        !   write(*,*) bending_springs_BC(0,5),bending_springs_BC(0,6), &
        !   bending_springs_BC(0,7),bending_springs_BC(0,8), kby_BC(0), kbz_BC(0)
        !   write(*,*) Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1
        !   write(*,*) Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2
        !   write(*,*) Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3

        ! Adding these terms to the Jacobian matrix
          ! Derivative of middle point with respect to middle point
           FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))) = &
           FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))) + Fxx2

    FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))+np) = &
    FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))+np) + Fyy2

    FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))+2*np) = &
    FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))+2*np) + Fzz2

       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+np) = &
       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+np) + Fxy2

     FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+2*np) = &
     FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,2))+2*np) + Fxz2

       FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))) = &
       FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,2))) + Fyx2

     FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))) = &
     FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,2))) + Fzx2
          ! Derivative of middle point with respect to right point
           FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))) = &
           FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))) + Fxx3

    FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,3))+np) = &
    FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,3))+np) + Fyy3

    FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,3))+2*np) = &
    FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,3))+2*np) + Fzz3

       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))+np) = &
       FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))+np) + Fxy3

     FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))+2*np) = &
     FJa(nint(bending_springs_BC(1,2)),nint(bending_springs_BC(1,3))+2*np) + Fxz3

       FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,3))) = &
       FJa(nint(bending_springs_BC(1,2))+np,nint(bending_springs_BC(1,3))) + Fyx3

     FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,3))) = &
     FJa(nint(bending_springs_BC(1,2))+2*np,nint(bending_springs_BC(1,3))) + Fzx3

      ! Calculate derivative of force of right part of the bending spring 0 acting
      ! on point 2 (right point correspongding to 4)
      call calculate_bending_spring_force_derivative(&
          x_BC(0),x(nint(bending_springs_BC(1,2))),x(nint(bending_springs_BC(1,3))), &
          y_BC(0),y(nint(bending_springs_BC(1,2))),y(nint(bending_springs_BC(1,3))), &
          z_BC(0),z(nint(bending_springs_BC(1,2))),z(nint(bending_springs_BC(1,3))), &
          bending_springs_BC(0,5),bending_springs_BC(0,6), &
          bending_springs_BC(0,7),bending_springs_BC(0,8),&
          kby_BC(0), kbz_BC(0), 4, &
          Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1, &
          Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2, &
          Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3)

        !  write(*,*) 'second spring on point 2 BC bending Jacobian Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1'
        !  write(*,*) x_BC(0),x(nint(bending_springs_BC(1,2))),x(nint(bending_springs_BC(1,3)))
        !  write(*,*) y_BC(0),y(nint(bending_springs_BC(1,2))),y(nint(bending_springs_BC(1,3)))
        !  write(*,*)  z_BC(0),z(nint(bending_springs_BC(1,2))),z(nint(bending_springs_BC(1,3)))
        !  write(*,*) bending_springs_BC(0,5),bending_springs_BC(0,6), &
        !  bending_springs_BC(0,7),bending_springs_BC(0,8), kby_BC(0), kbz_BC(0)
        !  write(*,*) Fxx1,Fyy1,Fzz1,Fxy1,Fxz1,Fyx1,Fzx1
        !  write(*,*) Fxx2,Fyy2,Fzz2,Fxy2,Fxz2,Fyx2,Fzx2
        !  write(*,*) Fxx3,Fyy3,Fzz3,Fxy3,Fxz3,Fyx3,Fzx3

        ! Adding these terms to the Jacobian matrix
          ! Derivative of right point with respect to middle point
           FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))) = &
           FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))) + Fxx2

    FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,2))+np) = &
    FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,2))+np) + Fyy2

    FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,2))+2*np) = &
    FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,2))+2*np) + Fzz2

       FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))+np) = &
       FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))+np) + Fxy2

     FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))+2*np) = &
     FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,2))+2*np) + Fxz2

       FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,2))) = &
       FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,2))) + Fyx2

     FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,2))) = &
     FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,2))) + Fzx2
          ! Derivative of right point with respect to right point
           FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))) = &
           FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))) + Fxx3

    FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,3))+np) = &
    FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,3))+np) + Fyy3

    FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,3))+2*np) = &
    FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,3))+2*np) + Fzz3

       FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))+np) = &
       FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))+np) + Fxy3

     FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))+2*np) = &
     FJa(nint(bending_springs_BC(1,3)),nint(bending_springs_BC(1,3))+2*np) + Fxz3

       FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,3))) = &
       FJa(nint(bending_springs_BC(1,3))+np,nint(bending_springs_BC(1,3))) + Fyx3

     FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,3))) = &
     FJa(nint(bending_springs_BC(1,3))+2*np,nint(bending_springs_BC(1,3))) + Fzx3

end subroutine internal_bending_force_BC_derivative

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
integer, intent(in) :: position
real(kind=pr), intent(out) :: Fxx1, Fyy1, Fzz1, Fxy1, Fxz1, Fyx1, Fzx1
real(kind=pr), intent(out) :: Fxx2, Fyy2, Fzz2, Fxy2, Fxz2, Fyx2, Fzx2
real(kind=pr), intent(out) :: Fxx3, Fyy3, Fzz3, Fxy3, Fxz3, Fyx3, Fzx3
real(kind=pr) :: xy1, xy2, xz1, xz2, xxy1, xxz1, yxy1, zxz1, xxy2, xxz2, yxy2, zxz2
real(kind=pr) :: Ay, By, Az, Bz, kbyT, kbzP

 xy1 = sqrt((x2-x1)**2+(y2-y1)**2)
 xy2 = sqrt((x3-x2)**2+(y3-y2)**2)
 xz1 = sqrt((x2-x1)**2+(z2-z1)**2)
 xz2 = sqrt((x3-x2)**2+(z3-z2)**2)

xxy1 = (x2-x1)/(xy1)**2
xxz1 = (x2-x1)/(xz1)**2
yxy1 = (y2-y1)/(xy1)**2
zxz1 = (z2-z1)/(xz1)**2

xxy2 = (x3-x2)/(xy2)**2
xxz2 = (x3-x2)/(xz2)**2
yxy2 = (y3-y2)/(xy2)**2
zxz2 = (z3-z2)/(xz2)**2

   Ay = ((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))/ &
       (((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))**2+((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))**2)
   By = ((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))/ &
       (((y3-y2)*(y2-y1)+(x3-x2)*(x2-x1))**2+((x3-x2)*(y2-y1)-(y3-y2)*(x2-x1))**2)
   Az = ((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))/ &
       (((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))**2+((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))**2)
   Bz = ((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))/ &
       (((z3-z2)*(z2-z1)+(x3-x2)*(x2-x1))**2+((x3-x2)*(z2-z1)-(z3-z2)*(x2-x1))**2)

kbyT = kby*(theta - theta0)
kbzP = kbz*(phi - phi0)

select case (position)
  case (2) !left point
      !Derivative with respect to the left point X1
        Fxx1 = kby*yxy1*((y3-y2)*Ay + (x3-x2)*By) - 2*kbyT*yxy1*xxy1 + &
               kbz*zxz1*((z3-z2)*Az + (x3-x2)*Bz) - 2*kbzP*zxz1*xxz1

        Fyy1 = kby*xxy1*((x3-x2)*Ay - (y3-y2)*By) + 2*kbyT*yxy1*xxy1

        Fzz1 = kbz*xxz1*((x3-x2)*Az - (z3-z2)*Bz) + 2*kbzP*zxz1*xxz1

        Fxy1 = - kby*yxy1*((x3-x2)*Ay - (y3-y2)*By) + kbyT*(1/(xy1)**2 - 2*(yxy1)**2)

        Fxz1 = - kbz*zxz1*((x3-x2)*Az - (z3-z2)*Bz) + kbzP*(1/(xz1)**2 - 2*(zxz1)**2)

        Fyx1 = - kby*xxy1*((x3-x2)*By + (y3-y2)*Ay) - kbyT*(1/(xy1)**2 - 2*(xxy1)**2)

        Fzx1 = - kbz*xxz1*((x3-x2)*Bz + (z3-z2)*Az) - kbzP*(1/(xz1)**2 - 2*(xxz1)**2)

        !Derivative with respect to the middle point X2
        Fxx2 = - kby*yxy1*((y3-y1)*Ay + (x3-(2*x2)+x1)*By) + 2*kbyT*yxy1*xxy1 &
               - kbz*zxz1*((z3-z1)*Az + (x3-(2*x2)+x1)*Bz) + 2*kbzP*zxz1*xxz1

        Fyy2 = - kby*xxy1*((x3-x1)*Ay - (y3-2*y2+y1)*By) - 2*kbyT*yxy1*xxy1

        Fzz2 = - kbz*xxz1*((x3-x1)*Az - (z3-2*z2+z1)*Bz) - 2*kbzP*zxz1*xxz1

        Fxy2 = kby*yxy1*((x3-x1)*Ay - (y3-2*y2+y1)*By) - kbyT*(1/(xy1)**2 - 2*(yxy1)**2)

        Fxz2 = kbz*zxz1*((x3-x1)*Az - (z3-2*z2+z1)*Bz) - kbzP*(1/(xz1)**2 - 2*(zxz1)**2)

        Fyx2 = kby*xxy1*((y3-y1)*Ay + (x3-2*x2+x1)*By) + kbyT*(1/(xy1)**2 - 2*(xxy1)**2)

        Fzx2 = kbz*xxz1*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + kbzP*(1/(xz1)**2 - 2*(xxz1)**2)

        !Derivative with respect to the third point X3
        Fxx3 = kby*yxy1*((y2-y1)*Ay - (x2-x1)*By) + kbz*zxz1*((z2-z1)*Az - (x2-x1)*Bz)

        Fyy3 = kby*xxy1*((x2-x1)*Ay + (y2-y1)*By)

        Fzz3 = kbz*xxz1*((x2-x1)*Az + (z2-z1)*Bz)

        Fxy3 = - kby*yxy1*((x2-x1)*Ay + (y2-y1)*By)

        Fxz3 = - kbz*zxz1*((x2-x1)*Az + (z2-z1)*Bz)

        Fyx3 = - kby*xxy1*((y2-y1)*Ay - (x2-x1)*By)

        Fzx3 = - kbz*xxz1*((z2-z1)*Az - (x2-x1)*Bz)

    case (3) !middle point
        !Derivative with respect to the left point X1
        Fxx1 = - kby*(yxy1+yxy2)*((y3-y2)*Ay + (x3-x2)*By) + 2*kbyT*yxy1*xxy1 &
               - kbz*(zxz1+zxz2)*((z3-z2)*Az + (x3-x2)*Bz) + 2*kbzP*zxz1*xxz1

        Fyy1 = - kby*(xxy1+xxy2)*((x3-x2)*Ay - (y3-y2)*By) - 2*kbyT*yxy1*xxy1

        Fzz1 = - kbz*(xxz1+xxz2)*((x3-x2)*Az - (z3-z2)*Bz) - 2*kbzP*zxz1*xxz1

        Fxy1 = kby*(yxy1+yxy2)*((x3-x2)*Ay - (y3-y2)*By) + kbyT*(2*(yxy1)**2 - 1/(xy1)**2)

        Fxz1 = kbz*(zxz1+zxz2)*((x3-x2)*Az - (z3-z2)*Bz) + kbzP*(2*(zxz1)**2 - 1/(xz1)**2)

        Fyx1 = kby*(xxy2+xxy1)*((y3-y2)*Ay + (x3-x2)*By) - kbyT*(2*(xxy1)**2 - 1/(xy1)**2)

        Fzx1 = kbz*(xxz1+xxz2)*((z3-z2)*Az + (x3-x2)*Bz) - kbzP*(2*(xxz1)**2 - 1/(xz1)**2)

        !Derivative with respect to the middle point X2
        Fxx2 = kby*(yxy1+yxy2)*((y3-y1)*Ay + (x3-2*x2+x1)*By) + 2*kbyT*(yxy2*xxy2 - yxy1*xxy1) + &
               kbz*(zxz1+zxz2)*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + 2*kbzP*(zxz2*xxz2 - zxz1*xxz1)

        Fyy2 = kby*(xxy2+xxy1)*((x3-x1)*Ay - (y3-2*y2+y1)*By) &
               - 2*kbyT*(yxy2*xxy2 - yxy1*xxy1)

        Fzz2 = kbz*(xxz2+xxz1)*((x3-x1)*Az - (z3-2*z2+z1)*Bz) &
               - 2*kbzP*(zxz2*xxz2 - zxz1*xxz1)

        Fxy2 = - kby*(yxy2+yxy1)*((x3-x1)*Ay - (y3-2*y2+y1)*By) &
               - kbyT*(1/(xy2)**2 - 1/(xy1)**2 - 2*(yxy2)**2 + 2*(yxy1)**2)

        Fxz2 = - kbz*(zxz2+zxz1)*((x3-x1)*Az - (z3-2*z2+z1)*Bz) &
               - kbzP*(1/(xz2)**2 - 1/(xz1)**2 - 2*(zxz2)**2 + 2*(zxz1)**2)

        Fyx2 = - kby*(xxy1+xxy2)*((y3-y1)*Ay + (x3-2*x2+x1)*By) + &
               kbyT*(1/(xy2)**2 - 1/(xy1)**2 + 2*(xxy1)**2 - 2*(xxy2)**2)

        Fzx2 = - kbz*(xxz1+xxz2)*((z3-z1)*Az + (x3-2*x2+x1)*Bz) + &
               kbzP*(1/(xz2)**2 - 1/(xz1)**2 + 2*(xxz1)**2 - 2*(xxz2)**2)

        !Derivative with respect to the third point X3
        Fxx3 = - kby*(yxy1+yxy2)*((y2-y1)*Ay - (x2-x1)*By) - 2*kbyT*yxy2*xxy2 &
               - kbz*(zxz1+zxz2)*((z2-z1)*Az - (x2-x1)*Bz) - 2*kbzP*zxz2*xxz2

        Fyy3 = - kby*(xxy1+xxy2)*((x2-x1)*Ay + (y2-y1)*By) + 2*kbyT*yxy2*xxy2

        Fzz3 = - kbz*(xxz1+xxz2)*((x2-x1)*Az + (z2-z1)*Bz) + 2*kbzP*zxz2*xxz2

        Fxy3 = kby*(yxy1+yxy2)*((x2-x1)*Ay + (y2-y1)*By) - kbyT*(2*(yxy2)**2 - 1/(xy2)**2)

        Fxz3 = kbz*(zxz1+zxz2)*((x2-x1)*Az + (z2-z1)*Bz) - kbzP*(2*(zxz2)**2 - 1/(xz2)**2)

        Fyx3 = kby*(xxy2+xxy1)*((y2-y1)*Ay - (x2-x1)*By) + kbyT*(2*(xxy2)**2 - 1/(xy2)**2)

        Fzx3 = kbz*(xxz1+xxz2)*((z2-z1)*Az - (x2-x1)*Bz) + kbzP*(2*(xxz2)**2 - 1/(xz2)**2)

    case (4) !right point
       !Derivative with respect to the left point X1
        Fxx1 = kby*yxy2*((y3-y2)*Ay + (x3-x2)*By) + &
               kbz*zxz2*((z3-z2)*Az + (x3-x2)*Bz)

        Fyy1 = kby*xxy2*((x3-x2)*Ay - (y3-y2)*By)

        Fzz1 = kbz*xxz2*((x3-x2)*Az - (z3-z2)*Bz)

        Fxy1 = - kby*yxy2*((x3-x2)*Ay - (y3-y2)*By)

        Fxz1 = - kbz*zxz2*((x3-x2)*Az - (z3-z2)*Bz)

        Fyx1 = - kby*xxy2*((y3-y2)*Ay + (x3-x2)*By)

        Fzx1 = - kbz*xxz2*((z3-z2)*Az + (x3-x2)*Bz)

        !Derivative with respect to the middle point X2
        Fxx2 = - kby*yxy2*((y3-y1)*Ay + (x3-2*x2+x1)*By) - 2*kbyT*yxy2*xxy2 &
               - kbz*zxz2*((z3-z1)*Az + (x3-2*x2+x1)*Bz) - 2*kbzP*zxz2*xxz2

        Fyy2 = - kby*xxy2*((x3-x1)*Ay - (y3-2*y2+y1)*By) + 2*kbyT*yxy2*xxy2

        Fzz2 = - kbz*xxz2*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + 2*kbzP*zxz2*xxz2

        Fxy2 = kby*yxy2*((x3-x1)*Ay - (y3-2*y2+y1)*By) + kbyT*(1/(xy2)**2 - 2*(yxy2)**2)

        Fxz2 = kbz*zxz2*((x3-x1)*Az - (z3-2*z2+z1)*Bz) + kbzP*(1/(xz2)**2 - 2*(zxz2)**2)

        Fyx2 = kby*xxy2*((y3-y1)*Ay + (x3-2*x2+x1)*By) - kbyT*(1/(xy2)**2 - 2*(xxy2)**2)

        Fzx2 = kbz*xxz2*((z3-z1)*Az + (x3-2*x2+x1)*Bz) - kbzP*(1/(xz2)**2 - 2*(xxz2)**2)


        !Derivative with respect to the third point X3
        Fxx3 = kby*yxy2*((y2-y1)*Ay - (x2-x1)*By) + 2*kbyT*yxy2*xxy2 + &
               kbz*zxz2*((z2-z1)*Az - (x2-x1)*Bz) + 2*kbzP*zxz2*xxz2

        Fyy3 = kby*xxy2*((x2-x1)*Ay + (y2-y1)*By) - 2*kbyT*yxy2*xxy2

        Fzz3 = kbz*xxz2*((x2-x1)*Az + (z2-z1)*Bz) - 2*kbzP*zxz2*xxz2

        Fxy3 = - kby*yxy2*((x2-x1)*Ay + (y2-y1)*By) - kbyT*(1/(xy2)**2 - 2*(yxy2)**2)

        Fxz3 = - kbz*zxz2*((x2-x1)*Az + (z2-z1)*Bz) - kbzP*(1/(xz2)**2 - 2*(zxz2)**2)

        Fyx3 = - kby*xxy2*((y2-y1)*Ay - (x2-x1)*By) + kbyT*(1/(xy2)**2 - 2*(xxy2)**2)

        Fzx3 = - kbz*xxz2*((z2-z1)*Az - (x2-x1)*Bz) + kbzP*(1/(xz2)**2 - 2*(xxz2)**2)

end select

end subroutine
