! Wing wrapper for different wing shapes
subroutine draw_wing(mask,mask_color,us,Insect,color_wing,M_body,M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)

  real(kind=pr) :: t1
  t1 = MPI_wtime()

  select case(Insect%WingShape)
  case ("mosquito_iams")
    call draw_wing_mosquito(mask,mask_color,us,Insect,color_wing,M_body,&
         M_wing,x_pivot,rot)
  case ("rectangular")
    call draw_wing_rectangular(mask,mask_color,us,Insect,color_wing,M_body,&
         M_wing,x_pivot,rot)
  case ("suzuki")
    ! this wing has a finite thickness
    call draw_wing_suzuki(mask,mask_color,us,Insect,color_wing,M_body,&
         M_wing,x_pivot,rot)
  case ("TwoEllipses")
    call draw_wing_twoellipses(mask,mask_color,us,Insect,color_wing,M_body,&
         M_wing,x_pivot,rot)
  case default
    ! we assume the default to be defined in fourier coefficients, the subroutine
    ! yells if it does not recongnize the wing.
    call draw_wing_fourier(mask,mask_color,us,Insect,color_wing,M_body,M_wing,&
         x_pivot,rot)
  end select

  time_insect_wings = time_insect_wings + MPI_wtime() - t1
end subroutine draw_wing

!-------------------------------------------------------------------------------

! Draws a wings that is given by a radius(theta), where the radius is given
! by a Fourier series. The Fourier coefficients are stored in the insect
! datastructure, so the function Set_Wing_Fourier_coefficients must be called
! before calling this subroutine. Fourier series is evaluated in
! Radius_Fourier
subroutine draw_wing_fourier(mask,mask_color,us,Insect,color_wing,M_body,&
           M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=pr) :: R, R0, R_tmp
  real(kind=pr) :: y_tmp, x_tmp, z_tmp
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta

  !-- define the wings fourier coeffients, but only once
  call Setup_Wing_Fourier_coefficients(Insect)


  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the various coordinate systems we are going to use
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x = periodize_coordinate(x - Insect%xc_body)
        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot)

        !-- first, check if the point lies inside the rectangle L_span x L_span
        !-- here we assume that the chordlength is NOT greater than the span
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
        if ((x_wing(1)>=-(Insect%L_span+Insect%safety)).and.(x_wing(1)<=Insect%L_span+Insect%safety)) then
        if (abs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then

          !-- get normalized angle (theta)
          theta = atan2 (x_wing(2)-Insect%yc,x_wing(1)-Insect%xc )
          theta = ( theta + pi ) / (2.d0*pi)

          !-- construct R by evaluating the fourier series
          R0 = Radius_Fourier(theta,Insect)

          !-- get smooth (radial) step function
          R = dsqrt ( (x_wing(1)-Insect%xc)**2 + (x_wing(2)-Insect%yc)**2 )
          R_tmp = steps(R,R0)

          !-- smooth also the thicknes
          z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness
          mask_tmp = z_tmp*R_tmp

          !-----------------------------------------
          ! set new value for mask and velocity us
          !-----------------------------------------
          if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then
            mask(ix,iy,iz) = mask_tmp
            mask_color(ix,iy,iz) = color_wing
            !------------------------------------------------
            ! solid body rotation
            ! Attention: the Matrix transpose(M) brings us back to the body
            ! coordinate system, not to the inertial frame. this is done in
            ! the main routine Draw_Insect
            !------------------------------------------------
            v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
            v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
            v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)
            ! note we set this only if it is a part of the wing
            us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
          endif
        endif
        endif
        endif

      enddo
    enddo
  enddo

end subroutine draw_wing_fourier

!-------------------------------------------------------------------------------
! Suzuki's rectangular wingshape as defined in section 3.4 of my thesis (Thomas,
! "Numerical modeling of fluid-structure interaction in bio-inspired propulsion")
! This wing has finite thickness.
!-------------------------------------------------------------------------------
subroutine draw_wing_suzuki(mask,mask_color,us,Insect,color_wing,M_body,&
           M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=pr) :: R, R0, R_tmp
  real(kind=pr) :: y_tmp, x_tmp, z_tmp, y_left, y_right
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta, x_top, x_bot

  ! wing shape (determine between which x-values (x_bot, x_top) the wing is
  ! these values depend on the spanwise direction (which is y)
  x_top = 0.0667d0
  x_bot = -0.35d0

  y_right = 1.0d0
  y_left = 0.1667d0

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the various coordinate systems we are going to use
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x = periodize_coordinate(x - Insect%xc_body)
        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot)

        ! spanwise length:
        if ((x_wing(2)>=y_left-Insect%safety).and.(x_wing(2)<=y_right+Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if ((x_wing(2)-y_left)<0.d0) then ! xs is chordlength coordinate
                y_tmp = steps(-(x_wing(2)-y_left),0.d0)
              else
                y_tmp = steps( (x_wing(2)-y_left),y_right-y_left)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot)
              else
                x_tmp = steps( x_wing(1), x_top)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
                v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
                v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_suzuki

!-------------------------------------------------------------------------------

subroutine draw_wing_rectangular(mask,mask_color,us,Insect,color_wing,M_body,&
           M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=pr) :: R, R0, R_tmp
  real(kind=pr) :: y_tmp, x_tmp, z_tmp
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot


  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the various coordinate systems we are going to use
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x = periodize_coordinate(x - Insect%xc_body)
        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            x_top = Insect%b_top
            x_bot =-Insect%b_bot
            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                y_tmp = steps(-x_wing(2),0.d0)
              else
                y_tmp = steps( x_wing(2),Insect%L_span)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot)
              else
                x_tmp = steps( x_wing(1), x_top)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
                v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
                v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_rectangular


!-------------------------------------------------------------------------------
! Draws a wing
! here, a wing is a rigid plate of constant thickness that differs from
! a rectangular plate only in the x-direction
!
! note to save a bit of computing time, we first check the easy
! conditions (thickness and spanwise length) and then the shape
! function since this saves many evaluations of the shape.
subroutine draw_wing_twoellipses(mask,mask_color,us,Insect,color_wing,M_body,&
           M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)


  integer :: ix,iy,iz
  real(kind=pr) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=pr) :: R, R0, R_tmp,a_body
  real(kind=pr) :: y_tmp, x_tmp, z_tmp
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  a_body = 0.5d0 * Insect%L_span

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the various coordinate systems we are going to use
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x = periodize_coordinate(x-Insect%xc_body)
        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            if ((1.d0 - ((x_wing(2)-a_body)**2)/(a_body**2)) >= 0.d0) then
              x_top =  dsqrt((Insect%b_top**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
              x_bot = -dsqrt((Insect%b_bot**2)*(1.d0-((x_wing(2)-a_body)**2)/(a_body**2)))
            else
              x_top = 0.d0
              x_bot = 0.d0
            endif

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                y_tmp = steps(-x_wing(2),0.d0)
              else
                y_tmp = steps( x_wing(2),Insect%L_span)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot)
              else
                x_tmp = steps( x_wing(1), x_top)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
                v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
                v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_twoellipses



!-------------------------------------------------------------------------------
! Draws a wing of a mosquito. it is a simple ellipse shape, as presented in
! [1] Iams "Flight stability of mosquitos: A reduced model" SIAM J. Appl. Math. 74(5) 1535--1550 (2014)
subroutine draw_wing_mosquito(mask,mask_color,us,Insect,color_wing,M_body,&
           M_wing,x_pivot,rot)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_wing
  real(kind=pr),intent(in)::M_body(1:3,1:3),M_wing(1:3,1:3),x_pivot(1:3),rot(1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: x_body(1:3),x_wing(1:3),x(1:3)
  real(kind=pr) :: R, R0, R_tmp,a_wing
  real(kind=pr) :: y_tmp, x_tmp, z_tmp
  real(kind=pr) :: v_tmp(1:3), mask_tmp, theta,x_top,x_bot

  a_wing = 0.5d0 * Insect%L_span
  Insect%b_top = 0.1474d0
  Insect%b_bot = 0.1474d0

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the various coordinate systems we are going to use
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x = periodize_coordinate(x-Insect%xc_body)
        x_body = matmul(M_body,x)
        x_wing = matmul(M_wing,x_body-x_pivot)

        ! spanwise length:
        if ((x_wing(2)>=-Insect%safety).and.(x_wing(2)<=Insect%L_span + Insect%safety)) then
          ! thickness: (note left and right wing have a different orientation of the z-axis
          ! but this does not matter since this is the same.
          if (dabs(x_wing(3))<=0.5*Insect%WingThickness + Insect%safety) then
            ! wing shape (determine between which x-values (x_bot, x_top) the wing is
            ! these values depend on the spanwise direction (which is y)
            if ((1.d0 - ((x_wing(2)-a_wing)**2)/(a_wing**2)) >= 0.d0) then
              x_top =  dsqrt((Insect%b_top**2)*(1.d0-((x_wing(2)-a_wing)**2)/(a_wing**2)))
              x_bot = -dsqrt((Insect%b_bot**2)*(1.d0-((x_wing(2)-a_wing)**2)/(a_wing**2)))
            else
              x_top = 0.d0
              x_bot = 0.d0
            endif

            ! in the x-direction, the actual wing shape plays.
            if ((x_wing(1)>x_bot-Insect%safety).and.(x_wing(1)<x_top+Insect%safety)) then
              !-- smooth length
              if (x_wing(2)<0.d0) then  ! xs is chordlength coordinate
                y_tmp = steps(-x_wing(2),0.d0)
              else
                y_tmp = steps( x_wing(2),Insect%L_span)
              endif

              !-- smooth height
              z_tmp = steps(dabs(x_wing(3)),0.5d0*Insect%WingThickness) ! thickness

              !-- smooth shape
              if (x_wing(1)<0.d0) then
                x_tmp = steps(-x_wing(1),-x_bot)
              else
                x_tmp = steps( x_wing(1), x_top)
              endif

              mask_tmp = z_tmp*y_tmp*x_tmp

              if ((mask(ix,iy,iz) <= mask_tmp).and.(mask_tmp>0.0)) then
                mask(ix,iy,iz) = mask_tmp
                mask_color(ix,iy,iz) = color_wing
                !------------------------------------------------
                ! solid body rotation
                ! Attention: the Matrix transpose(M) brings us back to the body
                ! coordinate system, not to the inertial frame. this is done in
                ! the main routine Draw_Insect
                !------------------------------------------------
                v_tmp(1) = rot(2)*x_wing(3)-rot(3)*x_wing(2)
                v_tmp(2) = rot(3)*x_wing(1)-rot(1)*x_wing(3)
                v_tmp(3) = rot(1)*x_wing(2)-rot(2)*x_wing(1)

                ! note we set this only if it is a part of the wing
                us(ix,iy,iz,1:3) = matmul(transpose(M_wing), v_tmp)
              endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
end subroutine draw_wing_mosquito



!-------------------------------------------------------------------------------
! evaluates the fourier series given in the ai, bi
!-------------------------------------------------------------------------------
real(kind=pr) function Radius_Fourier(theta,Insect)
  use vars
  implicit none
  integer :: i,j, n_radius
  real(kind=pr) :: R0, theta2, dphi
  type(diptera),intent(inout)::Insect
  real(kind=pr), intent(in) :: theta

  n_radius = 25000
  dphi = (2.d0*pi) / (dble(n_radius-1))

  ! evaluate the entire R(theta) once with very fine resolution, so when
  ! calling it for the second time we only need linear interpolation.
  if (.not.Insect%wings_radius_table_ready) then
    ! loop over all thetas
    do j = 1, n_radius
      R0 = Insect%a0_wings / 2.d0
      theta2 = dble(j-1) * dphi
      ! evaluate Fourier series
      do i = 1, Insect%nfft_wings
        R0 = R0 + Insect%ai_wings(i)*dcos(2.d0*pi*dble(i)*theta2) &
                + Insect%bi_wings(i)*dsin(2.d0*pi*dble(i)*theta2)
      enddo
      Insect%R0_table(j)=R0
    enddo
    ! skip setup on next call
    Insect%wings_radius_table_ready = .true.
  endif

  ! linear interpolation, if already stored the radius
  j = floor( theta / dphi ) + 1
  Radius_Fourier = Insect%R0_table(j) + ((theta-dble(j-1)*dphi) / dphi) &
                 * (Insect%R0_table(j+1)-Insect%R0_table(j))
end function

!-------------------------------------------------------------------------------
! Here all hard-coded fourier series coefficients for different wings shapes are
! collected. This routine is only called once per time step, and it doesn't do
! anything when called for the second time.
! In the first call, the arrays ai and bi, that hold the Fourier
! coefficients, are allocated. Then they are filled with the values corresponding
! to Insect%WingShape. If the routine is called with an unkown wing shape, it
! stops the code. This prevents errors for wings that are NOT given by Fourier
! series.
!-------------------------------------------------------------------------------
subroutine Setup_Wing_Fourier_coefficients(Insect)
  use vars
  implicit none
  real(kind=pr) :: xroot, yroot
  type(diptera),intent(inout)::Insect

  if (Insect%wingsetup_done) then
    ! the second call is just a return statement
    return
  endif

  !-----------------------------------------
  ! hard-coded Fourier coefficients for R(theta)
  !-----------------------------------------
  select case (Insect%WingShape)
  case ('drosophila')
    !********************************************
    ! Drosophila wing from Jan Gruber's png file
    !********************************************
    Insect%nfft_wings = 40
    Insect%a0_wings = 0.5140278
    Insect%ai_wings(1:Insect%nfft_wings) = &
        (/0.1276258,-0.1189758,-0.0389458,0.0525938,0.0151538,-0.0247938,&
          -0.0039188,0.0104848,-0.0030638,-0.0064578,0.0042208,0.0043248,&
          -0.0026878,-0.0021458,0.0017688,0.0006398,-0.0013538,-0.0002038,&
          0.0009738,0.0002508,-0.0003548,-0.0003668,-0.0002798,0.0000568,&
          0.0003358,0.0001408,-0.0002208,0.0000028,0.0004348,0.0001218,&
          -0.0006458,-0.0003498,0.0007168,0.0003288,-0.0007078,-0.0001368,&
          0.0007828,0.0001458,-0.0007078,-0.0001358/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
        (/-0.1072518,-0.0449318,0.0296558,0.0265668,-0.0043988,-0.0113218,&
          -0.0003278,0.0075028,0.0013598,-0.0057338,-0.0021228,0.0036178,&
          0.0013328,-0.0024128,-0.0007688,0.0011478,0.0003158,-0.0005528,&
          0.0000458,0.0003768,0.0002558,0.0000168,-0.0006018,-0.0006338,&
          0.0001718,0.0007758,0.0001328,-0.0005888,-0.0001088,0.0006298,&
          0.0000318,-0.0008668,-0.0000478,0.0009048,0.0001198,-0.0008248,&
          -0.0000788,0.0007028,-0.0000118,-0.0006608/)

    ! wing root point
    xroot =+0.1122
    yroot =-0.0157
    ! center of circle
    Insect%xc =-0.1206 + xroot
    Insect%yc = 0.3619 + yroot
  case ('drosophila_mutated')
    !********************************************
    ! mutated Drosophila wing from Jan Gruber's png file
    !********************************************
    Insect%nfft_wings = 70
    Insect%a0_wings = 0.4812548
    Insect%ai_wings(1:Insect%nfft_wings) = &
           (/0.1593968, -0.1056828, -0.0551518, 0.0508748, 0.0244538, -0.0264738,&
            -0.0080828, 0.0181228, 0.0023648, -0.0134578, -0.0037068, 0.0064508,&
            0.0028748, -0.0014258, -0.0006028, -0.0008898, -0.0020408, 0.0009218,&
            0.0029938, 0.0002768, -0.0026968, -0.0011518, 0.0017798, 0.0016538,&
            -0.0006098, -0.0012998, -0.0001918, 0.0003478, 0.0001408, 0.0003098,&
            0.0001078, -0.0005568, -0.0005998, 0.0006128, 0.0009078, -0.0003798,&
            -0.0009268, 0.0002128, 0.0009098, -0.0000598, -0.0010668, -0.0003428,&
            0.0009228, 0.0007688, -0.0003568, -0.0010458, -0.0004378, 0.0008738,&
            0.0009478, -0.0004108, -0.0012248, -0.0000638, 0.0013148, 0.0004978,&
            -0.0010638, -0.0007148, 0.0006338, 0.0007438, -0.0003278, -0.0006078,&
            0.0001838, 0.0003768, -0.0001698, -0.0002148, 0.0001318, 0.0001628,&
            -0.0000878, 0.0000068, 0.0001478, -0.0001128/)

    Insect%bi_wings(1:Insect%nfft_wings) = &
          (/-0.1132588, -0.0556428, 0.0272098, 0.0221478, -0.0063798, -0.0059078,&
              0.0043788, 0.0043208, -0.0003308, -0.0026598, -0.0013158, 0.0025178,&
              0.0022438, -0.0023798, -0.0037048, 0.0001528, 0.0031218, 0.0022248,&
              -0.0007428, -0.0027298, -0.0018298, 0.0014538, 0.0028888, 0.0000648,&
              -0.0023508, -0.0009418, 0.0017848, 0.0016578, -0.0008058, -0.0017348,&
              -0.0001368, 0.0011138, 0.0004218, -0.0005918, -0.0002798, 0.0002388,&
              0.0002148, 0.0001408, 0.0000218, -0.0005138, -0.0003458, 0.0008208,&
              0.0009888, -0.0007468, -0.0015298, 0.0002728, 0.0015588, 0.0002758,&
              -0.0012498, -0.0006908,0.0008718, 0.0008848, -0.0003038, -0.0008048,&
              -0.0001538, 0.0005418, 0.0003658, -0.0001988, -0.0003938, 0.0000048,&
              0.0003008, 0.0000538, -0.0002748, -0.0000598, 0.0002898, 0.0001398,&
              -0.0002108, -0.0001888, 0.0001838, 0.0001888 /)

    ! wing root point
    xroot =+0.1122
    yroot =-0.0157
    ! center of circle
    Insect%xc =-0.1206 + xroot
    Insect%yc = 0.3619 + yroot
  case ('drosophila_sandberg')
    !********************************************
    !  Drosophila wing from Ramamurti & Sandberg ( JEB 210, 881-896, 2007)
    !********************************************
    Insect%nfft_wings = 24
    Insect%a0_wings = 0.4995578
    Insect%ai_wings(1:Insect%nfft_wings) = &
    (/0.0164168,-0.1621518,0.0030938,0.0601108,-0.0083988,-0.0199988,&
    0.0049048,0.0047878,-0.0005648,-0.0001108,-0.0008638,-0.0006928,&
    0.0006608,0.0001978,0.0001558,0.0006878,-0.0007498,-0.0008018,&
    0.0003878,0.0007028,0.0000408,-0.0001108,-0.0001068,-0.0003958 &
    /)
    Insect%bi_wings(1:Insect%nfft_wings) = &
    (/-0.2083518,-0.0106488,0.0878308,-0.0018168,-0.0338278,0.0045768,&
    0.0113778,-0.0020678,-0.0026928,0.0002758,-0.0000838,-0.0001298,&
    0.0004118,0.0005638,-0.0001018,-0.0006918,-0.0002268,0.0005238,&
    0.0004008,-0.0001818,-0.0003038,-0.0000068,-0.0001218,0.0002008 &
    /)
    Insect%xc =-0.0235498
    Insect%yc = 0.1531398
  case ('drosophila_maeda')
    !********************************************
    !  Drosophila wing from Maeda and Liu, similar to Liu and Aono, BB2009
    !********************************************
    Insect%nfft_wings = 25
    !Insect%a0_wings = 0.591294836514357
    !ai = (/0.11389995408864588, -0.08814321795213981, -0.03495210456149335,&
    !0.024972085605453047, 0.009422293191002384, -0.01680813499169695,&
    !-0.006006435254421029, 0.012157932943676907, 0.00492283934032996,&
    !-0.009882103857127606, -0.005421102356676356, 0.007230876076797827,&
    !0.005272314598249222, -0.004519437431722127, -0.004658072133773225,&
    !0.0030795046767766853, 0.003970792618725898, -0.0016315879319092456,&
    !-0.002415442110272326, 0.0011118187761994598, 0.001811261693911865,&
    !-2.6496695842951815E-4, -0.0012472769174353662, -1.7427507835680091E-4,&
    !0.0010049640224536927/)
    !bi = (/0.0961275426181888, 0.049085916171592914, -0.022051083533094627,&
    !-0.014004783021121204, 0.012955446778711292, 0.006539648525493488,&
    !-0.011873438993933363, -0.00691719567010525, 0.008479044683798266,&
    !0.0045388280405204194, -0.008252172088956379, -0.005091347100627815,&
    !0.004626409662755484, 0.004445034936616318, -0.0030708884306814804,&
    !-0.004428808427471962, 0.0014113707529017868, 0.003061279043478891,&
    !-8.658653756413232E-4, -0.002153349816945423, 3.317570161883452E-4,&
    !0.001573518502682025, 2.14583094242007E-4, -0.0011299834277813852,&
    !-5.172854674801216E-4/)
    Insect%a0_wings = 0.585432698694358
    Insect%ai_wings(1:Insect%nfft_wings) = &
    (/0.113400475583443, -0.0862823485047213, -0.0346234482214816,&
    0.0237625254732323,0.00902498439287132,-0.0158926757445186,&
    -0.00549384372979449,0.0114928668063701,0.00431222381497978,&
    -0.00951270119733201,-0.00484045133879639,0.00706223174320460,&
    0.00473736389439926,-0.00449539769983697,-0.00418487169011745,&
    0.00320520052884641,0.00355631891057573,-0.00183155403463614,&
    -0.00191680264797099,0.00144768631289857,0.00135580122365068,&
    -0.000579638217642394,-0.000818378434108882,0.000132570375969864,&
    0.000683325977327827/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/0.0939265226506824,0.0486063180327962,-0.0206591129298861,&
    -0.0136085709392758,0.0118575265347540,0.00604510770670991,&
    -0.0110263907282936,-0.00636979352727611,0.00786779216718321,&
    0.00390493804324433,-0.00797763174406198,-0.00450123591642554,&
    0.00445099872769504,0.00387237248613979,-0.00305464314668877,&
    -0.00398381251524846,0.00144450353105449,0.00257445316700965,&
    -0.00104247508055041,-0.00167946127380679,0.000577428923826108,&
    0.00114016779684690,-2.63209684213992e-05,-0.000753899380930065,&
    -0.000294894042986087/)

    !Insect%xc = 0.0 ! original mesh
    Insect%xc = 0.0473 ! shifted towards t.e. to 1/4 of the root chord ("+" sign here)
    !Insect%xc = -0.0728 ! shifted towards l.e., to 0.2cmean from the l.e. (Liu and Aono BB 2009)
    !Insect%yc = 0.7
    !Insect%yc = 0.712 ! measured using kinematics snapshots
    Insect%yc = 0.702 ! According to Maeda's email, Jun 21, 2014
  case ('drosophila_sun')
    !********************************************
    !  Drosophila virilis wing from Chen and Sun, Acta Mech Sin 2014
    !********************************************
    Insect%nfft_wings = 25
    Insect%a0_wings = 0.5427789795180327
    Insect%ai_wings(1:Insect%nfft_wings) = &
    (/0.10879717599747996, -0.11445383382313232, -0.02023898255134523,&
    0.04903268079573884, 0.0012813019346402, -0.02397317767942499,&
    0.0013575396713610029, 0.0108149787395804, -0.001514114743464855,&
    -0.005364275911068656, 3.6505751634048205E-4, 0.002640180169907162,&
    -3.2673259786225535E-4, -0.0014323857426683313, 2.431115176929324E-4,&
    5.392319229992534E-4, -4.5833334881866856E-4, -1.3216432233072333E-4,&
    6.563502263270568E-4, 2.0750829321817808E-4, -4.807960800434886E-4,&
    -2.9006261005712504E-4, 2.7578746591965946E-4, 2.7519915193569374E-4,&
    -3.0570604954113513E-4/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/-0.09385487982296374, -0.010821846776797858, 0.030052970821579587,&
    0.005312859230387492, -0.006054695188204192, -0.0015421303479118118,&
    -0.002533264559802815, -0.0014806147599133366, 0.003640199794653037,&
    0.0020416212413267134, -0.0024948946435721206, -7.83017244422372E-4,&
    0.0021574389122894035, 2.6950667683726845E-4, -0.00131044444112179,&
    6.404390762251693E-5, 2.513250728448789E-4, -4.7634735716375334E-4,&
    -1.5949800516545527E-5, 5.001276053841919E-4, 8.445613796483002E-5,&
    -5.510759077970704E-4, -3.3722093938416713E-4, 3.524656540450335E-4,&
    2.9924999100355387E-4/)
    Insect%xc = 0.0
    Insect%yc = 0.399446382250523
  case ('bumblebee')
    !********************************************
    !  Bumblebee
    !  http://www.entomology.umn.edu/museum/links/coursefiles/JPEG%20images/Hymenoptera%20web%20jpeg/Bombus-wings.jpg
    !********************************************
    Insect%nfft_wings = 25
    Insect%a0_wings = 0.594557593733011d0
    Insect%ai_wings(1:Insect%nfft_wings) = &
     (/-0.0128037920989526d0,-0.106777418654552d0,0.0380851289321982d0,&
    0.0330548081081197d0,-0.0178496286355627d0,-0.00328588543359649d0,&
    0.0108246924336137d0,-0.00489302388329943d0,-0.00708808441172961d0,&
    0.00518244772891516d0,0.00445979960844562d0,-0.000108072056165527d0,&
    0.00204437845603716d0,0.00147176382618797d0,-0.00229559463098105d0,&
    -0.000514633972391526d0,0.00134150515430486d0,-0.000149860228261824d0,&
    9.01456938813568d-05,0.00150639712261487d0,0.000914624010720407d0,&
    -0.000737650894315551d0,-0.000843849923321745d0,-0.000354971670482499d0,&
    -0.000382956472432449d0/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/-0.0158061138788171d0,0.0308243584184200d0,-0.00903330410923372d0,&
    -0.0185758334697500d0,-0.000924452934252486d0,-0.00242101213359519d0,&
    -0.00204549530064489d0,0.00291468131401423d0,-0.000140755032337495d0,&
    -0.00135036427128534d0,0.00141285439042451d0,-0.000334215276598231d0,&
    -0.00161521722061879d0,-0.000164055684312904d0,-0.000256278551727569d0,&
    -0.000740258481681094d0,0.000847498161852221d0,0.00157442110960973d0,&
    -0.000559835622451578d0,-0.000617498559228280d0,0.00115413452523474d0,&
    0.000322564770099778d0,-0.000917375185844477d0,4.44819399488798d-05,&
    0.000710028654602170d0/)
    Insect%xc = -0.1d0
    Insect%yc = 0.501549263807117d0
  case ('b_ignitus')
    !********************************************
    !  Bumblebee B. ignitus
    !  Digitized from images taken at Liu Lab
    !********************************************
    Insect%nfft_wings = 25
    Insect%a0_wings = 0.536472532931637d0
    Insect%ai_wings(1:Insect%nfft_wings) = &
    (/-0.0447167394708177d0,-0.106357727795917d0,0.0504418160417239d0,&
    0.0217275689429364d0,-0.0259085955164794d0,0.00272535910748833d0,&
    0.00925289824790763d0,-0.00453010629382665d0,-0.000726647749565597d0,&
    0.00258280099999843d0,-0.00193033529765617d0,-0.00121090519402499d0,&
    0.00149872968653121d0,0.000716207684720514d0,-0.000205317764190544d0,&
    0.000120507537444963d0,-0.000381477942805165d0,-0.000364957961985063d0,&
    -6.70598716926467d-05,0.000166365788794039d0,0.000332993591840758d0,&
    -0.000225912231239784d0,-0.000554023819155716d0,0.000352735383706648d0,&
    0.000650085631908143d0/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/-0.0580660125663764d0,0.0271775529659247d0,0.0178916506228727d0,&
    -0.0196983386855655d0,-0.00865040473524334d0,0.0112078637630294d0,&
    0.00505882127179290d0,-0.00516874871678530d0,-0.000418585234573997d0,&
    0.00248996756589669d0,-0.00248081765717699d0,-0.00165307115885468d0,&
    0.00236884835642553d0,0.000920860396041608d0,-0.00160449459432319d0,&
    7.96078949775159d-05,0.000716588388745441d0,0.000306756717543478d0,&
    0.000310638954298390d0,-0.000523512353114016d0,-0.000773372382092419d0,&
    1.97258594500968d-05,0.000261943571939630d0,0.000262003935722642d0,&
    0.000278542046262820d0/)
    Insect%xc = -0.13d0
    Insect%yc = 0.434820393790595d0
  case ('flapper_sane')
    !********************************************
    !  Mechanical model from Sane and Dickinson, JEB 205, 2002
    !  'The aerodynamic effects...'
    !********************************************
    Insect%nfft_wings = 25
    Insect%a0_wings = 0.5379588906565078
    Insect%ai_wings(1:Insect%nfft_wings) = &
     (/0.135338653455782,-0.06793162622123261,-0.0398235167675977,&
    0.006442194893963269,0.0012783260416583853,-0.007014398516674715,&
    0.0017710765408983137,0.006401601802033519,-2.970619204124993E-4,&
    -0.0038483478773981405,-6.180958756568494E-4,8.015784831786756E-4,&
    -6.957513357109226E-4,-1.4028929172227943E-4,0.0013484885717868547,&
    4.827827498543977E-4,-9.747844462919694E-4,-5.838504331939134E-4,&
    2.72834004831554E-4,2.8152492682871664E-5,-1.2802199282558645E-4,&
    4.117887216124469E-4,3.364169982438278E-4,-3.33258003686823E-4,&
    -3.5615733035757616E-4/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/2.686408368800394E-4,0.01649582345310688,0.01288513083639708,&
    0.004711436946785864,-0.0035725088809005073,-0.00898640397179334,&
    -0.003856509905612652,0.004536524572892801,0.004849677692836578,&
    2.9194421255236984E-4,-7.512780802871473E-4,7.12685261783966E-4,&
    -1.5519932673320404E-4,-0.0012695469974603026,2.2861692091158138E-4,&
    0.0016461316319681953,5.257476721137781E-4,-7.686482830046961E-4,&
    -3.108879176661735E-4,2.2437540206568518E-4,-2.578427217327782E-4,&
    -2.5120263516966855E-4,4.1693453021778877E-4,3.9290173948150096E-4,&
    -1.9762601237675826E-4/)
    Insect%xc = 0.0
    Insect%yc = 0.6
  case ('flapper_dickinsonII')
    !********************************************
    ! Digitized from Dickinson et al 1999 Science, figure 1A, drawing
    ! of the mechanical robot
    !********************************************
    Insect%nfft_wings = 20
    Insect%a0_wings = 0.6442788
    Insect%ai_wings(1:Insect%nfft_wings) = &
     (/0.0482978,-0.1208378,0.0061008,0.0356718,-0.0148328,-0.0109958,&
    0.0110268,0.0018538,-0.0061998,0.0015458,0.0025508,-0.0017538,&
    -0.0002578,0.0015018,-0.0003158,-0.0006048,0.0007168,-0.0001568,&
    -0.0005018,0.0004118/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/-0.0521708,0.0051828,0.0369428,-0.0002868,-0.0177448,0.0023218,&
    0.0081378,-0.0036288,-0.0038168,0.0031348,0.0011858,-0.0023828,&
    -0.0001638,0.0016098,-0.0004768,-0.0007188,0.0007228,0.0002278,&
    -0.0005798,0.0001228/)
    Insect%yc = 0.5282438
    Insect%xc = -0.1184548

    ! overwrite this because otherwise the wing is not entirely in the
    ! bounding box. (I actually do not recall why the L_span and L_chord are
    ! in the params file at all...)
    Insect%L_span = 1.2
  case ('robofly_dickinson')
    !********************************************
    ! Digitized from the hand drawn figure M. Dickinson sent via email, which
    ! contained the exact location of the pivot point. He also sent a CAD drawing
    ! which looks slightly different, and had no pivot point marked.
    !********************************************
    Insect%nfft_wings = 28
    Insect%a0_wings = 0.5313628
    Insect%ai_wings(1:Insect%nfft_wings) = &
    (/-0.0245658,-0.0842918,0.0218028,0.0105418,-0.0095288,0.0012928,&
    0.0021928,0.0000328,-0.0007648,-0.0015808,0.0013808,0.0013068,&
    -0.0010748,0.0002408,-0.0000378,-0.0010888,0.0008248,0.0004708,&
    -0.0003988,0.0002658,-0.0003178,-0.0004218,0.0002768,0.0000818,&
    0.0000318,0.0001228,-0.0001918,-0.0000558/)
    Insect%bi_wings(1:Insect%nfft_wings) = &
     (/-0.0905448,0.0278058,0.0392558,-0.0125248,-0.0159598,0.0048268,&
    0.0038898,-0.0028828,0.0012618,0.0012998,-0.0019058,0.0003118,&
    0.0003198,-0.0004298,0.0006388,-0.0000648,-0.0002308,0.0002518,&
    -0.0003948,0.0000928,0.0004478,-0.0003078,-0.0000888,0.0001638,&
    -0.0002348,0.0001398,0.0001398,-0.0002358/)
    Insect%yc = 0.4645238
    Insect%xc = -0.0716018

    ! overwrite this because otherwise the wing is not entirely in the
    ! bounding box. (I actually do not recall why the L_span and L_chord are
    ! in the params file at all...)
    Insect%L_span = 1.05
  case ('hawkmoth1')
    ! this wingshape is digitized from figure 1 from Kim et al. "Hovering and forward flight of the hawkmoth
    ! M. sexta: trim search and 6DOF dynamic stability characterization (Bioinspir Biomim. 10 (2015) 056012)"
    ! its area is about 0.30, which is lower than most references found. in their paper, they state A=0.3788,
    ! which raises question if fig 1 is to-scale drawing or sketch.
    Insect%nfft_wings = 28
    Insect%a0_wings = 0.5860758
    Insect%ai_wings(1:Insect%nfft_wings) = (/0.0219308,-0.1252418,0.0154668,0.0356038,-0.0203008,-0.0061968,&
    0.0178288,0.0002728,-0.0089908,0.0022758,0.0022948,-0.0046148,&
    -0.0008808,0.0032598,-0.0003708,-0.0019528,0.0006858,0.0008268,&
    -0.0008358,-0.0000718,0.0008938,0.0000348,-0.0004598,0.0004428,&
    0.0003158,-0.0003108,-0.0000658,0.0002798/)
    Insect%bi_wings(1:Insect%nfft_wings) = (/0.0062418,0.0452798,0.0303808,-0.0184998,-0.0179088,0.0068018,&
    0.0030268,-0.0058108,0.0017748,0.0040588,-0.0033678,-0.0030638,&
    0.0021898,0.0007208,-0.0015538,0.0007128,0.0016948,-0.0003828,&
    -0.0005898,0.0006388,0.0002888,-0.0005258,0.0000808,0.0002248,&
    -0.0004308,-0.0002758,0.0002298,-0.0000548/)
    Insect%yc = 0.4171918
    Insect%xc = -0.0395258
  case ('hawkmoth2')
    ! this wingshape is digitized from https://en.wikipedia.org/wiki/Manduca_sexta#/media/File:Manduca_sexta_female_sjh.JPG
    ! it has a greater aerea (A=0.40), but the original image is tricky since it is rotated. we therefore used a bit of modeling
    ! for this wing shape.
    Insect%nfft_wings = 18
    Insect%a0_wings = 0.6617728
    Insect%ai_wings(1:Insect%nfft_wings) = (/-0.0837648,-0.0802108,0.0703808,0.0069808,-0.0183478,0.0156518,&
    -0.0000308,-0.0153718,-0.0011538,0.0032378,-0.0005008,0.0020798,&
    0.0019888,-0.0009568,-0.0016378,-0.0010208,0.0005658,0.0009028/)
    Insect%bi_wings(1:Insect%nfft_wings) = (/0.0086968,0.0763208,0.0216658,-0.0322558,-0.0125988,0.0042128,&
    -0.0066278,-0.0040288,0.0093858,0.0045358,-0.0043238,0.0006298,&
    0.0010848,-0.0028958,0.0007268,0.0022578,-0.0013068,-0.0003538/)
    Insect%yc = 0.3946798
    Insect%xc = -0.2157968

  case default
    write (*,*) "Insect module: trying to set up fourier descriptors for wing&
                & shape but the type Insect%WingShape is unknown! :: "// Insect%WingShape
    call abort()
  end select

end subroutine Setup_Wing_Fourier_coefficients
