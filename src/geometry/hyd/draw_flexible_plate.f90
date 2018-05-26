!-------------------------------------------------------------------------------
!-- Draws a flexible plate, where the deflection line is solved using the beam
!-- solver.
!-------------------------------------------------------------------------------
subroutine Draw_flexible_plate (time, beam)
  use mpi
  use vars
  use penalization ! mask array etc
  !-- use global variables for solid solver, as well as routines
  use solid_model
  implicit none

  real(kind=pr), intent(in) :: time
  type(Solid), intent(inout) :: beam
  real(kind=pr) :: psi,gamma,tmp,tmp2,psi_dt,beta_dt,gamma_dt
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate,u_tmp,rot_body_b,v_tmp,v0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate
  !-- for the triangles:
  real(kind=pr) :: a,b,c,alpha,beta,h,safety, s,s1,s2, ux,uy, R ,dmax, love
  real(kind=pr) :: c1,c2,u1,u2,uu
  !-- for leading edge state:
  real(kind=pr) :: alpha_t, alpha_tt, smoothing
  real(kind=pr), dimension(1:6) :: LeadingEdge
  integer :: ix,iy,iz,is

  ! thickness of the smoothing layer and safety distance for mask generation
  if (nx==1) then
    smoothing = N_smooth*max(dy,dz)
    safety = 2.d0*N_smooth*max(dy,dz)
  else
    smoothing = N_smooth*max(dx,dy,dz)
    safety = 2.d0*N_smooth*max(dx,dy,dz)
  endif


  !-- get relative coordinate system
  call plate_coordinate_system( time,x0_plate,v0_plate,psi,beta,&
       gamma,psi_dt,beta_dt,gamma_dt,M_plate)

  !-- fetch leading edge motion state (the cylinder may rotate)
  call mouvement ( time, alpha, alpha_t, alpha_tt, LeadingEdge, beam)

  ! angular velocity of moving relative frame
  ! ISSUE #9: https://github.com/pseudospectators/FLUSI/issues/9
  ! FIXME: the following line IS WRONG and needs to be corrected
  rot_body_b = (/psi_dt, beta_dt, gamma_dt/)

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0



  ! largest admissible sides (a,b) for triangles
  dmax = dsqrt( (ds)**2 + (t_beam+safety)**2 )

  !-- For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)
      !-- global coordinates
      x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz/)
      !-- in the plate system
      x_plate = matmul( M_plate, x-x0_plate )

      !-- check z-size in plate coordinate sytem (thus spanwise)
      if ((x_plate(3)>=-(0.5*L_span+t_beam+safety)).and.(x_plate(3)<=(0.5*L_span+t_beam+safety))) then
      !-- check x-size in plate system (so length, which is unity)
      if ((x_plate(1)>=-(0.0+t_beam+safety)).and.(x_plate(1)<=(1.0+t_beam+safety))) then
      !-- the beam bends in the x-y plane, so its maximum y-coordinate is -+1
      !-- if its handing straight down or up
      if ((x_plate(2)>=-(1.0+t_beam+safety)).and.(x_plate(2)<=(1.0+t_beam+safety))) then


      !-- initialize mask (locally) as very far away and velocity as zero
      mask(ix,iy,iz) = 5000.d0
      ux = 0.d0
      uy = 0.d0
      tmp2=0.d0

      !-------------------------------------------------------------------------
      !-- loop over points on the beam
      !-------------------------------------------------------------------------
      do is = 0,ns-2
        !-- a,b,c: sides of triangle
        a = dsqrt( (x_plate(1)-beam%x(is))**2 + (x_plate(2)-beam%y(is))**2 )
        b = dsqrt( (x_plate(1)-beam%x(is+1))**2 + (x_plate(2)-beam%y(is+1))**2 )

        ! at least one marker has to be close enough, otherwise we skip.
        if ((a<=dmax).or.(b<=dmax)) then
          !-- c is the distance between two markers, thus ds
          c = dsqrt( (beam%x(is)-beam%x(is+1))**2 + (beam%y(is)-beam%y(is+1))**2 )
          !-- height of the triangle: this is what we were looking for! (distance)
          h = dsqrt( max((a**2 - ((a**2+c**2-b**2)/(2.d0*c))**2),0.d0) )
          ! two parts on base side (i.e. on the centerline)
          c1 = dsqrt( a**2-h**2 )
          c2 = dsqrt( b**2-h**2 )

          !-- are both sides a,b small enough for the point to lie within the beam?
          if ((a<=dmax).and.(b<=dmax).and.(c1<c).and.(c2<c).and.(c1>0.d0).and.(c2>0.d0)) then
            !-----------------------------------------------------------------
            !-- we're in the area of the line segment, where the
            !-- height of the triangle defines the distance function
            !-----------------------------------------------------------------
            if ( h <= mask(ix,iy,iz) ) then
              mask(ix,iy,iz) = h

              !-- linear interpolation of velocity (along the element)
              s1 = dble(is)*ds
              s2 = dble(is+1)*ds
              s = s1 + c1
              ux = beam%vx(is) + ((s-s1)/(s2-s1))*(beam%vx(is+1)-beam%vx(is))
              uy = beam%vy(is) + ((s-s1)/(s2-s1))*(beam%vy(is+1)-beam%vy(is))

            endif

          else
            !-----------------------------------------------------------------
            !-- not on the segment, but maybe in a hinge?
            !-- which is the closer hinge?
            !-----------------------------------------------------------------
            if (a<b) then
              !-- we're in the region of the LEFT hinge
              if ( a <= mask(ix,iy,iz) ) then
                mask(ix,iy,iz) = a
                s = dble(is)*ds
                ! if this is the first point on the beam, assign an s coordinate
                ! that may be smaller than 0.0. we need the s-value for non-rectangular
                ! shapes (value of "tmp2")
                if (is==0) then
                  s = x_plate(1)*dcos(alpha) + x_plate(2)*dsin(alpha)
                endif
                ux = beam%vx(is)
                uy = beam%vy(is)
              endif
            elseif (b<a) then
              !-- we're in the region of the RIGHt hinge
              if ( b <= mask(ix,iy,iz) ) then
                mask(ix,iy,iz) = b
                s = dble(is+1)*ds
                ! if this is the last point on the beam, assign an "s" coordinate
                ! that may be greater than 1.0 we need the s-value for non-rectangular
                ! shapes (value of "tmp2")
                if ( is+1 == ns-1) then
                  !-- use the vector u that connects the last 2 points on the beam
                  uu = dsqrt( (beam%x(ns-1)-beam%x(ns-2))**2 +(beam%y(ns-1)-beam%y(ns-2))**2 )
                  u1 = (beam%x(ns-1)-beam%x(ns-2)) / uu
                  u2 = (beam%y(ns-1)-beam%y(ns-2)) / uu
                  s = 1.d0+((x_plate(1)-beam%x(is+1))*u1 +(x_plate(2)-beam%y(is+1))*u2)
                endif
                ux = beam%vx(is+1)
                uy = beam%vy(is+1)
              endif
            endif
          endif
        endif
      enddo

      ! mask is now the distance function from the centerline
      ! convert distance function to mask function
      call smoothstep( tmp, mask(ix,iy,iz), t_beam, smoothing )

      !-------------------------------------------------------------------------
      !-- make plate finite in z-direction, possibly with non-rectangular shapes
      !-------------------------------------------------------------------------
      if (infinite=="yes" .or. nx==1) then
        ! plate is very large - span is infinite, or flow is two-dimensional
        ! in these cases, the rigid direction is constant (spanwise) and we do not
        ! need to make the plate finite
        tmp2 = 1.d0
      else
        ! plate is finite, we have to do something about the finite length
        if (x_plate(3)<0.d0) then
          call smoothstep( tmp2, -x_plate(3), z_bottom(s), smoothing )
        else
          call smoothstep( tmp2,  x_plate(3), z_top(s), smoothing )
        endif
      endif

      !-- final value of mask function at this point
      mask(ix,iy,iz) = tmp*tmp2

      !-- assign mask color
      if (mask(ix,iy,iz) > 0.d0) mask_color(ix,iy,iz)=1

      !-- this is the velocity in the relative system
      u_tmp = (/ux,uy,0.d0/)
      !-- add solid body rotation to the velocity field of the beam
      v_tmp(1) = rot_body_b(2)*x_plate(3)-rot_body_b(3)*x_plate(2)
      v_tmp(2) = rot_body_b(3)*x_plate(1)-rot_body_b(1)*x_plate(3)
      v_tmp(3) = rot_body_b(1)*x_plate(2)-rot_body_b(2)*x_plate(1)

      ! FIXME: make sure that upon correction above, you're sure that v_tmp
      ! at this point is in BODY reference frame

      ! up to now, all velocities are in the relative system, now bring them back
      ! note v0_plate is already described in the global system
      us(ix,iy,iz,1:3) = v0_plate + matmul(transpose(M_plate),u_tmp+v_tmp)

      endif !-- end of bounding box checks
      endif
      endif
    enddo
   enddo
  enddo

  !-----------------------------------------------------------------------------
  ! Add cylinder add leading edge, if desired. Required for Turek's validation
  ! test case, and other cases. Note the angle ALPHA, which describes the leading
  ! edge angle WITHIN THE RELATIVE SYSTEM
  !-----------------------------------------------------------------------------
  if (has_cylinder=="yes") then


    !-- For all grid points of this subdomain
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          !-- global coordinates
          x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz/)
          !-- in the plate system
          x_plate = matmul( M_plate, x-x0_plate )
          !-- move the cylinder (in the direction normal to alpha, along L.E.)
          x_plate(1) = x_plate(1) + R_cylinder*dcos(alpha)
          x_plate(2) = x_plate(2) + R_cylinder*dsin(alpha)

          !-- check z-size in plate coordinate sytem (thus spanwise)
          if ((x_plate(3)>=-(0.5*L_span+safety)).and.(x_plate(3)<=(0.5*L_span+safety))) then
          !-- check if inside cylinder
          if ((x_plate(1)>=-(R_cylinder+safety)).and.(x_plate(1)<=(R_cylinder+safety))) then
          if ((x_plate(2)>=-(R_cylinder+safety)).and.(x_plate(2)<=(R_cylinder+safety))) then
            ! add smoothed cylinder
            R = dsqrt( x_plate(1)**2 + x_plate(2)**2 )
            call smoothstep ( tmp, R, R_cylinder, smoothing )

            !-- make plate finite in z-direction
            if ((nx>1).and.(infinite/="yes")) then
              call smoothstep( tmp2, abs(x_plate(3)), 0.5*L_span, smoothing )
            else
              !-- 2D runs have infinite span
              tmp2=1.d0
            endif

            if (infinite=="yes") tmp2 = 1.d0

            !-- velocity field of the cylinder (in relative system)
            v_tmp = cross( x_plate, (/0.d0,0.d0,-alpha_t/) )
            u_tmp = matmul(transpose(M_plate),v_tmp)

            !-- override mask if old value is smaller
            if (mask(ix,iy,iz)<tmp*tmp2) then
              mask(ix,iy,iz) = tmp*tmp2
              mask_color(ix,iy,iz) = 1
              us(ix,iy,iz,:) = u_tmp
            endif
          endif
          endif
          endif
        enddo
      enddo
    enddo
  endif
end subroutine Draw_flexible_plate
