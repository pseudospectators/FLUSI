! FSI wrapper for different (possibly time-dependend) mask functions 
subroutine Create_Mask_fsi(time)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr) :: t1

  t1 = MPI_wtime() 

  ! Actual mask functions:
  select case (iMask)
  case ("sphere")
     call Draw_Sphere()
  case ("Jerry")      
     call Draw_Jerry ( time )      
  case ("Flapper")    
     call Flapper ( time )    
  case default    
     if (mpirank == 0) then
        write (*,*) "iMask not properly set. Suicide!"
        stop
     endif
  end select

  ! -- for global timing.
  time_mask = time_mask + MPI_wtime() - t1
end subroutine Create_Mask_fsi


! Draws Jerry, the first insect with rigid wings.  Jerry is symmetric
! w.r.t. the x-axis, but since we do MPI we currently do not take
! advantage of this. (unlike the experimental matlab file) Jerry has
! two wings, a body and a head with eyes.
subroutine Draw_Jerry (time)
  use fsi_vars
  use mpi_header
  implicit none

  real(kind=pr), intent(in) :: time
  integer :: ix, iy, iz
  real (kind=pr) :: phi_max, alpha_max, phase,phi,phi_dt,alpha,alpha_dt, y0_wing, x0_wing, z0_wing
  real (kind=pr) :: x, y, z, xs, ys, zs, L_body, a_body, b_body, x0_body, y0_body, z0_body, x0_head,y0_head,z0_head
  real (kind=pr) :: R_head, R_eye, x0_eye, y0_eye, z0_eye, sin_45
  real (kind=pr) :: b_top, b_bot, xx,yy, zz, cos_a, sin_a, cos_p, sin_p,&
       L_chord, delta, N_smooth, R, R_limit, tmp, x_tmp, y_tmp, z_tmp, y_bot,y_top

  ! parameters for motion protocoll
  phi_max     = 60.d0*pi/180.d0  ! phi is up/down angle (flapping)
  alpha_max   = 45.d0*pi/180.d0  ! alpha is tethering
  phase       = 10.d0*pi/180.d0  ! phase shift between flapping and tethering
  
  phi         = phi_max*dcos(time)
  phi_dt      =-phi_max*dsin(time)
  alpha       = alpha_max*dsin(time+phase)
  alpha_dt    = alpha_max*dcos(time+phase)

  L_body = 1.d0  ! length of Jerrys body
  b_body = 0.1d0 ! smaller ellipse axis length
  a_body = L_body/2.d0 ! larger ellipse axis length
  L_chord= 1.d0  ! chord length of Jerry's wing

  ! center of Jerry's body
  x0_body = xl/2.d0
  y0_body = yl/2.d0
  z0_body = zl/2.d0

  ! head root point
  R_head = 1.25d0*b_body
  x0_head = x0_body
  y0_head = y0_body + 0.9d0*a_body + R_head
  z0_head = z0_body

  ! wing shape ( top and bottom ellipses )
  b_top = 0.1d0
  b_bot = 0.3d0

  ! save a bit computing time
  cos_a = dcos(alpha)
  sin_a = dsin(alpha)
  cos_p = dcos(phi)
  sin_p = dsin(phi)
  sin_45 = dsin(45.d0*pi/180.d0)
  
  ! width of smoothing zone
  N_smooth = 2.0
  delta = 2.0*N_smooth*max(dx,dy,dz)  ! delta is the safety distance to account for smoothing layer

  ! initialize fields
  mask = 0.d0
  us = 0.d0

  do iz = ra(3), rb(3)
     do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
           ! global coordinates
           x = dble(ix)*dx
           y = dble(iy)*dy
           z = dble(iz)*dz

           !!! First wing !!!
           ! wing root point
           y0_wing = y0_body + 0.50*a_body ! wing root point can be off the center of the body
           x0_wing = x0_body + dsqrt(b_body**2 *(1.d0- ((y0_wing-y0_body)/a_body)**2 ) ) ! get the body thickness at this point
           z0_wing = z0_body

           ! wing root system
           xx = x - x0_wing
           yy = y - y0_wing
           zz = z - z0_wing

           ! wing root system, rotated!
           xs = xx*cos_p - zz*sin_p
           ys = cos_a*yy + sin_a*(+xx*sin_p+zz*cos_p)
           zs =-sin_a*yy + cos_a*(+xx*sin_p+zz*cos_p)

           ! wing shape
           if ( ((xs-a_body)**2)/(a_body**2) <= 1.d0 ) then
              y_top = sqrt( (b_top**2)*(1.d0 - ((xs-a_body)**2)/(a_body**2)) )
              y_bot =-sqrt( (b_bot**2)*(1.d0 - ((xs-a_body)**2)/(a_body**2)) )
           else
              y_top = 0.d0
              y_bot = 0.d0
           endif
           !  y_top = b_top
           !  y_bot =-b_bot

           ! WING
           if ( ((xs>=0.d0-delta).and.(xs<=L_chord+delta)) &
                .and. (dabs(zs)<=3.d0*max(dx,dy,dz)+delta) &
                .and. ((ys>y_bot-delta).and.(ys<y_top+delta)) ) then !is it inside the rectangle?
              ! smooth length
              if (xs<0.d0) then  ! xs is chordlength coordinate  
                 call SmoothStep ( x_tmp, -xs, 0.d0, N_smooth*max(dx,dy,dz) )
              else
                 call SmoothStep ( x_tmp,  xs, L_chord, N_smooth*max(dx,dy,dz) )
              endif

              ! smooth height
              call SmoothStep ( z_tmp, dabs(zs), 3.d0*max(dx,dy,dz), N_smooth*max(dx,dy,dz) ) ! note we define 4px thickness of the wings here


              ! smooth shape
              if (ys<0.d0) then
                 call SmoothStep ( y_tmp, -ys, -y_bot, N_smooth*max(dx,dy,dz) )
              else  
                 call SmoothStep ( y_tmp,  ys,  y_top, N_smooth*max(dx,dy,dz) )
              endif

              mask(ix,iy,iz) = z_tmp*y_tmp*x_tmp ! note we draw the wings first, the mask is empty. no need to check if there is already something

              if (mask(ix,iy,iz) > 1.d-4) then
                 us(ix,iy,iz,1) =    phi_dt*zz
                 us(ix,iy,iz,2) = -alpha_dt*zz
                 us(ix,iy,iz,3) =  alpha_dt*yy-phi_dt*xx
              endif
           endif

           !!! Second wing !!!
           ! wing root point
           x0_wing = x0_body - dsqrt(b_body**2 *(1.d0- ((y0_wing-y0_body)/a_body)**2 ) ) ! get the body thickness at this point

           ! wing root system
           xx = x0_wing -x  ! note change in sign here

           ! wing root system, rotated!
           xs = xx*cos_p - zz*sin_p
           ys = cos_a*yy + sin_a*(+xx*sin_p+zz*cos_p)
           zs =-sin_a*yy + cos_a*(+xx*sin_p+zz*cos_p)

           ! wing shape
           if ( ((xs-a_body)**2)/(a_body**2) <= 1.d0 ) then
              y_top = sqrt( (b_top**2)*(1.d0 - ((xs-a_body)**2)/(a_body**2)) )
              y_bot =-sqrt( (b_bot**2)*(1.d0 - ((xs-a_body)**2)/(a_body**2)) )
           else
              y_top = 0.d0
              y_bot = 0.d0
           endif
           !  y_top = b_top
           !  y_bot =-b_bot 
           !--------------------------------------------
           ! WING
           !--------------------------------------------
           if ( ((xs>=0.d0-delta).and.(xs<=L_chord+delta)) &
                .and. (dabs(zs)<=3.d0*max(dx,dy,dz)+delta) &
                .and. ((ys>y_bot-delta).and.(ys<y_top+delta)) ) then !is it inside the rectangle?
              ! smooth length
              if (xs<0.d0) then  ! xs is chordlength coordinate  
                 call SmoothStep ( x_tmp, -xs, 0.d0, N_smooth*max(dx,dy,dz) )
              else
                 call SmoothStep ( x_tmp,  xs, L_chord, N_smooth*max(dx,dy,dz) )
              endif

              ! smooth height
              call SmoothStep ( z_tmp, dabs(zs), 3.d0*max(dx,dy,dz), N_smooth*max(dx,dy,dz) ) ! note we define 4px thickness of the wings here


              ! smooth shape
              if (ys<0.d0) then
                 call SmoothStep ( y_tmp, -ys, -y_bot, N_smooth*max(dx,dy,dz) )
              else  
                 call SmoothStep ( y_tmp,  ys,  y_top, N_smooth*max(dx,dy,dz) )
              endif
              mask(ix,iy,iz) = z_tmp*y_tmp*x_tmp;

              if (mask(ix,iy,iz) > 1.d-4) then
                 us(ix,iy,iz,1) =              -phi_dt*zz
                 us(ix,iy,iz,2) = -alpha_dt*zz
                 us(ix,iy,iz,3) =  alpha_dt*yy-phi_dt*xx
              endif
           endif

           ! Body (ellipsoid)
           if ( dabs(y-y0_body) < L_body/2.d0 + delta ) then ! are we possibly inside the body?
              ! radius
              R = dsqrt( (x-x0_body)**2 + (z-z0_body)**2 )     
              ! compute R as a function of the position (y-y0_body)
              if ( ((y-y0_body)/a_body)**2 <= 1.d0) then
                 R_limit = dsqrt(b_body**2 *(1.d0- ((y-y0_body)/a_body)**2 ) )
                 if ( R < R_limit + delta ) then
                    call SmoothStep ( tmp,R,R_limit,N_smooth*max(dx,dy,dz) )
                    if ( tmp > mask(ix,iy,iz) ) then
                       mask(ix,iy,iz) = tmp ! note now that the wings are drawn, be careful not to overwrite existing values.
                       us (ix,iz,iz,:) = 0.d0
                    endif
                 endif
              endif
           endif

           ! Head (sphere)
           if ( dabs(y-y0_head) < R_head+delta ) then
              R = dsqrt ( (x-x0_head)**2 + (z-z0_head)**2 + (y-y0_head)**2 )
              if ( R < R_head + delta ) then
                 call SmoothStep ( tmp,R,R_head,N_smooth*max(dx,dy,dz) )
                 if ( tmp > mask(ix,iy,iz) ) then
                    mask(ix,iy,iz) = tmp ! note now that the wings are drawn, be careful not to overwrite existing values.
                    us (ix,iz,iz,:) = 0.d0
                 endif
              endif
           endif

           !!! eye 1 (sphere)
           ! eye root point
           R_eye = 0.50d0*R_head;
           x0_eye = x0_head + sin_45*R_head*0.80d0 ! the 0.8 moves the eye into the head
           y0_eye = y0_head + sin_45*R_head*0.80d0
           z0_eye = z0_head + sin_45*R_head*0.80d0

           if ( dabs(y-y0_eye) < R_eye+delta ) then
              R = dsqrt ( (x-x0_eye)**2 + (z-z0_eye)**2 + (y-y0_eye)**2 )
              if ( R < R_eye + delta ) then
                 call SmoothStep ( tmp,R,R_eye,N_smooth*max(dx,dy,dz) )
                 if ( tmp > mask(ix,iy,iz) ) then
                    mask(ix,iy,iz) = tmp ! note now that the wings are drawn, be careful not to overwrite existing values.
                    us (ix,iz,iz,:) = 0.d0
                 endif
              endif
           endif

           !!! eye 2 (sphere)
           ! eye root point
           R_eye = 0.50d0*R_head;
           x0_eye = x0_head - sin_45*R_head*0.80d0 ! note sign change for symmetry
           if ( dabs(y-y0_eye) < R_eye+delta ) then
              R = dsqrt ( (x-x0_eye)**2 + (z-z0_eye)**2 + (y-y0_eye)**2 )
              if ( R < R_eye + delta ) then
                 call SmoothStep ( tmp,R,R_eye,N_smooth*max(dx,dy,dz) )
                 if ( tmp > mask(ix,iy,iz) ) then
                    mask(ix,iy,iz) = tmp ! note now that the wings are drawn, be careful not to overwrite existing values.
                    us (ix,iz,iz,:) = 0.d0
                 endif
              endif
           endif

        enddo
     enddo
  enddo
end subroutine Draw_Jerry


subroutine Flapper (time)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  integer :: iy, iz
  real (kind=pr) :: R, alpha_t, un, alpha_max
  real (kind=pr) :: y, z, ys,zs, alpha,L,H, tmp,N, eps_inv

  alpha_max = 30.d0*pi/180.d0
  alpha   =           alpha_max*dsin(time*2.d0*pi)
  alpha_t = (2.d0*pi)*alpha_max*dcos(time*2.d0*pi)
  eps_inv = 1.d0/eps
  
  y0 = 1.d0
  z0 = 1.d0

  L = 1.d0
  H = 0.0625

  ! initialize fields
  mask = 0.d0
  us = 0.d0

  N=3.0 ! smoothing coefficient

  do iz = ra(3), rb(3)
     do iy = ra(2), rb(2)
        y = dble(iy)*dy - y0
        z = dble(iz)*dz - z0

        ! transformed
        ys =  dcos(alpha)*y + dsin(alpha)*z
        zs = -dsin(alpha)*y + dcos(alpha)*z

        if ( (ys>=0.d0) .and. (ys<=L) )  then

           call SmoothStep (tmp, abs(zs), H, N*max(dx,dy,dz))
           mask(:,iy,iz) = tmp*eps_inv

           if (dabs(zs)<=2.d0*H) then
              R = dsqrt( y**2 + z**2  )
              un = R*alpha_t
              us(:,iy,iz,2) = -dsin(alpha)*un
              us(:,iy,iz,3) = +dcos(alpha)*un
           endif

        endif
     enddo
  enddo

  ! draw the endpoints
  do iz = ra(3), rb(3)
     do iy = ra(2), rb(2)
        y = dble(iy)*dy - y0
        z = dble(iz)*dz - z0

        ! origin
        if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
           R = dsqrt( y**2 + z**2 )
           call SmoothStep (tmp, R, H, N*max(dx,dy,dz))
           if (mask(1,iy,iz)<tmp) then
              mask(:,iy,iz) = tmp*eps_inv

              R = dsqrt( y**2 + z**2  )
              un = R*alpha_t
              us(:,iy,iz,2) = -dsin(alpha)*un
              us(:,iy,iz,3) = +dcos(alpha)*un
           endif
        endif

        y = dble(iy)*dy - (y0 + L*dcos(alpha))
        z = dble(iz)*dz - (z0 + L*dsin(alpha))

        if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
           R = dsqrt( y**2 + z**2 )
           call SmoothStep (tmp, R, H, N*max(dx,dy,dz))
           if (mask(1,iy,iz)<tmp) then
              mask(:,iy,iz) = tmp*eps_inv

              y = dble(iy)*dy - y0
              z = dble(iz)*dz - z0
              R = dsqrt( y**2 + z**2  )
              un = R*alpha_t
              us(:,iy,iz,2) = -dsin(alpha)*un
              us(:,iy,iz,3) = +dcos(alpha)*un
           endif
        endif

     enddo
  enddo
end subroutine Flapper
