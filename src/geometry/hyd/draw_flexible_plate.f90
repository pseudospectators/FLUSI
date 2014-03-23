!-------------------------------------------------------------------------------
!-- Draws a flexible plate, where the deflection line is solved using the beam 
!-- solver. 
!-------------------------------------------------------------------------------
subroutine Draw_flexible_plate (time, beam)
  use mpi
  use fsi_vars
  !-- use global variables for solid solver, as well as routines
  use solid_model 
  implicit none

  real(kind=pr), intent(in) :: time
  type(Solid), intent(inout) :: beam
  real(kind=pr) :: psi,gamma,tmp,tmp2,psi_dt,beta_dt,gamma_dt
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate,u_tmp,rot_body,v_tmp,v0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate
  !-- for the triangles:
  real(kind=pr) :: a,b,c,alpha,beta,h, safety, s,s1,s2, ux,uy, R
  integer :: ix,iy,iz,is
  
  N_smooth = 3.d0  
  safety = 4.d0*N_smooth*max(dx,dy,dz)
    
  !-- get relative coordinate system
  call plate_coordinate_system( time,x0_plate,v0_plate,psi,beta,&
                               gamma,psi_dt,beta_dt,gamma_dt,M_plate)

  ! angular velocity of moving relative frame
  rot_body = (/psi_dt, beta_dt, gamma_dt/)

  mask = 0.d0
  us = 0.d0
  
  !-- For all grid points of this subdomain
  do iz = ra(3), rb(3)
   do iy = ra(2), rb(2)
    do ix = ra(1), rb(1)
      !-- global coordinates
      x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz/)
      !-- in the plate system
      x_plate = matmul( M_plate, x-x0_plate )
      
      
      !-- check z-size in plate coordinate sytem (thus spanwise)
      if ((x_plate(3)>=-(0.5*L_span+safety)).and.(x_plate(3)<=(0.5*L_span+safety))) then
      !-- check x-size in plate system (so length, which is unity)
      if ((x_plate(1)>=-(0.0+safety)).and.(x_plate(1)<=(1.0+safety))) then
      !-- the beam bends in the x-y plane, so its maximum y-coordinate is -+1
      !-- if its handing straight down or up
      if ((x_plate(2)>=-(1.0+safety)).and.(x_plate(2)<=(1.0+safety))) then
      
      
      !-- initialize mask (locally) as very far away and velocity as zero
      mask(ix,iy,iz) = 5000.d0
      ux = 0.d0
      uy = 0.d0
      
      !-- loop over points on the beam
      do is = 0,ns-2
          !-- a,b,c: sides of triangle
          a = dsqrt( (x_plate(1)-beam%x(is))**2 + (x_plate(2)-beam%y(is))**2 )
          !-- if the first lagrangian marker is too far away, the second
          !-- one is as well, and we can skip the whole point, greatly
          !-- reducing the computational complexity
          if ( a < t_beam+safety ) then
            b = dsqrt( (x_plate(1)-beam%x(is+1))**2 + (x_plate(2)-beam%y(is+1))**2 )
            !-- c is the distance between two markers, thus ds
            c = dsqrt( (beam%x(is)-beam%x(is+1))**2 + (beam%y(is)-beam%y(is+1))**2 )
            !-- angles in the triangle (Law of cosines)
            alpha = acos ( (b*b+c*c-a*a)/(2.d0*b*c) )
            beta  = acos ( (a*a+c*c-b*b)/(2.d0*a*c) )
            
            !-- height of the triangle: this is what we were looking for!
            h = dsin(alpha)*b
            
            
            if ((abs(alpha)<=pi/2.d0).and.(abs(beta)<=pi/2.d0)) then
                !-- we're in the area of the line segment, where the
                !-- height of the triangle defines the distance function
                if ( h<= mask(ix,iy,iz) ) then
                  mask(ix,iy,iz) = h
                
                  !-- linear interpolation of velocity (along the element)
                  s1 = dble(is)*ds
                  s2 = dble(is+1)*ds
                  s = s1 + dsqrt(b*b - h*h)
                  ux = beam%vx(is) + ((s-s1)/(s2-s1))*(beam%vx(is+1)-beam%vx(is))
                  uy = beam%vy(is) + ((s-s1)/(s2-s1))*(beam%vy(is+1)-beam%vy(is))
                endif
            else
                !-- we're in a hinge zone, where it is the closest
                !-- lagrangian marker that defines the distance (circular hinges)
                mask(ix,iy,iz) = min(mask(ix,iy,iz),min(a,b))
                !-- assign the velocity of the closer hinge (if any)
                if (a==mask(ix,iy,iz)) then
                  ux = beam%vx(is)
                  uy = beam%vy(is)
                elseif (b==mask(ix,iy,iz) ) then
                  ux = beam%vx(is+1)
                  uy = beam%vy(is+1)
                endif
                
            endif          
          endif
      enddo
  
      !-- mask is now the distance function from the centerline
      !-- convert distance function to mask function
      call smoothstep( tmp, mask(ix,iy,iz)-t_beam, 0.d0, N_smooth*max(dx,dy,dz) )
      !-- make plate finite in z-direction
      call smoothstep( tmp2, abs(x_plate(3)), 0.5*L_span, N_smooth*max(dx,dy,dz) )
      !-- final value
      mask(ix,iy,iz) = tmp*tmp2
      
      
      !-- this is the velocity in the relative system
      u_tmp = (/ux,uy,0.d0/)
      !-- add solid body rotation to the velocity field of the beam
      v_tmp(1) = rot_body(2)*x_plate(3)-rot_body(3)*x_plate(2)
      v_tmp(2) = rot_body(3)*x_plate(1)-rot_body(1)*x_plate(3)
      v_tmp(3) = rot_body(1)*x_plate(2)-rot_body(2)*x_plate(1)
      !-- up to now, all velocities are in the relative system, now bring them back
      !-- note v0_plate is already described in the ablosute system
      us(ix,iy,iz,1:3)=v0_plate + matmul(transpose(M_plate),u_tmp+v_tmp)
      
      endif !-- end of bounding box checks
      endif
      endif
    enddo
   enddo
  enddo
  
  
  !-----------------------------------------------------------------------------
  ! Add cylinder add leading edge, if desired. Required for Turek's validation
  ! test case.
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
          !-- move the cylinder
          x_plate(1) = x_plate(1)+R_cylinder
          
          !-- check z-size in plate coordinate sytem (thus spanwise)
          if ((x_plate(3)>=-(0.5*L_span+safety)).and.(x_plate(3)<=(0.5*L_span+safety))) then
          !-- check if inside cylinder 
          if ((x_plate(1)>=-(R_cylinder+safety)).and.(x_plate(1)<=(R_cylinder+safety))) then
          if ((x_plate(2)>=-(R_cylinder+safety)).and.(x_plate(2)<=(R_cylinder+safety))) then
            ! add smoothed cylinder
            R = dsqrt( x_plate(1)**2 + x_plate(2)**2 )
            call smoothstep ( tmp, R, R_cylinder, N_smooth*max(dx,dy,dz) )
            
            !-- make plate finite in z-direction
            call smoothstep( tmp2, abs(x_plate(3)), 0.5*L_span, N_smooth*max(dx,dy,dz) )
            
            !-- override mask if old value is smaller
            if (mask(ix,iy,iz)<=tmp*tmp2) then
              mask(ix,iy,iz) = tmp*tmp2
              us(ix,iy,iz,:) = 0.d0
            endif
          endif
          endif
          endif
        enddo
      enddo
    enddo
  endif
end subroutine Draw_flexible_plate


