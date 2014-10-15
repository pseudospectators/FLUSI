! this routine draws a rigid flapping plate, that is infinite in the x-direction
! and rotates around the x0-axis by the angle alpha.
subroutine Flapper (time, mask, mask_color, us)
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  integer :: iy,iz,ix
  real (kind=pr) :: R, alpha_t, un, alpha_max
  real (kind=pr) :: x,y,z,ys,zs, alpha,L,H,W, tmp1, N, tmp2
  real (kind=pr) :: safety

  ! motion protocoll (pitching around the x-axis with point y0,z0)
  alpha_max = deg2rad(14.d0)
  alpha   =           alpha_max*dsin(time*2.d0*pi)
  alpha_t = (2.d0*pi)*alpha_max*dcos(time*2.d0*pi)

  ! length of plate
  L = 1.d0
  ! half the plate thickness 
  H = 2.0d0*dx
  ! width of plate (Aspect ratio)
  W = 0.54d0
  ! smoothing coefficient
  N = 1.5d0
  ! safety
  safety = 2.d0*N*max(dx,dy,dz)+H

  do ix = ra(1), rb(1)
   do iy = ra(2), rb(2)
    do iz = ra(3), rb(3)
      x = dble(ix)*dx - x0
      y = dble(iy)*dy - y0
      z = dble(iz)*dz - z0

      ! transformed
      ys =  dcos(alpha)*y + dsin(alpha)*z
      zs = -dsin(alpha)*y + dcos(alpha)*z

      if ( (ys>=0.d0) .and. (ys<=L) )  then
        if ( (zs>=-H-safety) .and. (zs<=H+safety) )  then
          if ( (x>=-0.5d0*W-safety) .and. (x<=0.5d0*W+safety) )  then
            call SmoothStep (tmp1, abs(zs), H, N*max(dx,dy,dz))
            call SmoothStep (tmp2, abs(x), 0.5d0*W, N*max(dx,dy,dz))
            mask(ix,iy,iz) = tmp1*tmp2

            ! assign color "1" where >0 indicates something "useful"
            if (mask(ix,iy,iz) > 1.0d-12) then 
              mask_color(ix,iy,iz) = 1
              R = dsqrt( y**2 + z**2  )
              ! normal velocity
              un = R*alpha_t
              ! transform to fixed coordinate system
              us(:,iy,iz,2) = -dsin(alpha)*un
              us(:,iy,iz,3) = +dcos(alpha)*un
            endif
          endif
        endif
      endif
      
    enddo
   enddo
  enddo

  ! draw the endpoints
  do ix = ra(1), rb(1)
   do iy = ra(2), rb(2)
    do iz = ra(3), rb(3)
        x = dble(ix)*dx - x0
        y = dble(iy)*dy - y0
        z = dble(iz)*dz - z0

        ! leading edge (fixed endpoint)
        if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
           R = dsqrt( y**2 + z**2 )
           call SmoothStep (tmp1, R, H, N*max(dx,dy,dz))
           call SmoothStep (tmp2, dabs(x), 0.5d0*W, N*max(dx,dy,dz))
           
           ! overwrite if new value is larger
           if (mask(ix,iy,iz)<tmp1*tmp2) then
              mask(ix,iy,iz) = tmp1*tmp2
              ! assign color "1" where >0 indicates something "useful"
              if (mask(ix,iy,iz) > 1.0e-12) then
                mask_color(:,iy,iz) = 1
                ! velocity (remember: rotation around x0,y0,z0)
                y = dble(iy)*dy - y0
                z = dble(iz)*dz - z0
                R = dsqrt( y**2 + z**2  )
                un = R*alpha_t
                us(ix,iy,iz,2) = -dsin(alpha)*un
                us(ix,iy,iz,3) = +dcos(alpha)*un
              endif
           endif
        endif

        ! trailing edge (moving endpoint)
        x = dble(ix)*dx - x0
        y = dble(iy)*dy - (y0 + L*dcos(alpha))
        z = dble(iz)*dz - (z0 + L*dsin(alpha))
        
        if ( ((y>=-15.d0*dy).and.(y<=+15.d0*dy)) .and.((z>=-15.d0*dz).and.(z<=+15.d0*dz))  )then
           R = dsqrt( y**2 + z**2 )
           call SmoothStep (tmp1, R, H, N*max(dx,dy,dz))
           call SmoothStep (tmp2, dabs(x), 0.5d0*W, N*max(dx,dy,dz))
           ! overwrite if new value is larger
           if (mask(ix,iy,iz)<tmp1*tmp2) then
              mask(ix,iy,iz) = tmp1*tmp2
              if (mask(ix,iy,iz) > 1.0e-12) then 
                ! assign color "1" where >0 indicates something "useful"
                mask_color(ix,iy,iz) = 1
                ! velocity (remember: rotation around x0,y0,z0)
                y = dble(iy)*dy - y0
                z = dble(iz)*dz - z0
                R = dsqrt( y**2 + z**2  )
                un = R*alpha_t
                us(ix,iy,iz,2) = -dsin(alpha)*un
                us(ix,iy,iz,3) = +dcos(alpha)*un
              endif
           endif
        endif

     enddo
  enddo
  enddo
end subroutine Flapper


!-------------------------------------------------------------------------------
! cavity as used by romain for "open cavity" tests. the wall is from 
! -2 ... -1 and +1 ... +2
! xxxx------------xxxx
! constant in all other directions. this test was used for reproducing the 
! results of Benjamin's JCP for the Neumann BC
!-------------------------------------------------------------------------------
subroutine romain_open_cavity(mask, mask_color, us)
  use vars
  implicit none
  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x=dble(ix)*dx - 0.5*xl
        if (x<=-1.d0 .or. x>=+1.d0 ) then
          mask (ix, iy, iz) = 1.d0
          us (ix,iy,iz,1:3) = 0.d0
          ! assign color "1" where >0 indicates something "useful"
          mask_color(ix,iy,iz) = 1
        endif
      enddo
    enddo
  enddo
end subroutine romain_open_cavity


subroutine taylor_couette(mask, mask_color, us)
  use vars
  implicit none
  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, R,R1,R2,omega
  
  R1=0.5d0
  R2=1.0d0
  omega=1.25d0
  
  do iz=ra(3),rb(3)
    z = dble(iz)*dz - 0.5d0*zl
    do iy=ra(2),rb(2)
      y = dble(iy)*dy - 0.5d0*yl
      
      R=dsqrt(z*z+y*y)
      
      ! inner cylinder
      if ( R<=R1) then
        mask (ix, iy, iz) = 1.d0
        us (ix,iy,iz,1) = 0.d0
        us (ix,iy,iz,2) = +omega * z
        us (ix,iy,iz,3) = -omega * y 
        mask_color(ix,iy,iz) = 0
      endif
      
      ! outer cylinder
      if (R>=R2) then
        mask (ix, iy, iz) = 1.d0
        us (ix,iy,iz,1:3) = 0.d0
        mask_color(ix,iy,iz) = 0
      endif
    enddo
  enddo
  
end subroutine