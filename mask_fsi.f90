! FSI wrapper for different (possibly time-dependend) mask functions 
subroutine Create_Mask_fsi(time)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr) :: t1

  t1 = MPI_wtime() 
  
  ! do not create any mask when not using penalization
  if (iPenalization==1) then
    ! reset solid velocity
    if (iMoving==1) us = 0.d0
    ! Actual mask functions:
    select case (iMask)
    case ("sphere")    
      call Draw_Sphere()
    case ("Flapper")    
      call Flapper (time)    
    case ("Insect")
      call Draw_Insect (time)
    case default    
      if (mpirank == 0) then
          write (*,*) "iMask="//iMask//" not properly set; stopping."
          stop
      endif
    end select
  endif
  ! -- for global timing.
  time_mask = time_mask + MPI_wtime() - t1
end subroutine Create_Mask_fsi



subroutine update_us_fsi(ub)
  use fsi_vars
  implicit none
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  ! this is a stub in the FSI case since we always create the mask and the solid
  ! velocity field simultaneously
end subroutine update_us_fsi


subroutine Flapper (time)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  integer :: iy, iz
  real (kind=pr) :: R, alpha_t, un, alpha_max
  real (kind=pr) :: y, z, ys,zs, alpha,L,H, tmp, N

  real(kind=pr) :: eps_inv
  
  eps_inv = 1.d0 / eps


  alpha_max = 30.d0*pi/180.d0
  alpha   =           alpha_max*dsin(time*2.d0*pi)
  alpha_t = (2.d0*pi)*alpha_max*dcos(time*2.d0*pi)
  
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
           mask(:,iy,iz) = tmp

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
