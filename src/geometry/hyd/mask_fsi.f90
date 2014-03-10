! FSI wrapper for different (possibly time-dependend) mask functions 
subroutine create_mask_fsi(time)
  use mpi
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr) :: t1, eps_inv

  t1 = MPI_wtime() 
  
  !-------------------------------------------------------------
  ! create obstacle mask
  !-------------------------------------------------------------  
  ! do not create any mask when not using penalization
  if (iPenalization==1) then
    ! reset solid velocity
    if (iMoving==1) us = 0.d0
    ! Actual mask functions:
    select case (iMask)
    case ("sphere","Sphere")    
      call Draw_Sphere()
    case ("Flapper")    
      call Flapper (time)    
    case ("Insect")
      call Draw_Insect (time)
    case ("plate","Plate")
      call Draw_Plate (time) ! 2d plate, etc (Dmitry, 25 Oct 2013)
    case default    
      if (mpirank == 0) then
          write (*,*) "iMask="//iMask//" not properly set; stopping."
          stop
      endif
    end select
  endif

  !-------------------------------------------------------------
  ! cavity mask
  !-------------------------------------------------------------  
  ! if desired, add cavity mask surrounding the domain
  if ((iCavity=="yes").and.(iPenalization==1)) then
    call Add_Cavity ()
  endif
  
  
  !------------------------------------------------------------      
  ! For forces on wings/body
  !------------------------------------------------------------
  if (iMask=='Insect') then 
    eps_inv = 1.d0/eps
    maskpart = maskpart*eps_inv
  endif
  
  ! -- for global timing.
  time_mask = time_mask + MPI_wtime() - t1
end subroutine create_mask_fsi



subroutine update_us_fsi(ub)
  use fsi_vars
  implicit none
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  ! this is a stub in the FSI case since we always create the mask and the solid
  ! velocity field simultaneously
end subroutine update_us_fsi


subroutine Flapper (time)
  use mpi
  use fsi_vars
  implicit none

  real(kind=pr), intent(in) :: time
  integer :: iy, iz
  real (kind=pr) :: R, alpha_t, un, alpha_max
  real (kind=pr) :: y, z, ys,zs, alpha,L,H, tmp, N


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
              mask(:,iy,iz) = tmp

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
              mask(:,iy,iz) = tmp

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


! add a cavity mask around the domain
! note: if the obstacle extends in this zone (e.g. the flying insect touching
! the lower end of the domain), it will be overwritten.
subroutine Add_Cavity()
  use mpi
  use fsi_vars
  implicit none
  integer :: ix,iy,iz
  
  do ix = ra(1), rb(1)
    do iy = ra(2), rb(2)
       do iz = ra(3), rb(3)
       
          if ((ix<=cavity_size-1).or.(ix>=nx-1-cavity_size+1)) then            
            mask(ix,iy,iz) = 1.d0
            if (iMoving == 1) then
              us(ix,iy,iz,:) = 0.d0
            endif
          endif
          
          if ((iy<=cavity_size-1).or.(iy>=ny-1-cavity_size+1)) then       
            mask(ix,iy,iz) = 1.d0
            if (iMoving == 1) then
              us(ix,iy,iz,:) = 0.d0
            endif
          endif     
          
          if ((iz<=cavity_size-1).or.(iz>=nz-1-cavity_size+1)) then            
            mask(ix,iy,iz) = 1.d0
            if (iMoving == 1) then
              us(ix,iy,iz,:) = 0.d0
            endif
          endif
       enddo
    enddo
  enddo
 
end subroutine Add_Cavity
