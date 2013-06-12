subroutine init_fields (n1, time,it, dt0, dt1, uk, work_nlk, vort, workvis)
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  integer, intent (inout) :: n1,it
  real (kind=pr), intent (inout) :: time, dt1, dt0
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (out) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent (out) :: work_nlk
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (inout) :: vort
  real (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent(inout) :: workvis
  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, r, a, gamma0, x00,r00, omega


  !-- Assign zero values
  time = 0.0d0
  dt1 = 0.0d0
  uk = dcmplx(0.0d0,0.0d0)
  work_nlk = dcmplx(0.0d0,0.0d0)
  workvis = 0.0
  it = 0  
  vort = 0.0d0


  if (inicond == 1) then
     if (mpirank==0) write (*,*) "*** inicond: vortex ring initial condition"
     !--------------------------------------------------
     ! Vortex ring
     !--------------------------------------------------
     r00 = yl/8.d0
     a   = 0.4131d0 * r00 
     a   = 0.82d0 * r00 
     gamma0 = 12.0d0
     x00 = 0.5d0 * xl


     ! define vorticity in phy space
     do iz = ra(3), rb(3)
        z = dble(iz) * zl / dble(nz)
        do iy = ra(2), rb(2)
           y = dble(iy) * yl / dble(ny)
           do ix = ra(1), rb(1)
              x = dble(ix) * xl / dble(nx)
              r = dsqrt( (y-0.5d0*yl)**2 + (z-0.5d0*zl)**2 )
              omega = (gamma0 / (pi*a**2)) * dexp(-( (x-x00)**2 + (r-r00)**2 )/(a**2) )

              if ( dabs(r)> 1.0d-12) then
                 vort (ix,iy,iz,2) = -omega * ( (z-0.5d0*zl)/r) ! sin(theta)
                 vort (ix,iy,iz,3) =  omega * ( (y-0.5d0*yl)/r) ! cos(theta)
              else
                 vort (ix,iy,iz,2) = 0.d0
                 vort (ix,iy,iz,3) = 0.d0
              endif

           end do
        end do
     end do

     call Vorticity2Velocity (uk, work_nlk(:,:,:,:,0), vort)


  elseif (inicond == 2) then
     if (mpirank==0) write (*,*) "*** inicond: turbulence (random vorticity) initial condition"
     !--------------------------------------------------
     ! random vorticity
     !--------------------------------------------------
     call random_seed
     do iz = ra(3), rb(3)
        do iy = ra(2), rb(2)
           do ix = ra(1), rb(1)
              call RANDOM_NUMBER(r)
              vort (ix,iy,iz,1) = 50.d0*(2.0d0*r - 1.d0) !* (1.d0-eps*mask(ix,iy,iz))
              call RANDOM_NUMBER(r)
              vort (ix,iy,iz,2) = 50.d0*(2.0d0*r - 1.d0) !* (1.d0-eps*mask(ix,iy,iz))
              call RANDOM_NUMBER(r)!
              vort (ix,iy,iz,3) = 50.d0*(2.0d0*r - 1.d0) !* (1.d0-eps*mask(ix,iy,iz))
           end do
        end do
     end do

     call Vorticity2Velocity (uk, work_nlk(:,:,:,:,0), vort)

  elseif (inicond == 30) then
     do iz=ra(3),rb(3)
        z=zl*(dble(iz)/dble(nz) -0.5d0)
        do iy=ra(2),rb(2)
           y=yl*(dble(iy)/dble(ny) -0.5d0)
           do ix=ra(1),rb(1)
              x=xl*(dble(ix)/dble(nx) -0.5d0)
              vort(ix,iy,iz,1)=-2.d0*dsin(y) ! vort is u in physical space
              vort(ix,iy,iz,2)=2.d0*dsin(x)
              vort(ix,iy,iz,3)=0.d0
           enddo
        enddo
     enddo
     call coftxyz(vort(:,:,:,1),uk(:,:,:,1))
     call coftxyz(vort(:,:,:,2),uk(:,:,:,2))
     call coftxyz(vort(:,:,:,3),uk(:,:,:,3))

  elseif (inicond == 0) then
     if (mpirank==0) write (*,*) "*** inicond: mean flow"
     !--------------------------------------------------
     ! mean flow only
     !--------------------------------------------------
     uk = dcmplx(0.0d0,0.0d0)
     if ( iMeanFlow == 1) then
        if ( ( ca(1) == 0 ) .and. ( ca(2) == 0 ) .and. ( ca(3) == 0 ) ) then
           uk(0, 0, 0,1) = Ux + Ax * time;                
           uk(0, 0, 0,2) = Uy + Ay * time;               
           uk(0, 0, 0,3) = Uz + Az * time;              
        endif
     endif
  elseif (inicond == 99) then 
     !--------------------------------------------------
     ! read from backup
     !--------------------------------------------------  
     call Read_Runtime_Backup ( time, dt0, dt1, n1, it, uk, work_nlk, workvis)
  elseif (inicond == 3) then 
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     !--------------------------------------------------
     ! fluid at rest
     !--------------------------------------------------  
     uk = dcmplx(0.0d0,0.0d0)        
  else
     if (mpirank == 0) then
        write (*,'(A)') '??? ERROR: Invalid initial condition'
     endif
     stop
  endif

end subroutine init_fields







subroutine Vorticity2Velocity (uk, work, vort)
  !----------------------------------------------------------------
  ! computes the divergence-free velocity in Fourier space u
  ! given vort in phys. space. 
  ! work is a work array
  !----------------------------------------------------------------
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (out) :: uk
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent (inout) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent (in) :: vort
  integer :: ix, iy, iz
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2

  !-------------------------------------------------
  ! compute vort in Fourier space
  !-------------------------------------------------
  call coftxyz (vort(:,:,:,1), work(:,:,:,1))
  call coftxyz (vort(:,:,:,2), work(:,:,:,2))
  call coftxyz (vort(:,:,:,3), work(:,:,:,3)) 

  !-------------------------------------------------
  ! compute streamfunction in Fourier space
  ! work(:,:,:,1:3, 1) will contain the three components of streamfunction
  !------------------------------------------------- 

  do iy = ca(3), cb(3)    ! ky : 0..ny/2-1 ,then, -ny/2..-1     
     ky = scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
     ky2 = ky*ky
     do ix = ca(2), cb(2)  ! kx : 0..nx/2
        kx = scalex*dble(ix)                
        kx2 = kx*kx
        do iz = ca(1),cb(1)  ! kz : 0..nz/2-1 ,then, -nz/2..-1           
           kz      = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           kz2     = kz*kz
           k_abs_2 = kx2+ky2+kz2
           if (abs(k_abs_2) .ne. 0.0) then  
              work(iz,ix,iy,1) = -work(iz,ix,iy,1) / k_abs_2
              work(iz,ix,iy,2) = -work(iz,ix,iy,2) / k_abs_2
              work(iz,ix,iy,3) = -work(iz,ix,iy,3) / k_abs_2
           else
              work(iz,ix,iy,1) = dcmplx(0.d0,0.d0)
              work(iz,ix,iy,2) = dcmplx(0.d0,0.d0)
              work(iz,ix,iy,3) = dcmplx(0.d0,0.d0)
           endif
        enddo
     enddo
  enddo

  !-----------------------------------------------
  !-- compute velocity as curl of streamfunction
  !-----------------------------------------------
  do iy = ca(3), cb(3) ! ky : 0..ny/2-1 ,then, -ny/2..-1     
     ky = scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
     do ix = ca(2), cb(2) ! kx : 0..nx/2
        kx = scalex*dble(ix)
        do iz = ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1
           kz = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           uk(iz,ix,iy,1) = dcmplx(0d0,1d0)*( &
                ky*work(iz,ix,iy,3) - kz*work(iz,ix,iy,2) )
           uk(iz,ix,iy,2) = dcmplx(0d0,1d0)*( &
                kz*work(iz,ix,iy,1) - kx*work(iz,ix,iy,3) )
           uk(iz,ix,iy,3) = dcmplx(0d0,1d0)*( &
                kx*work(iz,ix,iy,2) - ky*work(iz,ix,iy,1) )
        enddo
     enddo
  enddo


end subroutine Vorticity2Velocity
