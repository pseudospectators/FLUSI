subroutine cal_vis (dt, vis)
  !---------------------------------------------------------------
  !  Calculate visfusive term for time advancement
  !  exp (-nu*k^2*dt)
  !  It is real valued, 
  !  its global size is 0:nz-1, 0:nx/2, 0:ny-1
  !  In newer versions, this is computed only if the time step
  !  changes. Is does so rarely, because we round it to one digit.
  !  also, dealiasing is done here (multiply aliased avenumbers by zero)
  !---------------------------------------------------------------
  use share_vars
  use mpi_header
  implicit none
  integer :: ix, iy, iz
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),intent(out)::vis
  real(kind=pr),intent(in) :: dt
  real(kind=pr) :: kx, ky, kz, t1, kx2, ky2, kz2, kx_trunc, ky_trunc, kz_trunc
  t1 = MPI_wtime()
  
  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)  
  
  do iy = ca(3), cb(3)
     ! ky - y-wavenumber: 0..ny/2-1 ,then, -ny/2..-1
     ky  = (scaley*dble(modulo(iy+ny/2,ny)-ny/2))**2
     ky2 = dble(modulo(iy+ny/2,ny)-ny/2)/ky_trunc
     ky2 = ky2*ky2
     
     do ix = ca(2), cb(2)
        ! kx - x-wavenumber: 0..nx/2
        kx  = (scalex*dble(ix))**2
        kx2 = dble(ix)/kx_trunc
	kx2 = kx2*kx2
	
        do iz = ca(1), cb(1)
           kz  = (scalez*dble(modulo(iz+nz/2,nz)-nz/2))**2
           kz2 = dble(modulo(iz+nz/2,nz)-nz/2)/kz_trunc
	   kz2 = kz2*kz2
	   
	   ! we can do dealiasing here, this is most efficient
	   if ((kx2 + ky2 + kz2  .ge. 1.d0).and.(iDealias==1)) then
              vis(iz,ix,iy) = dcmplx(0.d0,0.d0)
	   else
	      vis(iz,ix,iy) = dexp(-dt*nu*(kx +ky +kz))
           endif        
           
        enddo
     enddo
  enddo

  time_vis = time_vis + MPI_wtime() - t1
 
end subroutine cal_vis