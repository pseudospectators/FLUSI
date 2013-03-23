subroutine cal_vis (dt, vis)
  !---------------------------------------------------------------
  !  Calculate visfusive term for time advancement
  !  exp (-nu*k^2*dt)
  !  It is real valued, 
  !  its global size is 0:nz-1, 0:nx/2, 0:ny-1
  !---------------------------------------------------------------
  use share_vars
  use mpi_header
  implicit none
  integer :: ix, iy, iz
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),intent(out)::vis
  real(kind=pr),intent(in) :: dt
  real(kind=pr) :: kx, ky, kz, t1
  t1 = MPI_wtime()
  
  
  do iy = ca(3), cb(3)
     ! ky - y-wavenumber: 0..ny/2-1 ,then, -ny/2..-1
     ky=(scaley*dble(modulo(iy+ny/2,ny)-ny/2))**2
     do ix = ca(2), cb(2)
        ! kx - x-wavenumber: 0..nx/2
        kx=(scalex*dble(ix))**2
        do iz = ca(1), cb(1)
           kz=(scalez*dble(modulo(iz+nz/2,nz)-nz/2))**2
           vis(iz,ix,iy)=dexp(-dt*nu*(kx +ky +kz))
        enddo
     enddo
  enddo

  time_vis = time_vis + MPI_wtime() - t1
end subroutine cal_vis
