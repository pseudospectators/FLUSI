subroutine poisson (fk)
!---------------------------------------------------------------
! Calculate solution to the Poisson equation fk_in = -\nabla^2 fk_out 
! in Fourier space
!---------------------------------------------------------------
  use share_vars
  implicit none

  integer :: kx, ky, kz, iy, iz
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: fk

  do iy = ca(3), cb(3)
     ! ky   - y-wavenumber: first positive, 0..ny/2-1; then negative, -ny/2..-1
     ky = modulo( iy+ny/2, ny ) - ny/2

     do kx = ca(2), cb(2)
        ! kx   - x-wavenumber: only positive, 0..nx/2
        
        do iz = ca(1),cb(1)
           ! kz   - z-wavenumber: first positive, 0..nz/2-1; then negative, -nz/2..-1
           kz = modulo( iz+nz/2, nz ) - nz/2

           if ( (kx == 0) .and. (ky == 0) .and. (kz == 0) ) then
              fk(iz, kx, iy) = 0.0d0
           else
              fk(iz, kx, iy) = fk(iz, kx, iy) / &
                 ( dble(kx)**2 *scalex**2 + dble(ky)**2 *scaley**2 + dble(kz)**2 *scalez**2 )
           end if
        enddo
     end do
  end do

end subroutine poisson
