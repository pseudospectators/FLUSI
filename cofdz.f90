subroutine cofdz (fk)
!---------------------------------------------------------------
!     calculation of d/dz in the fourier-space ==> *ik
!     for the third index
!     scaling included
!     P3DFFT version
!     with STRIDE-1 ordering (X,Y,Z)_physical->(Z,X,Y)_Fourier
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer iz, ix, iy
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: fk
  real (kind=pr) :: scale1

  do iy = ca(3), cb(3)
     ! modulo(...)   - z-wavenumber: first positive, 0..nz/2-1; then negative, -nz/2..-1 
     do ix = ca(2), cb(2) 
        fk(:, ix, iy) = dcmplx(0d0,1d0) * fk(:, ix, iy) * scalez &
                     * dble( modulo( (/(iz+nz/2, iz=ca(1),cb(1))/), nz ) - nz/2 )
     end do
  end do

end subroutine cofdz
