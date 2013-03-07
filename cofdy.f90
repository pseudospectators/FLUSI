subroutine cofdy (fk)
!---------------------------------------------------------------
!     calculation of d/dy in the fourier-space ==> *ik
!     for the second index
!     scaling included
!     P3DFFT version
!     with STRIDE-1 ordering (X,Y,Z)_physical->(Z,X,Y)_Fourier
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer k, ix, iy
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: fk
  real (kind=pr) :: scale1

  do iy = ca(3), cb(3)
     ! k   - y-wavenumber: first positive, 0..ny/2-1; then negative, -ny/2..-1
     k = modulo( iy+ny/2, ny ) - ny/2
     do ix = ca(2), cb(2)
        fk(:, ix, iy) = dcmplx(0d0,1d0) * dble(k) * fk(:, ix, iy) * scaley
     end do
  end do

end subroutine cofdy
