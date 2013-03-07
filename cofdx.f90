! this routine is actually mood
subroutine cofdx (fk)
!---------------------------------------------------------------
!     calculation of d/dx in the fourier-space ==> *ik
!     for the first index
!     scaling included
!     P3DFFT version 
!     with STRIDE-1 ordering (X,Y,Z)_physical->(Z,X,Y)_Fourier
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer k, iy
  complex (kind=pr), dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)), intent (inout) :: fk      
  real (kind=pr) :: scale1

  do iy = ca(3), cb(3)
     do k = ca(2), cb(2)  
        ! k   - x-wavenumber: only positive, 0..nx/2
        fk(:, k, iy) = dcmplx(0d0,1d0) * dble(k) * fk(:, k, iy) * scalex
     end do
  end do

end subroutine cofdx
