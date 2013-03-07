!--------------------------------------
! Compute extents of transposed arrays
!--------------------------------------
subroutine trextents ( idir, n, mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  implicit none

  integer, intent(in) :: idir
  integer, dimension (1:3), intent(in) :: n
  integer, dimension (2), intent(in) :: mpidims, mpicoords
  integer, dimension (1:3), intent(out) :: ka, kb, ks, kat, kbt, kst

  ka(1) = 0
  kb(1) = n(1)-1
  ka(2) = mpicoords(2)*n(2)/mpidims(2)
  kb(2) = (mpicoords(2)+1)*n(2)/mpidims(2)-1
  ka(3) = mpicoords(1)*n(3)/mpidims(1)
  kb(3) = (mpicoords(1)+1)*n(3)/mpidims(1)-1
  ks(:) = kb(:)-ka(:)+1

  kat(1) = 0
  kbt(1) = n(idir+1)-1
  kat(1+idir) = mpicoords(3-idir)*n(1)/mpidims(3-idir)
  kbt(1+idir) = (mpicoords(3-idir)+1)*n(1)/mpidims(3-idir)-1
  kat(4-idir) = ka(4-idir)
  kbt(4-idir) = kb(4-idir)
  kst(:) = kbt(:)-kat(:)+1

end subroutine trextents