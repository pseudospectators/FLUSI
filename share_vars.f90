module share_vars
  implicit none

  integer, parameter :: pr = 8
  integer, save :: nx, ny, nz, nt
  integer, save :: mpisize, mpirank, mpicommcart, mpiinteger, mpireal, mpicomplex
  integer, save :: iDealias, inicond, imask
  integer, save :: iSaveVelocity, iSaveVorticity, iSavePress, iSaveMask
  integer, save :: iTimeMethodFluid, iMoving, iPenalization, iMeanFlow, iDoBackup, iSaveSolidVelocity
  integer, dimension (2), save :: mpidims, mpicoords, mpicommslab
  integer, dimension (1:3), save :: ra, rb, rs, ca, cb, cs
  integer, dimension (:,:), allocatable, save :: ra_table, rb_table
  real (kind=pr), save :: tmax, cfl, nu, eps, pi, scalex, scaley, scalez, length, xl, yl, zl, dx, dy ,dz
  real (kind=pr), save :: Ux, Uy, Uz, Ax, Ay, Az, tstart, tsave, tdrag
  real (kind=pr), save :: x0, y0, z0
  real (kind=pr), dimension (:,:,:), allocatable, save :: dealiase, mask
  real (kind=pr), dimension (:,:,:,:), allocatable, save :: us

  integer, parameter :: nlines=256 ! maximum number of lines in PARAMS-file
end module share_vars

