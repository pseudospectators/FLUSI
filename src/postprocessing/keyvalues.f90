!-------------------------------------------------------------------------------
! ./flusi --postprocess --keyvalues mask_00000.h5
!-------------------------------------------------------------------------------
! load the specified *.h5 file and creates a *.key file that contains
! min / max / mean / L2 norm of the field data. This is used for unit testing
! so that we don't need to store entire fields but rather the *.key only
subroutine keyvalues(filename)
  use vars
  use mpi
  implicit none
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time, npoints, q, x,y,z
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer :: ix,iy,iz

  if (mpisize>1) then
    write (*,*) "--keyvalues is currently a serial version only, run it on 1CPU"
    call abort()
  endif

  call check_file_exists( filename )

  write (*,*) "analyzing file "//trim(adjustl(filename))//" for keyvalues"

  !---------------------------------------------------------
  ! in a first step, we fetch the attributes from the dataset
  ! namely the resolution is whats important
  ! this routine was created in the mpi2vis repo -> convert_hdf2xmf.f90
  !---------------------------------------------------------
  call fetch_attributes( filename, nx, ny, nz, xl, yl, zl, time , nu )
  write(*,'("File is at time=",es12.4)') time
  allocate ( field(0:nx-1,0:ny-1,0:nz-1) )

  ra=(/0,0,0/)
  rb=(/nx-1,ny-1,nz-1/)
  call read_single_file (filename, field)

  npoints = dble(nx)*dble(ny)*dble(nz)

  ! compute an additional quantity that depends also on the position
  ! (the others are translation invariant)
  q=0.d0
  do iz = 0, nz-1
   do iy = 0, ny-1
    do ix = 0, nx-1
      x = dble(ix)*xl/dble(nx)
      y = dble(iy)*yl/dble(ny)
      z = dble(iz)*zl/dble(nz)

      q = q + x*field(ix,iy,iz) + y*field(ix,iy,iz) + z*field(ix,iy,iz)
      enddo
    enddo
  enddo

  open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
  write (14,'(6(es15.8,1x))') time, maxval(field), minval(field),&
   sum(field)/npoints, sum(field**2)/npoints, q/npoints
  write (*,'(6(es15.8,1x))') time, maxval(field), minval(field),&
   sum(field)/npoints, sum(field**2)/npoints, q/npoints
  close (14)

  deallocate (field)
end subroutine keyvalues
