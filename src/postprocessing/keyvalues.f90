!-------------------------------------------------------------------------------
! ./flusi --postprocess --keyvalues mask_00000.h5
!-------------------------------------------------------------------------------
! load the specified *.h5 file and creates a *.key file that contains
! min / max / mean / L2 norm of the field data. This is used for unit testing
! so that we don't need to store entire fields but rather the *.key only
subroutine keyvalues(filename)
  use vars
  use mpi
  use p3dfft_wrapper
  use module_helpers
  
  implicit none
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time, npoints, q, x,y,z
  real(kind=pr) :: maxi,mini,squari,meani,qi
  real(kind=pr) :: maxl,minl,squarl,meanl,ql
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer :: ix,iy,iz,mpicode

  call check_file_exists( filename )

  if (root) write (*,*) "analyzing file "//trim(adjustl(filename))//" for keyvalues"

  !---------------------------------------------------------
  ! in a first step, we fetch the attributes from the dataset
  ! namely the resolution is whats important
  ! this routine was created in the mpi2vis repo -> convert_hdf2xmf.f90
  !---------------------------------------------------------
  call fetch_attributes( filename, nx, ny, nz, xl, yl, zl, time , nu, origin )
  ! initialize code and scaling factors for derivatives
  call decomposition_initialize()
  allocate(field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file (filename, field)

  npoints = dble(nx)*dble(ny)*dble(nz)

  ! compute an additional quantity that depends also on the position
  ! (the others are translation invariant)
  ql = 0.d0
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
      x = dble(ix)*xl/dble(nx)
      y = dble(iy)*yl/dble(ny)
      z = dble(iz)*zl/dble(nz)

      ql = ql + x*field(ix,iy,iz) + y*field(ix,iy,iz) + z*field(ix,iy,iz)
      enddo
    enddo
  enddo

  maxl = maxval(field)
  minl = minval(field)
  squarl = sum(field**2)
  meanl  = sum(field)

  call MPI_ALLREDUCE (ql,qi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (maxl,maxi,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (minl,mini,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (squarl,squari,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (meanl,meani,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)

  qi = qi / npoints
  squari = squari / npoints
  meani = meani / npoints

  if (root) then
    open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
    write (14,'(6(es15.8,1x))') time, maxi, mini, meani, squari, qi
    write (*,'(A)') "Result:"
    write (* ,'(6(A15,1x))') "time","maxval","minval","meanval","sumsquares","Q-integral"
    write (* ,'(6(es15.8,1x))') time, maxi, mini, meani, squari, qi
    write (*,'(A)') "These values can be used to compare two HDF5 files"
    close (14)
  endif

  deallocate (field)
end subroutine keyvalues
