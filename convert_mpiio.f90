!=========================================================
!      Collect the data from multiple files
!      produced with HIVE parallel code
!      and save them in one file per time step
!      in a way suitable for visualization using Vapor
!=========================================================



!----------------------------------
! Version 2 (Thomas):
! now writes into *.binary and *.ascii because I hated the "old_*" nomenclature
!----------------------------------


program convert_mpiio
  use mpi
  implicit none

  integer, parameter :: pr_in = 8, pr_out = 8

  ! single precision array for output
  integer, parameter :: mpireal = MPI_DOUBLE_PRECISION 

  integer :: nx, ny, nz, filedesc, mpicode
  integer :: ix, iy, iz
  integer :: record_length
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus
  real (kind=pr_in) :: allmin, allmax
  real (kind=pr_out), dimension (:,:,:), allocatable :: field
  real (kind=4), dimension (:,:,:), allocatable :: field2
  real (kind=pr_in) :: time
  character (len=16) :: fname
  character (len=80) :: ddname, basename
  character (len=16) :: endian_type

  write (*,*) '--------------------------------------'
  
  read *, nx
  read *, ny
  read *, nz
  read *, time
  read *, basename
  read *, endian_type

  print *, 'Start collect at t=', time

  write (fname, '(es10.4)') time

  call MPI_INIT (mpicode)

  allocate ( field(0:(nx-1),0:(ny-1),0:(nz-1)) )

  call MPI_FILE_OPEN (MPI_COMM_WORLD,basename(1:len_trim(basename))//'mpiio_'//fname,MPI_MODE_RDONLY,MPI_INFO_NULL,filedesc,mpicode)
  call MPI_FILE_READ_ORDERED (filedesc,field,nx*ny*nz,mpireal,mpistatus,mpicode)
  call MPI_FILE_CLOSE (filedesc,mpicode)


  allmax = maxval (field)
  allmin = minval (field)

  write(*,'("Min:Max= ",es12.4,1x,es12.4,1x,"nx:ny:nz = ",i3,1x,i3,1x,i3)') allmin, allmax,nx,ny,nz

  allocate ( field2(0:(nx-1),0:(ny-1),0:(nz-1)) )
  field2 = field
  ! ------------------------------------
  ! Write *,binary file (for vapor et al)
  ! ------------------------------------
  write (fname, '(es10.4)') time
  inquire (iolength=record_length) field
  open (12, file = 'old_'//basename(1:len_trim(basename))//trim(fname), form='unformatted', status='new')
  write (12) (((field2 (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
  close (12)

  ! ------------------------------------
  ! Write *,ascii file (for matlab)
  ! ------------------------------------
  open (11, file = basename(1:len_trim(basename))//trim(fname)//'.ascii', form='formatted', status='new')
  do ix = 0,nx-1
  do iy = 0,ny-1
  do iz = 0,nz-1
  write (11,'(es15.8)') field2 (ix,iy,iz) !es15.8 defines output precision
  enddo
  enddo
  enddo
  close (11)

  deallocate ( field )

  print *, 'File ','full_'//basename(1:len_trim(basename))//fname, ' written'
  
  write (*,*) '--------------------------------------'
  
  call MPI_FINALIZE (mpicode)

end program convert_mpiio
