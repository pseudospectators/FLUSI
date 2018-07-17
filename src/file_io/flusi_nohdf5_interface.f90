
! Write the field field_out to file filename.
! Replaces hdf5 wrapper routine with the same name,
! if compiled with HDF5FLAG = no
subroutine save_field_hdf5(time,filename,field_out)
  use mpi
  use vars
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  real(kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
  integer :: sz_out(1:3) ! local array size
  character(len=5) :: suffix
  character(len=256) :: fullname
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  integer :: mpicode

  t1 = MPI_wtime()

  !--Set up file name
  write(suffix,'(i5.5)') mpirank
  suffix = trim(adjustl(suffix))
  fullname = trim(adjustl(filename))//'.np'//suffix

  ! Array bounds and sizes
  sz_out(1) = rb(1)-ra(1) +1
  sz_out(2) = rb(2)-ra(2) +1
  sz_out(3) = rb(3)-ra(3) +1

  ! Write Fortran binary file
  if (sz_out(3) > 0) then ! fields with nz==1 are only allocated on some of the mpi ranks
    open(10, file = trim(adjustl(fullname)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3)*pr)
    write (10,rec=1) field_out
    close (10)
  endif

  call MPI_barrier(MPI_COMM_WORLD,mpicode)
  time_hdf5=time_hdf5 + MPI_wtime() - t1 ! performance analysis
end subroutine save_field_hdf5


! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine Read_Single_File ( filename, field )
  use mpi
  use vars
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),&
  dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
  intent (out) :: field

  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  integer :: sz_out(1:3) ! local array size
  character(len=5) :: suffix
  character(len=256) :: fullname
  integer :: mpicode

  !--Set up file name
  write(suffix,'(i5.5)') mpirank
  suffix = trim(adjustl(suffix))
  fullname = trim(adjustl(filename))//'.np'//suffix

  ! Array bounds and sizes
  sz_out(1) = rb(1)-ra(1) +1
  sz_out(2) = rb(2)-ra(2) +1
  sz_out(3) = rb(3)-ra(3) +1

  ! Verbose
  if (mpirank==0) then
    write (*,'("Reading file ",A,"  .....")',advance='no') trim(adjustl(filename))
  endif

  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )

  !-----------------------------------------------------------------------------
  ! load the file
  !-----------------------------------------------------------------------------
  if (sz_out(3) > 0) then
    open(10, file = trim(adjustl(fullname)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3)*pr)
    read (10,rec=1) field
    close (10)
  endif

  call MPI_barrier(MPI_COMM_WORLD,mpicode)
  if (mpirank==0) then
    write (*,'("...DONE! ")',advance='yes')
  endif

end subroutine Read_Single_File


!-------------------------------------------------------------------------------
! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "filename". Attributes are stores in one attribute
! "bckp" which contains 8 values
!-------------------------------------------------------------------------------
subroutine dump_field_backup(filename,field,dsetname,time,dt0,dt1,n1,it)
  use vars
  use module_helpers
  implicit none

  integer,intent(in) :: n1,it
  real(kind=pr), intent (in) :: time,dt1,dt0
  real(kind=pr),intent(in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*), intent (in) :: dsetname, filename
  real(kind=pr) :: t1, mbyte
  integer :: mpicode

  t1 = MPI_wtime()

  if (rb(3)-ra(3)+1 > 0) then
    if (root) then
      write(*,'("Writing to ",A," dset=",A," ...")',advance='no') trim(adjustl(filename)), trim(adjustl(dsetname))
    endif

    write(11) field
  endif

  ! footer
  call MPI_barrier(MPI_COMM_WORLD,mpicode)
  t1 = MPI_wtime() - t1
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1024.0d0/1024.0d0
  if (root) write(*,'(".. wrote ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1
end subroutine dump_field_backup


!-------------------------------------------------------------------------------
! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
!-------------------------------------------------------------------------------
subroutine read_field_backup(filename,dsetname,field)
  use vars
  use basic_operators, only : checknan
  implicit none

  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(out) :: field
  character(len=*), intent (in) :: dsetname, filename
  real(kind=pr)::mbyte,t1
  integer :: mpicode

  t1 = MPI_wtime()

  if (rb(3)-ra(3)+1 > 0) then
    if (mpirank==0) then
      write(*,'("Reading ",A," from backup file ",A)',advance='no') trim(adjustl(dsetname)),trim(adjustl(filename))
    endif

    read(11) field
!    call checknan(field,'recently read backup file!!')
  endif

  ! footer
  call MPI_barrier(MPI_COMM_WORLD,mpicode)
  mbyte = dble(nx)*dble(ny)*dble(nz)*8.d0/1024.0d0/1024.0d0
  t1 = MPI_wtime() -t1
  if (root) write(*,'(".. read ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1

end subroutine read_field_backup


!----------------------------------------------------
! This routine fetches the resolution, the domain size and the time
! THIS IS A STUB. NOT YET IMPLEMENTED WITHOUT HDF5
!----------------------------------------------------
! filename: a *.h5 file to read from.
! dsetname: a dataset inside the file. in our system, we only have one per file
!           and this matches the prefix: mask_00010.h5  --> dsetname = "mask"
! note:
!           the file must contain the dataset
!           but especially the attributes "nxyz", "time", "domain_size"
!----------------------------------------------------
subroutine fetch_attributes( filename, nx, ny, nz, xl, yl ,zl, time, viscosity, origin )
  use module_helpers, only : get_dsetname
  use vars, only : pr,mpirank
  use mpi
  implicit none

  character(len=*), intent(in) :: filename  ! file name
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time, viscosity, origin(1:3)

  real(kind=pr),dimension(1) :: attr_data1, attr_data0
  real(kind=pr),dimension(1:3) :: attr_data2
  integer,dimension(1:3) :: attr_data3

  call check_file_exists ( filename )
  !call read_attribute( filename, get_dsetname(filename), "time", attr_data0)
  !call read_attribute( filename, get_dsetname(filename), "viscosity", attr_data1)
  !call read_attribute( filename, get_dsetname(filename), "domain_size", attr_data2)
  !call read_attribute( filename, get_dsetname(filename), "nxyz", attr_data3)

  if (mpirank==0) write(*,'(A)',advance='no') "WARNING(fetch_attributes): Function not implemented without HDF5 support"

  time = 0
  viscosity = 0
  xl = 0
  yl = 0
  zl = 0
  nx = 0
  ny = 0
  nz = 0
end subroutine Fetch_attributes
