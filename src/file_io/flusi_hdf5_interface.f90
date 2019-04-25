!*******************************************************************************
! FLUSI layer for the hdf interface
! The basic HDF5 routines are in the hdf5_wrapper module (like writing a field to file or storing attributes)

! this compatibility layer only provides a flusi-specific interface for saving and reading
! a field from file.
! Our hdf5-file follow the convention that a file ux_00005.h5
!   * contains only one underscore _ in the name
!   * the first part until the underscore is the datasetname we expect to find in that file
!     i.e. ux_00005.h5/ux
!   * the second part of the filename until the dot can be anything, but usually is the time stamp
!   * each dataset contains the attributes nxyz (resolution), time, domain size, viscosity and penalization parameter
!     (the latter two are not indispensable)
!*******************************************************************************

!-------------------------------------------------------------------------------
! Standart wrapper for the HDF5 library, saves a single array (e.g. one vector component)
! to a single HDF5 file. A bunch of useful attributes (resolution, domain size,
! penalization, viscosity, time) are stored as well.
!-------------------------------------------------------------------------------
subroutine save_field_hdf5(time,filename,field_out)
  use vars
  use hdf5_wrapper
  use module_helpers
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(inout) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
  character(len=strlen) :: fname,fname2
  integer, dimension(1:3) :: tmp
  real(kind=pr) :: t1, mbyte
  integer :: mpicode

  t1 = MPI_wtime()
  ! check if the file name contains the suffix *.h5
  ! if not, add it
  if (index(filename,'.h5')==0 ) then
    ! file does not contain *.h5 ending -> add suffix
    fname = trim(adjustl(filename))//'.h5'
  else
    fname = trim(adjustl(filename))
  endif

  ! header
  if (root) then
    write(*,'("Writing to ",A," dset=",A," stride=",i1," ...")',advance='no') &
    trim(adjustl(fname)), trim(adjustl(get_dsetname(fname))), striding
  endif

  if (striding<1) call abort(8881, "Striding value is bad, exit!")

  if (striding==1) then
    ! save the entire field to disk (no striding)
    ! write actual field to the file
    call write_field_hdf5( fname, get_dsetname(fname), ra, rb, field_out)
    ! append some useful attributes to the field in the file
    call write_attribute( fname, get_dsetname(fname), "time",(/time/))
    call write_attribute( fname, get_dsetname(fname), "viscosity",(/nu/))
    call write_attribute( fname, get_dsetname(fname), "epsi",(/eps/))
    call write_attribute( fname, get_dsetname(fname), "domain_size",(/xl,yl,zl/))
    call write_attribute( fname, get_dsetname(fname), "origin",origin)
    call write_attribute( fname, get_dsetname(fname), "nxyz",(/nx,ny,nz/))
  else
    ! save strided field to disk
    call save_field_hdf5_strided(time,fname,field_out)
  endif

  ! footer
  t1 = MPI_wtime() - t1
  mbyte = dble(nx/striding)*dble(ny/striding)*dble(nz/striding)*4.d0/1024.0d0/1024.0d0
  if (root) write(*,'(".. wrote ",f9.2," MB in ",f9.2," s (",f9.2,"MB/s)")') &
  mbyte, t1, mbyte/t1

  ! if the last thing the code does is writing, and afterwards an ABORT is called,
  ! the writing may get interrupted. it is not a big performance loss, so here we have
  ! a barrier.
  call MPI_barrier(MPI_COMM_WORLD, mpicode)

end subroutine save_field_hdf5





!-------------------------------------------------------------------------------
! save a strided field to HDF5. this is certainly not the most elegant way to do
! it, nor the most general, but it works.
!-------------------------------------------------------------------------------
! we figure out what the array bounds of the downsampled array on the local CPU
! are, then copy the data to the smaller field, and then pass both to the HDF5
! wrapper and write it to disk
!-------------------------------------------------------------------------------
subroutine save_field_hdf5_strided(time,fname,field_out)
  use module_helpers
  use vars
  use hdf5_wrapper
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: fname
  real(kind=pr),dimension(:,:,:), allocatable :: field_red

  integer :: ix,iy,iz,ixred,iyred,izred
  integer,dimension(1:3) :: rared, rbred
  integer :: ixmin,ixmax,ixstride,iystride,izstride,iymin,izmin,iymax,izmax

  ! do not touch lower/upper bounds, code does NOT work on arbitrary subsets.
  ! ONLY striding of the FULL field is possible.
  ixmin=0; ixmax=nx-1; ixstride=striding;
  iymin=0; iymax=ny-1; iystride=striding;
  izmin=0; izmax=nz-1; izstride=striding;

  ixred=0; rared(1) = nx*99; rbred(1) = -99
  iyred=0; rared(2) = ny*99; rbred(2) = -99
  izred=0; rared(3) = nz*99; rbred(3) = -99

  do ix = ixmin,ixmax,ixstride
    ! check if this x-coordinate is on my local memory
    if ( on_proc((/ix,ra(2),ra(3)/)) ) then
      rared(1) = min(rared(1),ixred)
      rbred(1) = max(rbred(1),ixred)
    endif
    ixred=ixred+1
  enddo

  do iy = iymin,iymax,iystride
    ! check if this y-coordinate is on my local memory
    if ( on_proc((/ra(1),iy,ra(3)/)) ) then
      rared(2) = min(rared(2),iyred)
      rbred(2) = max(rbred(2),iyred)
    endif
    iyred=iyred+1
  enddo

  do iz = izmin,izmax,izstride
    ! check if this z-coordinate is on my local memory
    if ( on_proc((/ra(1),ra(2),iz/)) ) then
      rared(3) = min(rared(3),izred)
      rbred(3) = max(rbred(3),izred)
    endif
    izred=izred+1
  enddo

  allocate( field_red(rared(1):rbred(1),rared(2):rbred(2),rared(3):rbred(3)) )

  ! copy
  do ixred = rared(1),rbred(1)
    do iyred = rared(2),rbred(2)
      do izred = rared(3),rbred(3)
        field_red(ixred,iyred,izred) = field_out(ixstride*ixred,iystride*iyred,izstride*izred)
      enddo
    enddo
  enddo

  call write_field_hdf5( fname, get_dsetname(fname), rared, rbred, field_red)
  ! append some useful attributes to the field in the file
  call write_attribute( fname, get_dsetname(fname), "time",(/time/))
  call write_attribute( fname, get_dsetname(fname), "viscosity",(/nu/))
  call write_attribute( fname, get_dsetname(fname), "epsi",(/eps/))
  call write_attribute( fname, get_dsetname(fname), "domain_size",(/xl,yl,zl/))
  call write_attribute( fname, get_dsetname(fname), "origin",origin)
  call write_attribute( fname, get_dsetname(fname), "nxyz", (/nx,ny,nz/) / striding )

  deallocate (field_red)
end subroutine save_field_hdf5_strided


!----------------------------------------------------
! This routine fetches the resolution, the domain size and the time
! form a *.h5 file
!----------------------------------------------------
! filename: a *.h5 file to read from.
! dsetname: a dataset inside the file. in our system, we only have one per file
!           and this matches the prefix: mask_00010.h5  --> dsetname = "mask"
! note:
!           the file must contain the dataset
!           but especially the attributes "nxyz", "time", "domain_size"
!----------------------------------------------------
subroutine Fetch_attributes( filename, nx, ny, nz, xl, yl ,zl, time, viscosity, origin )
  use hdf5_wrapper
  use module_helpers
  use mpi
  implicit none

  character(len=*), intent(in) :: filename  ! file name
  integer, intent (out) :: nx, ny, nz
  real(kind=pr), intent(out) :: xl,yl,zl, time, viscosity, origin(1:3)

  real(kind=pr),dimension(1) :: attr_data1, attr_data0
  real(kind=pr),dimension(1:3) :: attr_data2
  integer,dimension(1:3) :: attr_data3

  call check_file_exists ( filename )
  call read_attribute( filename, get_dsetname(filename), "time", attr_data0)
  call read_attribute( filename, get_dsetname(filename), "viscosity", attr_data1)
  call read_attribute( filename, get_dsetname(filename), "domain_size", attr_data2)
  call read_attribute( filename, get_dsetname(filename), "origin", origin)
  call read_attribute( filename, get_dsetname(filename), "nxyz", attr_data3)

  time = attr_data0(1)
  viscosity = attr_data1(1)
  xl = attr_data2(1)
  yl = attr_data2(2)
  zl = attr_data2(3)
  nx = attr_data3(1)
  ny = attr_data3(2)
  nz = attr_data3(3)
end subroutine Fetch_attributes



!-------------------------------------------------------------------------------
! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
!-------------------------------------------------------------------------------
subroutine Read_Single_File ( filename, field )
  use vars
  use hdf5_wrapper
  use basic_operators, only : fieldmax, fieldmin, fieldmean, checknan
  use module_helpers
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(out) :: field
  integer, dimension(1:3) :: nxyz
  real(kind=pr), dimension(1:3) :: domain
  real(kind=pr), dimension(1) :: ttime, viscosity_dummy
  real(kind=pr) :: fmax,fmin,favg,t1,mbyte,t2

  t1 = MPI_wtime()

  call check_file_exists ( filename )
  call read_attribute( filename,get_dsetname(filename),"nxyz",nxyz)
  call read_attribute( filename,get_dsetname(filename),"domain_size",domain)
  call read_attribute( filename,get_dsetname(filename),"time",ttime)
  call read_attribute( filename,get_dsetname(filename),"viscosity",viscosity_dummy)
  call read_attribute( filename,get_dsetname(filename), "origin", origin)

  if (mpirank==0) then
    write(*,'(40("~"))')
    write(*,'("Reading from file ",A)') trim(adjustl(filename))
    write(*,'("dsetname=",A)') trim(adjustl(get_dsetname(filename)))
    write(*,'("nx=",i4," ny=",i4," nz=",i4," time=",g12.4," viscosity=",g16.4)') nxyz,ttime(1),viscosity_dummy(1)
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') domain
    write(*,'("x0=",g12.4," y0=",g12.4," z0=",g12.4," (origin)")') origin

    ! if the domain size doesn't match, proceed, but yell.
    if ((xl.ne.domain(1)).or.(yl.ne.domain(2)).or.(zl.ne.domain(3))) then
        write (*,'(A)') " WARNING! Domain size mismatch."
        write (*,'("in memory:   xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') xl,yl,zl
        write (*,'("but in file: xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') domain
        write (*,'(A)') "proceed, with fingers crossed."
    endif

    ! if the resolutions do not match, yell and hang yourself
    if ((nx/=nxyz(1)).or.(ny/=nxyz(2)).or.(nz/=nxyz(3))) then
      write (*,'(A)') "ERROR! Resolution mismatch"
      write (*,'(A)') "This happens if ra(:) and rb(:) are not properly initialized."
      write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
      write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nxyz
      call abort(125,"ERROR! Read_single_file: Resolution mismatch")
    endif
  endif


  ! actual reading of file
  call read_field_hdf5 ( filename, get_dsetname(filename), ra,rb, field)
  ! check if field contains NaN
  call checknan(field,"recently loaded field")

  fmax = fieldmax(field)
  fmin = fieldmin(field)
  favg = fieldmean(field)
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1024.0d0/1024.0d0
  t2 = MPI_wtime() - t1

  if (mpirank==0) then
    write (*,'("read ",f9.2," MB in ",f9.2," s (",f9.2,"MB/s)")') mbyte,t2,mbyte/t2
    write (*,'("max=",g12.4," min=",g12.4," mean=",g12.4)') fmax,fmin,favg
    write (*,'("Done reading file, Elapsed time=",g12.4,"s")') MPI_wtime() - t1
    write(*,'(40("~"))')
  endif

  if (nxyz(1)>1) then
    !3d run
    if ( mod(nxyz(1),2) /= 0 .or. mod(nxyz(2),2) /= 0 .or. mod(nxyz(2),2) /= 0 ) then
      write(*,*) "Warning-Warning-Warning"
      write(*,*) "WARNING UNEVEN NUMBER OF POINTS!! CODE MAY CRASH.", nxyz
    endif
  else
    ! 2d run
    if (mod(nxyz(2),2) /= 0 .or. mod(nxyz(2),2) /= 0 ) then
      write(*,*) "Warning-Warning-Warning"
      write(*,*) "WARNING UNEVEN NUMBER OF POINTS!! CODE MAY CRASH.", nxyz
    endif
  endif

end subroutine Read_Single_File



!-------------------------------------------------------------------------------
! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "filename". Attributes are stores in one attribute
! "bckp" which contains 8 values
!-------------------------------------------------------------------------------
subroutine dump_field_backup(filename,field,dsetname,time,dt0,dt1,n1,it)
  use vars
  use hdf5_wrapper
  use module_helpers
  implicit none

  integer,intent(in) :: n1,it
  real(kind=pr), intent (in) :: time,dt1,dt0
  real(kind=pr),intent(in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*), intent (in) :: dsetname, filename
  real(kind=pr) :: t1, mbyte
  t1 = MPI_wtime()

  if (backup_type == "one-file-backup") then
    !---------------------------------------------------------------------------
    ! all fields go to one HDF5 file
    !---------------------------------------------------------------------------
    ! header
    if (root) then
      write(*,'("Writing to ",A," dset=",A," ...")',advance='no') filename//".h5", dsetname
    endif

    call write_field_hdf5( filename//".h5",dsetname, ra, rb, field, overwrite=.false.)
    call write_attribute( filename//".h5",dsetname,"bckp",&
    (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/) )
  else
    !---------------------------------------------------------------------------
    ! each field goes to one hdf file
    !---------------------------------------------------------------------------
    ! header
    if (root) then
      write(*,'("Writing to ",A," dset=",A," ...")',advance='no') &
      filename//"_"//dsetname//".h5",dsetname
    endif

    call write_field_hdf5( filename//"_"//dsetname//".h5",dsetname, ra, rb, field, overwrite=.true.)
    call write_attribute( filename//"_"//dsetname//".h5",dsetname,"bckp",&
    (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/) )
  endif

  ! footer
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
  use hdf5_wrapper
  use basic_operators, only : checknan
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(out) :: field
  character(len=*), intent (in) :: dsetname, filename
  real(kind=pr)::mbyte,t1
  t1 = MPI_wtime()

  if (mpirank==0) then
    write(*,'("Reading ",A," from backup file ",A)',advance='no') trim(adjustl(dsetname)),trim(adjustl(filename))
  endif

  if (backup_type == "one-file-backup") then
    call read_field_hdf5( filename//".h5", dsetname, ra, rb, field )
  else
    call read_field_hdf5( filename//"_"//dsetname//".h5", dsetname, ra, rb, field )
  endif

  ! check if we read crap
  call checknan(field,'recently read backup file!!')

  mbyte = dble(nx)*dble(ny)*dble(nz)*8.d0/1024.0d0/1024.0d0
  t1 = MPI_wtime() -t1
  if (root) write(*,'(".. read ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1
end subroutine read_field_backup
