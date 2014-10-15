!-------------------------------------------------------------------------------
! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files. 
!-------------------------------------------------------------------------------
subroutine save_fields(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use insect_module
  implicit none

  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(in)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect
  
  real(kind=pr):: volume
  character(len=6) :: name

  !--Set up file name base    
  if ( save_only_one_period == "yes" ) then
    ! overwrite files from last period to save disk space
    ! i.e. t=1.05 is written to t=0.05, as well as 2.05 and 3.05
    write(name,'(i6.6)') floor( (time%time-real(floor(time%time/tsave_period)))*1000.d0 )
  else
    ! name is just the time
    write(name,'(i6.6)') floor(time%time*1000.d0) 
  endif
  
  if (mpirank == 0 ) then
    write(*,'("Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A," ...")',advance='no') & 
    time%time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
  endif

  ! Save the velocity
  if (isaveVelocity == 1) then
    call save_field_hdf5(time,"./ux_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"ux")
    call save_field_hdf5(time,"./uy_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"uy")
    call save_field_hdf5(time,"./uz_"//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"uz")
  endif  
  
  ! Save the pressure
  if (isavePress == 1) then  
    call save_field_hdf5(time,'./p_'//name,u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),"p")
  endif
     
  ! Save the Vorticity
  if (isaveVorticity==1) then
    !-- compute vorticity:
    call curl_x( u(:,:,:,1:3), nlk(:,:,:,1:3,1) )
    call save_field_hdf5(time,"./vorx_"//name,nlk(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1,1),"vorx")
    call save_field_hdf5(time,"./vory_"//name,nlk(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2,1),"vory")
    call save_field_hdf5(time,"./vorz_"//name,nlk(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3,1),"vorz")
  endif
      
  ! Save the Mask
  if (isaveMask == 1 .and. iPenalization == 1) then
    ! create mask at current time
    if(iMoving==1) call create_mask( time%time, mask, mask_color, us, Insect, beams )
    ! make sure mask is between 0 and 1 
    mask = mask*eps
    ! warn if we're about to store an empty mask
    call compute_mask_volume(mask,volume)
    if ((mpirank==0).and.(volume<1e-10)) write(*,*) "WARNING: saving empty mask"
    ! save the mask
    call save_field_hdf5(time,'./mask_'//name,mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),"mask")
    mask = mask/eps
  endif
  
  ! save solid velocity
  if (isaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call save_field_hdf5(time,'./usx_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"usx")
    call save_field_hdf5(time,'./usy_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"usy")
    call save_field_hdf5(time,'./usz_'//name,us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"usz")
  endif
  
!   !------------
!   ! passive scalar
!   !-------------
!   if (use_passive_scalar==1) then
!     call ifft ( ink=uk(:,:,:,4), outx=work(:,:,:,1) )
!     call save_field_hdf5(time,'./scalar_'//name,work(:,:,:,1),"scalar")
!   endif
  
  
  if (mpirank==0) write(*,*) " ...DONE!"
end subroutine save_fields


! Write the field field_out to file filename, saving the name of the
! field, dsetname, (as well as time, eps, and resolution) to the
! metadata of the hdf file .
subroutine save_field_hdf5(time,filename,field_out,dsetname)
  use mpi
  use vars
  use hdf5
  implicit none
  
  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  type(timetype), intent (in) :: time
  character(len=*), intent (in) :: filename
  character(len=*), intent (in) :: dsetname

  integer(hid_t) :: file_id   ! file identifier
  integer(hid_t) :: dset_id   ! dataset identifier
  integer(hid_t) :: filespace ! dataspace identifier in file
  integer(hid_t) :: memspace  ! dataspace identifier in memory
  integer(hid_t) :: plist_id  ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  ! hyperslab dimensions
  integer(hsize_t), dimension(rank) :: dimensions_local
  ! chunk dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims
  ! how many blocks to select from dataspace
  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  ! stride is spacing between elements, this is one here
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! HDF attribute variables
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension

  integer :: mpierror, i
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()
  
  !!! Tell HDF5 how our  data is organized:
  dimensions_file = (/nx,ny,nz/)
  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1
  ! each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
  do i = 1, 3
     call MPI_REDUCE(dimensions_local(i),chunking_dims(i),1, &
          MPI_INTEGER8,MPI_MAX,0,&
          MPI_COMM_WORLD,mpierror)
     call MPI_BCAST(chunking_dims(i),1,MPI_INTEGER8,0, &
          MPI_COMM_WORLD,mpierror )
  enddo

  !!! Set up the HDF data structures:

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)
  ! Setup file access property list with parallel I/O access.
  ! this sets up a property list (plist_id) with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store the MPI IO comminucator
  ! information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
#ifdef TURING  
  call h5fcreate_f('bglockless:'//trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
#else
  call h5fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
#endif 
  ! this closes the property list plist_id (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset.
  ! NB: chunking and hyperslab are unrelated
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
  ! Output files are single-precition
  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, &
       dset_id, error, plist_id)
  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error, stride, dimensions_local)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively, double precision in memory
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
       error, file_space_id = filespace, mem_space_id = memspace,&
       xfer_prp = plist_id)

  !!! Write the attributes to the HDF files.
  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
  adims = (/1/)
  call write_attribute_dble(adims,"time",(/time%time/),1,dset_id)
  call write_attribute_dble(adims,"epsi",(/eps/),1,dset_id)
  adims = (/3/)
  call write_attribute_dble(adims,"domain_size",(/xl,yl,zl/),3,dset_id)
  call write_attribute_int(adims,"nxyz",(/nx,ny,nz/),3,dset_id)

  !!! Close dataspaces:
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error) ! Close the dataset.
  call h5pclose_f(plist_id, error) ! Close the property list.
  call h5fclose_f(file_id, error) ! Close the file.
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

  ! write the XMF data for all of the saved fields
  if ((mpirank==0).and.(isaveXMF==1)) then
     ! the filename contains a leading "./" which we must remove
     call Write_XMF(time%time,&
          trim(adjustl(filename(3:len(filename)))),&
          trim(adjustl(dsetname))&
          )
  endif

  time_hdf5=time_hdf5 + MPI_wtime() - t1 ! performance analysis
end subroutine save_field_hdf5


! Generate an XMF file for paraview.  Note: this is a single scalar
! field, no time-stepping or vectors are available but this allows to
! directly copy-paste a single field and load it into paraview without
! any effort.
subroutine Write_XMF(time,filename,dsetname)
  use vars
  implicit none

  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename, dsetname
  character(len=128) :: tmp_time, tmp_nxyz

  write(tmp_time,'(es15.8)') time
  ! note index changes. paraview requirement.
  write(tmp_nxyz,'(3(i4,1x))') nz,ny,nx

  ! note the XMF file goes also in the fields/ directory
  open (14, file='./'//trim(adjustl(filename))//'.xmf', status='replace')

  write(14,'(A)') '<?xml version="1.0" ?>'
  write(14,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(14,'(A)') '<Xdmf Version="2.0">'
  write(14,'(A)') '<Domain>'
  write(14,'(A)') '<Grid Name="FLUSI_cartesian_grid" GridType="Uniform">'
  write(14,'(A)') '    <Time Value="'//trim(adjustl(tmp_time))//'" />'
  write(14,'(A)') '    <Topology TopologyType="3DCoRectMesh" Dimensions="'&
                //trim(adjustl(tmp_nxyz))//'"/>'
  write(14,'(A)') ' '
  write(14,'(A)') '    <Geometry GeometryType="Origin_DxDyDz">'
  write(14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write(14,'(A)') '    0 0 0'
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write(14,'(4x,3(es15.8,1x))') dx, dy, dz
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    </Geometry>'
  write(14,'(A)') ' '
  write(14,'(A)') '    <Attribute Name="'//trim(adjustl(dsetname)) &
                //'" AttributeType="Scalar" Center="Node">'
  write(14,'(A)') '    <DataItem Dimensions="'//trim(adjustl(tmp_nxyz))&
                //'" NumberType="Float" Format="HDF">'
  write(14,'(A)') '    '//trim(adjustl(filename))//'.h5:/'&
                //trim(adjustl(dsetname))
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    </Attribute>'
  write(14,'(A)') '</Grid>'
  write(14,'(A)') '</Domain>'
  write(14,'(A)') '</Xdmf>'

  close (14)
end subroutine Write_XMF


! write a runtime backup to disk. We're saving the fluid state vector 
! u=(/ux,uy,uz,p/) in one HDF5 file and some FSI/insect information in separate 
! files. The hdf5 file is written in double precision, to allow RK2/RK4 schemes
! to resume exactly (without loss of precision) where it left. AB2 schemes are
! different - they apply the startup step again (this associated with loos of 
! precision)
subroutine dump_runtime_backup(time,nbackup,u,Insect,beams)
  use vars
  use hdf5
  use solid_model
  use insect_module
  implicit none
  
  type(timetype),intent(inout) :: time
  integer,intent(inout) :: nbackup
  real(kind=pr),intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect

  character(len=18) :: filename
  real(kind=pr) :: t1
  integer :: error  ! error flags
  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: plist_id      ! property list identifier

  t1=MPI_wtime() ! performance diagnostic

  ! Create name for the backup file. We keep at any time at most 2 sets of 
  ! backups, runtime_backupX.h5 with X=0,1  
  if(mpirank == 0) then
     write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")',&
     advance='no') nbackup, time%time
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".h5")') nbackup

  ! Initialize HDF5 library and Fortran interfaces:
  call h5open_f(error)

  !!! Setup file access property list with parallel I/O access.
  ! Set up a property list ("plist_id") with standard values for
  ! FILE_ACCESS:
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store MPI IO comminucator information
  ! in the file access property list:
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
  call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, &
       access_prp = plist_id)
  ! Close the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Write the fluid backup field:
  call dump_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"ux",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call dump_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"uy",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call dump_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"uz",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call dump_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),"p",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)

       
  ! Close the file:
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library:
  call h5close_f(error)
  
  
  !-------------------------------------------------------------------------
  ! backup for the rigid body solver (free-flight insect)
  !-------------------------------------------------------------------------
  if ((method=="fsi").and.(mpirank==0)) then
  if (iMask=="Insect") then
  if ((Insect%BodyMotion=="takeoff").and.(Insect%KineFromFile=="simplified_dynamic")) then
    write (*,'(A)',advance="no") "insect bckp in "//filename//".rigidsolver"
    open(10, file=filename//".rigidsolver", form='formatted', status='replace') 
    write(10, *) SolidDyn%var_new, SolidDyn%var_this,&
                 SolidDyn%rhs_this, SolidDyn%rhs_old
    close(10)
  endif
  endif
  endif
  
  !-------------------------------------------------------------------------
  !-- backup solid solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_solid_model=="yes").and.(method=="fsi")) then
    call dump_solid_backup( time%time, beams, nbackup )
  endif    
  
  nbackup = 1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1 ! Performance diagnostic

  if(mpirank == 0) write(*,'(A)') "...DONE!"
end subroutine dump_runtime_backup


! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "file_id". Attributes are stores in one attribute
! "bckp" which contains 8 values
subroutine dump_field_backup(field,dsetname,time,dt0,dt1,n1,it,file_id)
  use mpi
  use vars
  use hdf5
  implicit none

  real(kind=pr),intent(in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname
  integer,intent(in) :: n1,it
  integer(hid_t), intent(in) :: file_id       ! file identifier

  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)

  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error, mpierror, i  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes

  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  ! Offsets
  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)

  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
  do i = 1, 3
     call MPI_REDUCE(dimensions_local(i), chunking_dims(i),1, &
          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST(chunking_dims(i), 1, MPI_INTEGER8, 0, &
          MPI_COMM_WORLD, mpierror)
  enddo

  ! -----------------------------------
  ! Create the data space for the  dataset.
  ! -----------------------------------
  ! Dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)
  ! Create chunked dataset.
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
  call h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, plist_id)
  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively.
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_file, error, &
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

  ! ------
  ! Attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
  ! ------
  adims = (/8/)
  allocate (attributes(1:8))
  aname = "bckp"
  attributes = (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/)
  call write_attribute_dble(adims,aname,attributes,8,dset_id)

  ! Close dataspaces, dataset and property list
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error)
  call h5pclose_f(plist_id, error)

  deallocate(attributes)
end subroutine dump_field_backup





! Read in a single file that follows the naming convention
! this is a serial routine (parallel version below)
! note you need to know what dimension the file has,
! call fetch_attributes first
subroutine read_single_file_serial(filename,field)
  use mpi
  use vars
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr), intent (out) :: field(0:nx-1,0:ny-1,0:nz-1)
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  character(len=80)             :: dsetname

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  
  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )

  if (mpisize>1) then
    write (*,*) "this routine is currently serial only"
    call abort()
  endif
  
  ! check if file exist
  call check_file_exists( filename )

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)
  
  ! Definition of memory distribution
  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = nx
  dimensions_local(2) = ny
  dimensions_local(3) = nz

  offset(1) = 0
  offset(2) = 0
  offset(3) = 0
  
  chunking_dims(1) = nx
  chunking_dims(2) = ny
  chunking_dims(3) = nz

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id,error)
  call H5close_f(error)
  
  
end subroutine read_single_file_serial




! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine Read_Single_File ( filename, field )
  use mpi
  use vars
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),&
  dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
  intent (out) :: field
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: nx_file, ny_file, nz_file, mpierror, i  

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1

  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )
  if (mpirank==0) then
    write (*,'("Reading file ",A,"  .....")',advance='no') trim(adjustl(filename))
  endif
  
  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )
  
  ! fetch attributes from file to see if it is a good idea to load it
  call Fetch_attributes( filename, dsetname,nx_file,ny_file,nz_file,& 
                         xl_file,yl_file ,zl_file,time )
         
  ! if the resolutions do not match, yell and hang yourself       
  if ((nx.ne.nx_file).or.(ny.ne.ny_file).or.(nz.ne.nz_file)) then                         
    if (mpirank == 0) then
    write (*,'(A)') "read_single_file: ERROR " // trim(filename)
    write (*,'("nx=",i4,"ny=",i4,"nz=",i4)') nx,ny,nz
    write (*,'("but in file: nx=",i4,"ny=",i4,"nz=",i4)') nx_file,ny_file,nz_file
    call abort()
    endif
  endif
  
  ! if the domain size doesn't match, proceed, but yell.
  if ((xl.ne.xl_file).or.(yl.ne.yl_file).or.(zl.ne.zl_file)) then                         
    if (mpirank == 0) then
    write (*,'(A)') "read_single_file: WARNING " // trim(filename)
    write (*,'("xl=",es12.4,"yl=",es12.4,"zl=",es12.4)')&
      xl,yl,zl
    write (*,'("but in file: xl=",es12.4,"yl=",es12.4,"zl=",es12.4)') & 
      xl_file,yl_file,zl_file
    write (*,'(A)') "proceed, with fingers crossed."
    endif
  endif
  
  !-----------------------------------------------------------------------------
  ! load the file
  !-----------------------------------------------------------------------------  
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)
  
  ! Definition of memory distribution
  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)
  
  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
  do i = 1, 3
     call MPI_REDUCE ( dimensions_local(i), chunking_dims(i),1, &
          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST  ( chunking_dims(i), 1, MPI_INTEGER8, 0, &
          MPI_COMM_WORLD, mpierror )
  enddo

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id,error)
  call H5close_f(error)
  
  if (mpirank==0) then
    write (*,'("...DONE! ")',advance='yes')
  endif
  
end subroutine Read_Single_File




! Load backup data from disk to initialize run for restart
subroutine read_runtime_backup(filename,time,u,Insect,beams)
  use vars
  use hdf5
  use solid_model
  use insect_module
  implicit none

  character(len=*),intent(in) :: filename
  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  type(solid),dimension(1:nBeams),intent(in) :: beams
  type(diptera),intent(in) :: Insect

  integer :: error  ! Error flag
  integer(hid_t) :: file_id       ! File identifier
  integer(hid_t) :: plist_id      ! Property list identifier

  if(mpirank == 0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! I'm trying to resume a backup file: "//filename
  endif
  
  call check_file_exists ( filename )

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Read fluid backup field:
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),"ux",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),"uy",&
       time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),"uz",&
        time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)
  call read_field_backup(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),"p",&
        time%time,time%dt_old,time%dt_new,time%n1,time%it,file_id)       

  ! Close the file:
  call h5fclose_f(file_id,error)
  ! Close FORTRAN interfaces and HDF5 library:
  call H5close_f(error)

  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time%time, time%dt_old
     write(*,'("!!! DONE READING BACKUP (thats good news!)")')
     write(*,'("---------")')
  endif

end subroutine read_runtime_backup


! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
subroutine read_field_backup(field,dsetname,time,dt0,dt1,n1,it,file_id)
  use mpi
  use vars
  use hdf5
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), &
       intent(out) :: field
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (out)  :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname
  integer,intent(out)           :: n1,it
  integer                       :: nx_file,ny_file,nz_file, mpierror, i

  integer(hid_t), intent(in) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  !integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes

  ! Definition of memory distribution
  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)

  ! Each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
  do i = 1, 3
     call MPI_REDUCE(dimensions_local(i),chunking_dims(i),1,MPI_INTEGER8,&
          MPI_MAX,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST(chunking_dims(i),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierror)
  enddo

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error) ! note the dataset remains opened

  ! attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
  adims = (/8/)
  allocate (attributes(1:8))
  aname = "bckp"
  call h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
  call h5aget_space_f(attr_id, aspace_id, error)
  call h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attributes, adims, error)
  call h5aclose_f(attr_id, error) ! Close the attribute.
  call h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  time    = attributes(1)
  dt1     = attributes(2)
  dt0     = attributes(3)
  n1      = int(attributes(4))
  it      = int(attributes(5))
  nx_file = int(attributes(6))
  ny_file = int(attributes(7))
  nz_file = int(attributes(8))

  if ( (nx_file.ne.nx).or.(nx_file.ne.nx).or.(nx_file.ne.nx)) then
     write (*,'(A)') "!!! Thats odd...the backup you're trying to resume doesn't have the same nx,ny,nz"
     write (*,'(A)') "I'll leave you crying and commit suicide here."
     call abort()
  endif

  deallocate (attributes)

  ! Close dataset
  call h5dclose_f(dset_id, error)
end subroutine read_field_backup


! Write a given attribute with attribute name aname and dimensions
! adims/dims to a given dataset identifier dset_id. Double version.
subroutine write_attribute_dble(adims,aname,attribute,dim,dset_id)
  use mpi
  use vars
  use hdf5
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  real (kind=pr), DIMENSION(dim), intent (in) :: attribute
  character(len=*), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id  ! dataset identifier
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier

  ! Determine the dataspace identifier aspace_id
  call h5screate_simple_f(arank,adims,aspace_id,error)

  ! set attr_id, ie create an attribute attached to the object dset_id
  call h5acreate_f(dset_id,aname,H5T_NATIVE_DOUBLE,aspace_id,attr_id,error)

  ! Write the attribute data attribute to the attribute identifierattr_id
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)

  call h5aclose_f(attr_id, error) ! Close the attribute.
  call h5sclose_f(aspace_id, error) ! Terminate access to the data space.
end subroutine write_attribute_dble

! Write a given attribute with attribute name aname and dimensions
! adims/dims to a given dataset identifier dset_id. Integer version.
subroutine write_attribute_int(adims,aname,attribute,dim,dset_id)
  use mpi
  use vars
  use hdf5
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  integer, DIMENSION(dim), intent (in) :: attribute
  character(len=*), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id  ! dataset identifier
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id   ! Attribute identifier

  ! Determine the dataspace identifier aspace_id
  call h5screate_simple_f(arank,adims,aspace_id,error)

  ! set attr_id, ie create an attribute attached to the object dset_id
  call h5acreate_f(dset_id,aname,H5T_NATIVE_INTEGER,aspace_id,attr_id,error)

  ! Write the attribute data attribute to the attribute identifier attr_id.
  call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attribute,adims,error)

  call h5aclose_f(attr_id,error) ! Close the attribute.
  call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
end subroutine write_attribute_int

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
subroutine Fetch_attributes( filename, dsetname,  nx, ny, nz, xl, yl ,zl, time )
  use hdf5
  implicit none

  integer, parameter :: pr = 8
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time

  character(len=*) :: filename  ! file name
  character(len=*) :: dsetname  ! dataset name
  character(len=4) :: aname     ! attribute name
  character(len=11) :: aname2

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: attr_id       ! attribute identifier
  integer(hid_t) :: aspace_id     ! attribute dataspace identifier

  real (kind=pr) ::  attr_data  ! attribute data
  real (kind=pr), dimension (1:3) :: attr_data2
  integer, dimension (1:3) :: attr_data3

  integer     ::   error ! error flag
  integer(hsize_t), dimension(1) :: data_dims


  call check_file_exists ( filename )
  
  ! Initialize FORTRAN interface.
  CALL h5open_f(error)

  ! Open an existing file.
  CALL h5fopen_f (filename, H5F_ACC_RDWR_F, file_id, error)
  ! Open an existing dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)


  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (time)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  aname = "time"
  CALL h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 1
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

  time = attr_data
  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (domain_length)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  aname2 = "domain_size"
  CALL h5aopen_f(dset_id, aname2, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 3
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data2, data_dims, error)

  xl = attr_data2(1)
  yl = attr_data2(2)
  zl = attr_data2(3)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! open attribute (sizes)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  aname = "nxyz"
  CALL h5aopen_f(dset_id, aname, attr_id, error)

  ! Get dataspace and read
  CALL h5aget_space_f(attr_id, aspace_id, error)
  data_dims(1) = 3
  CALL h5aread_f( attr_id, H5T_NATIVE_INTEGER, attr_data3, data_dims, error)

  nx = attr_data3(1)
  ny = attr_data3(2)
  nz = attr_data3(3)

  CALL h5aclose_f(attr_id, error) ! Close the attribute.
  CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.

  CALL h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
  CALL h5fclose_f(file_id, error) ! Close the file.
  CALL h5close_f(error)  ! Close FORTRAN interface.

  !   write (*,'("time=",es12.4," domain=",3(es12.4,1x),2x,3(i3,1x))') time, xl,yl,zl, nx, ny, nz
end subroutine Fetch_attributes



! checks if a given file ("fname") exists. if not, code is stopped brutally
subroutine check_file_exists(fname)
  use vars
  use mpi
  implicit none
  
  character (len=*), intent(in) :: fname
  logical :: exist1
  integer :: mpicode

  if (mpirank == 0) then
    inquire ( file=fname, exist=exist1 )
    if ( exist1 .eqv. .false.) then
      write (*,'("ERROR! file: ",A," not found")') trim(fname) 
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    endif  
  endif
  
end subroutine check_file_exists


! overwrite and initialize file
subroutine init_empty_file( fname )
  use vars
  implicit none
  character (len=*), intent(in) :: fname
  
  open (15, file=fname,status='replace')
  close(15)
end subroutine