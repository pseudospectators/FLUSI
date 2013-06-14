subroutine Save_scalar_HDF5( time,  filename, field_out, dsetname )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  use HDF5
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(in) :: field_out
  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename   ! file name
  character(len=*), intent (in) :: dsetname 

  integer(hid_t) :: file_id       ! file identifier 
  integer(hid_t) :: dset_id       ! dataset identifier 
  integer(hid_t) :: filespace     ! dataspace identifier in file 
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier 

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file  
  
  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: dimensions_local  

  integer(hsize_t),  dimension(rank) :: count  
  integer(hssize_t), dimension(rank) :: offset 
  integer(hsize_t),  dimension(rank) :: stride
  integer(hsize_t),  dimension(rank) :: block
  integer :: error  ! error flags
  
  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  character(len=4) :: aname ! attribute name
  character(len=11) :: aname2="domain_size" ! attribute name
  
  ! -----------------------------------------------------------
  ! this first part is to tell HDF5 how our (real!) data is organized
  ! -----------------------------------------------------------

  dimensions_file = (/nx,ny,nz/)
  ! the Fortran HDF wrappers have 1-based arrays
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  ! stride is spacing between elements, this is one here
  stride = 1
  ! how many blocks to select from dataspace
  count =  1 

  ! the block contains how many data points to store
  block = dimensions_local

  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)
  ! -----------------------------------------------------------
  ! now we call the HDF subroutines from the "chunks" example
  ! -----------------------------------------------------------
  
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! -----------------------------------------------------------
  ! Setup file access property list with parallel I/O access.
  ! -----------------------------------------------------------
  ! this sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! ----------------------------------
  ! Create the file collectively. (existing files are overwritten)
  ! ---------------------------------
  call H5Fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
       file_id, error, access_prp = plist_id)
  ! this closes the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)

  ! -----------------------------------
  ! Create the data space for the  dataset. 
  ! -----------------------------------
  ! dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset.
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, dimensions_local, error)
  ! double precision in memory...
  call H5Dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, & 
       dset_id, error, plist_id)
  call H5Sclose_f(filespace, error)


  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error, stride, block)

  ! Create property list for collective dataset write
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


  ! Write the dataset collectively. 
  ! but single precision on disk..
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
       error, file_space_id = filespace, mem_space_id = memspace,&
       xfer_prp = plist_id)

       
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
  ! now the data is written, we take care of the ATTRIBUTE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
  adims = (/1/) 
  aname = "time"
  call write_attribute_dble(adims,aname,(/time/),1,dset_id)
  aname = "epsi"
  call write_attribute_dble(adims,aname,(/eps/),1,dset_id)
  adims = (/3/)
  aname2 = "domain_size"
  call write_attribute_dble(adims,aname2,(/xl,yl,zl/),3,dset_id)
  adims = (/3/)
  aname = "nxyz"
  call write_attribute_int(adims,aname,(/nx,ny,nz/),3,dset_id)
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! time to close everything  (EVERYTHING!)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  

  
  ! Close dataspaces.
  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  ! Close the dataset.
  call H5Dclose_f(dset_id, error)
  ! Close the property list.
  call H5Pclose_f(plist_id, error)
  ! Close the file.
  call H5Fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library.
  call h5close_f(error)

  if (mpirank==0) then
     call Write_XMF ( time, trim(adjustl(filename)) , trim(adjustl(dsetname)) )
  endif

end subroutine Save_scalar_HDF5






subroutine Write_XMF ( time, filename, dsetname )
  !--------------------------------------------
  ! this routine generates an XMF file for paraview note: this is a
  ! single scalar field, no time-stepping or vectors are available but
  ! this allows to directly copy-paste a single field and load it into
  ! paraview without any effort.
  !--------------------------------------------
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename, dsetname
  character(len=128) :: tmp_time, tmp_nxyz

  write(tmp_time,'(es15.8)') time
  write(tmp_nxyz,'(3(i4,1x))') nx,ny,nz

  open (14, file=trim(adjustl(filename))//'.xmf', status='replace')
  write(14,'(A)') '<?xml version="1.0" ?>'
  write(14,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(14,'(A)') '<Xdmf Version="2.0">'
  write(14,'(A)') '<Domain>'
  write(14,'(A)') '<Grid Name="FLUSI_cartesian_grid" GridType="Uniform">'  
  write(14,'(A)') '    <Time Value="'//trim(adjustl(tmp_time))//'" />'  
  write(14,'(A)') '    <Topology TopologyType="3DCoRectMesh" Dimensions="'//trim(adjustl(tmp_nxyz))//'"/>'
  write(14,'(A)') ' '
  write(14,'(A)') '    <Geometry GeometryType="Origin_DxDyDz">'
  write(14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write(14,'(A)') '    0 0 0'
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    <DataItem Dimensions="3" NumberType="Float" Format="XML">'
  write(14,'(4x,3(es15.8,1x))') dx, dy, dz!          0.25 0.25 0.25
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    </Geometry>'
  write(14,'(A)') ' '
  write(14,'(A)') '    <Attribute Name="'//trim(adjustl(dsetname))//'" AttributeType="Scalar" Center="Node">'
  write(14,'(A)') '    <DataItem Dimensions="'//trim(adjustl(tmp_nxyz))//'" NumberType="Float" Format="HDF">'
  write(14,'(A)') '    '//trim(adjustl(filename))//'.h5:/'//trim(adjustl(dsetname))
  write(14,'(A)') '    </DataItem>'
  write(14,'(A)') '    </Attribute>'
  write(14,'(A)') '</Grid>'
  write(14,'(A)') '</Domain>'
  write(14,'(A)') '</Xdmf>'

  close (14)

end subroutine Write_XMF




subroutine Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,workvis,work)
  use mpi_header 
  use share_vars
  use hdf5
  implicit none
  character(len=18) :: filename
  real(kind=pr),intent(inout) :: time,dt1,dt0
  integer,intent(inout) :: n1,nbackup,it
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3), intent(inout) :: uk !!!!!!!!!!!!!!!!!!!!!!!!!!!!§§§§§§§§§§§§§§§§§§§§§§§fixme INOUT
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1), intent(inout):: nlk
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),intent(inout) :: workvis
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(inout) :: work
  real(kind=pr) :: t1
  
  integer :: error  ! error flags  
  integer(hid_t) :: file_id       ! file identifier 
  integer(hid_t) :: plist_id      ! property list identifier 
  

  
  t1=MPI_wtime()

  if(mpirank ==0) then
     write(*,'("*** info: time=",es8.2," dumping runtime_backup",i1," to disk....")') time, nbackup
  endif

  ! create current filename
  write(filename,'("runtime_backup",i1,".h5")') nbackup 

  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error) 
  
  ! -----------------------------------------------------------
  ! Setup file access property list with parallel I/O access.
  ! -----------------------------------------------------------
  ! this sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! ----------------------------------
  ! Create the file collectively. (existing files are overwritten)
  ! ---------------------------------
  call H5Fcreate_f ( filename, H5F_ACC_TRUNC_F, file_id, error, &
       access_prp = plist_id)
  ! this closes the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)  
  
  
  call cofitxyz ( uk(:,:,:,1), work)
  call Dump_Field_Backup (work,"ux",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( uk(:,:,:,2), work)
  call Dump_Field_Backup (work,"uy",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( uk(:,:,:,3), work)
  call Dump_Field_Backup (work,"uz",time,dt0,dt1,n1,it,file_id  )
  
  call cofitxyz ( nlk(:,:,:,1,0), work)
  call Dump_Field_Backup (work,"nlkx0",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( nlk(:,:,:,2,0), work)
  call Dump_Field_Backup (work,"nlky0",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( nlk(:,:,:,3,0), work)
  call Dump_Field_Backup (work,"nlkz0",time,dt0,dt1,n1,it,file_id  )
  call cofitxyz ( nlk(:,:,:,1,1), work)
  call Dump_Field_Backup (work,"nlkx1",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( nlk(:,:,:,2,1), work)
  call Dump_Field_Backup (work,"nlky1",time,dt0,dt1,n1,it,file_id  )  
  call cofitxyz ( nlk(:,:,:,3,1), work)
  call Dump_Field_Backup (work,"nlkz1",time,dt0,dt1,n1,it,file_id  )  
  
  ! Close the file.
  call H5Fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library.
  call h5close_f(error)
  
  nbackup = 1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1
  
  
  
  if (mpirank==0) write(*,'("time=",es15.8," dt0=",es15.8)') time, dt0
  
end subroutine Dump_Runtime_Backup




subroutine Dump_Field_Backup (field,dsetname,time,dt0,dt1,n1,it,file_id  )
  use mpi_header
  use share_vars
  use hdf5
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(in) :: field
  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (in) :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname   
  integer,intent(in) :: n1,it

  integer(hid_t), intent(in) :: file_id       ! file identifier 
  integer(hid_t) :: dset_id       ! dataset identifier 
  integer(hid_t) :: filespace     ! dataspace identifier in file 
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier 

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file  
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  
  integer(hssize_t), dimension(rank) :: offset 
  integer(hsize_t),  dimension(rank) :: stride
  integer(hsize_t),  dimension(rank) :: block
  integer :: error  ! error flags
  
  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes
 
  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  ! stride is spacing between elements, this is one here
  stride = 1
  ! how many blocks to select from dataspace
  count =  1
  ! the block contains how many data points to store
  block = dimensions_local
  ! offsets
  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)
  
  
  ! -----------------------------------
  ! Create the data space for the  dataset. 
  ! -----------------------------------
  ! dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)  
  ! Create chunked dataset.
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, dimensions_local, error)
  call H5Dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, & 
       dset_id, error, plist_id)
  call H5Sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error , &
       stride, block)

  ! Create property list for collective dataset write
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively. 
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_file, error, & 
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
       
  ! ------
  ! attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
  ! ------
  adims = (/8/)
  allocate (attributes(1:8))
  aname = "bckp"
  attributes = (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/)
  call write_attribute_dble(adims,aname,attributes,8,dset_id)
  
  ! Close dataspaces, dataset and property list
  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Dclose_f(dset_id, error)
  call H5Pclose_f(plist_id, error)
  
  
  deallocate (attributes)
end subroutine




subroutine Read_Runtime_Backup(time,dt0,dt1,n1,it,uk,nlk,workvis, work)
  use mpi_header 
  use share_vars
  use hdf5
  implicit none
  real(kind=pr),intent(out) :: time,dt1,dt0
  integer,intent(out) :: n1,it
  complex(kind=pr), dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3),intent(out) :: uk
  complex(kind=pr), dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1),intent(out):: nlk
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),intent(out) :: workvis
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),intent(inout) :: work
  character(len=18) :: filename
  
  integer :: error  ! error flags  
  integer(hid_t) :: file_id       ! file identifier 
  integer(hid_t) :: plist_id      ! property list identifier 
  !integer(hid_t) :: dset_id       ! dataset identifier 
  
  if(mpirank==0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! warning: trying to resume a backup file"
  endif

  
  
  filename = 'runtime_backup9.h5'
  
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error) 
  
  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)   
  ! open the file in parallel
  call H5Fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)
  
  
  
  
  call Read_Field_Backup ( work,"ux",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, uk(:,:,:,1) )
  call Read_Field_Backup ( work,"uy",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, uk(:,:,:,2) )
  call Read_Field_Backup ( work,"uz",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, uk(:,:,:,3) )
  
  call Read_Field_Backup ( work,"nlkx0",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,1,0) )
  call Read_Field_Backup ( work,"nlky0",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,2,0) )
  call Read_Field_Backup ( work,"nlkz0",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,3,0) )
  call Read_Field_Backup ( work,"nlkx1",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,1,1) )
  call Read_Field_Backup ( work,"nlky1",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,2,1) )
  call Read_Field_Backup ( work,"nlkz1",time,dt0,dt1,n1,it,file_id )
  call coftxyz ( work, nlk(:,:,:,3,1) )
  
  call H5Fclose_f (file_id,error)
  call H5close_f(error)
 

  call cal_vis ( dt1, workvis )
 
 
  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time, dt0
     write(*,'("!!! DONE READING BACKUP (succes!)")') 
     write(*,'("---------")')
  endif

end subroutine Read_Runtime_Backup




subroutine Read_Field_Backup (field,dsetname,time,dt0,dt1,n1,it,file_id  )
  use mpi_header
  use share_vars
  use hdf5
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(out) :: field
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (out)  :: time,dt1,dt0
  character(len=*), intent (in) :: dsetname   
  integer,intent(out)           :: n1,it
  integer                       :: nx_file,ny_file,nz_file

  integer(hid_t), intent(in) :: file_id       ! file identifier 
  integer(hid_t) :: dset_id       ! dataset identifier 
  integer(hid_t) :: filespace     ! dataspace identifier in file 
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier 

! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file  
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  
  integer(hssize_t), dimension(rank) :: offset 
  integer(hsize_t),  dimension(rank) :: stride
  integer(hsize_t),  dimension(rank) :: block
  integer :: error  ! error flags
  
  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  !integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes
 
  !----------------------------------------------------------------------------
  ! definition of memory distribution
  !----------------------------------------------------------------------------
  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  ! stride is spacing between elements, this is one here
  stride = 1
  ! how many blocks to select from dataspace
  count =  1
  ! the block contains how many data points to store
  block = dimensions_local
  ! offsets
  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3) 

  
  !----------------------------------------------------------------------------
  ! read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)
  
  ! Create chunked dataset
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, dimensions_local, error)
  
  ! Open an existing dataset.
  call H5Dopen_f(file_id, dsetname, dset_id, error)
       
  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, block)

  ! Create property list for collective dataset read
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
 
  call H5Dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )
  
  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error) 
  call H5Pclose_f(plist_id, error) ! note the dataset remains opened
  
  !----------------------------------------------------------------------------
  ! attributes (we save everything in one, all double. to be converted
  ! when reading (to integer)
  !----------------------------------------------------------------------------
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
    stop
  endif
  
  deallocate (attributes)  
  
  
  ! Close dataset
  call H5Dclose_f(dset_id, error)
  
end subroutine



subroutine Dump_Runtime_Backup_OLD(time,dt0,dt1,n1,it,nbackup,uk,nlk,workvis)
  use mpi_header 
  use share_vars
  implicit none
  real(kind=pr),intent(in) :: time,dt1,dt0
  integer,intent(inout) :: n1,nbackup,it
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3),&
       intent(in) :: uk
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1),&
       intent(in):: nlk
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),&
       intent(in) :: workvis
  integer :: filedesc,mpicode
  real(kind=pr) :: t1
  integer,dimension(MPI_STATUS_SIZE) :: mpistatus
  character(len=1) :: name1

  t1=MPI_wtime()

  if(mpirank ==0) then
     write(*,'("*** info: time=",es8.2," dumping runtime_backup",i1," to disk....")') time,nbackup
  endif

  write(name1,'(I1)') nbackup

  ! ---------------------------------------------------------------------------
  ! first part: delete existing file, create new one, and save scalars
  ! ----------------------------------------------------------------------------
  call MPI_FILE_DELETE('runtime_backup'//name1,MPI_INFO_NULL,mpicode)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'runtime_backup'//name1,&
       MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,filedesc,mpicode)
  ! dump time 
  call MPI_FILE_WRITE_ALL(filedesc,time,1,mpireal,mpistatus,mpicode)
  ! dump n1(important when running AB2 scheme!!)
  call MPI_FILE_WRITE_ALL(filedesc,n1,1,mpiinteger,mpistatus,mpicode)
  call MPI_FILE_WRITE_ALL(filedesc,it,1,mpiinteger,mpistatus,mpicode)
  call MPI_FILE_WRITE_ALL(filedesc,nx,1,mpiinteger,mpistatus,mpicode)
  call MPI_FILE_WRITE_ALL(filedesc,ny,1,mpiinteger,mpistatus,mpicode)
  call MPI_FILE_WRITE_ALL(filedesc,nz,1,mpiinteger,mpistatus,mpicode)
  ! dump a few other parameters
  call MPI_FILE_WRITE_ALL(filedesc,(/dt0,dt1/),2,mpireal,mpistatus,mpicode)  
  call MPI_FILE_CLOSE(filedesc,mpicode) 
  ! close file(I really don't yet know why)


  !-----------------------------------------------------------------------------
  ! 2nd part: open file again and read all the fields, one by one.
  ! NOTE: you cannot store uk(:,:,:,1:3) directly, since the ordering
  ! won't match if you use different #cpu for writing and reading
  ! ----------------------------------------------------------------------------

  call MPI_FILE_OPEN(MPI_COMM_WORLD,'runtime_backup'//name1,&
       MPI_MODE_WRONLY+MPI_MODE_APPEND,MPI_INFO_NULL,filedesc,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,uk(:,:,:,1),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,uk(:,:,:,2),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,uk(:,:,:,3),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,1,0),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,2,0),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,3,0),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,1,1),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,2,1),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,nlk(:,:,:,3,1),product(cs),&
       mpicomplex,mpistatus,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,workvis,product(cs),mpireal,&
       mpistatus,mpicode)     
  call MPI_FILE_CLOSE(filedesc,mpicode)


  nbackup=1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1
end subroutine Dump_Runtime_Backup_OLD




subroutine Read_Runtime_Backup_OLD(time,dt0,dt1,n1,it,uk,nlk,workvis)
  use mpi_header 
  use share_vars
  implicit none
  real(kind=pr),intent(out) :: time,dt1,dt0
  integer,intent(out) :: n1,it
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3),&
       intent(out) :: uk
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1),&
       intent(out):: nlk
  real(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),&
       intent(out) :: workvis
  integer :: filedesc,mpicode,ibackup,nx_file,ny_file,nz_file
  real(kind=pr) :: time1
  integer,dimension(MPI_STATUS_SIZE) :: mpistatus
  integer(kind=MPI_OFFSET_KIND) :: mpioffset
  character(len=1) :: name1 
  time=0.d0

  if(mpirank==0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! warning: trying to resume a backup file"
  endif

  ! Read from backup file   
  do ibackup=0,1
     write(name1,'(I1)') ibackup

     call MPI_FILE_OPEN(MPI_COMM_WORLD,'runtime_backup'//name1,&
          MPI_MODE_RDONLY,MPI_INFO_NULL,filedesc,mpicode)

     if(mpicode == 0) then
        ! read time from backup file
        call MPI_FILE_READ_ALL(filedesc,time1,1,mpireal,mpistatus,mpicode)
        !------------------------------------
        ! if the backup is newer, read it
        !------------------------------------
        if(time1 > time) then
           time=time1
           if(mpirank == 0) then
              write(*,'("*** runtime_backup",i1," is at time=",es8.2)') &
                   ibackup,time
           endif

           call MPI_FILE_READ_ALL(filedesc,n1,1,mpiinteger,mpistatus,mpicode)
           call MPI_FILE_READ_ALL(filedesc,it,1,mpiinteger,mpistatus,mpicode)
           call MPI_FILE_READ_ALL(filedesc,nx_file,1,mpiinteger,mpistatus,&
                mpicode)
           call MPI_FILE_READ_ALL(filedesc,ny_file,1,mpiinteger,mpistatus,&
                mpicode)
           call MPI_FILE_READ_ALL(filedesc,nz_file,1,mpiinteger,mpistatus,&
                mpicode)
           call MPI_FILE_READ_ALL(filedesc,dt0,1,mpireal,mpistatus,mpicode)
           call MPI_FILE_READ_ALL(filedesc,dt1,1,mpireal,mpistatus,mpicode)

           if((nx_file/=nx).or.(ny_file/=ny).or.(nz_file/=nz) ) then
              call MPI_FILE_CLOSE(filedesc,mpicode)
              if(mpirank==0) write(*,*) &
                   "resolution of backup file and PARAMS.ini file do not match."
              stop
           endif

           if(mpirank == 0) then
              write(*,'("*** read n1=",i1," it=",i5," dt0=",es12.4," dt1=",es12.4)') n1,it,dt0,dt1
           endif

           call MPI_FILE_GET_POSITION(filedesc,mpioffset,mpicode)
           call MPI_FILE_SET_VIEW(filedesc,mpioffset,MPI_INTEGER,MPI_INTEGER,&
                "native",MPI_INFO_NULL,mpicode)
           !----------------------------
           ! read all the fields
           !----------------------------   
           call MPI_FILE_READ_ORDERED(filedesc,uk(:,:,:,1),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,uk(:,:,:,2),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,uk(:,:,:,3),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,1,0),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,2,0),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,3,0),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,1,1),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,2,1),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,nlk(:,:,:,3,1),product(cs),&
                mpicomplex,mpistatus,mpicode)
           call MPI_FILE_READ_ORDERED(filedesc,workvis,product(cs),mpireal,&
                mpistatus,mpicode)   
        endif
     endif
     call MPI_FILE_CLOSE(filedesc,mpicode)
  enddo


  if(time1 == 0) then
     if(mpirank == 0) then
        write(*,*) 'Unable to resume'
     endif
     stop
  endif

  if(mpirank == 0) then
     write(*,'("!!! DONE READING BACKUP (succes!)")') 
     write(*,'("---------")')
  endif

end subroutine Read_Runtime_Backup_OLD






subroutine save_fields_new(time,dt1,uk,u,vort,nlk,work)
  use mpi_header 
  use share_vars
  implicit none
  real(kind=pr),intent(in) :: time,dt1  
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3),&
       intent(in) :: uk
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3),&
       intent(out):: nlk
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
       intent(inout) :: work
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3),&
       intent(inout) :: vort,u  
  integer :: ix,iy,iz
  character(len=17) :: name
  real(kind=pr) :: t1
  real(kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  

  t1=MPI_wtime()

  !--Set up file name base 
  write(name,'(i5.5)') floor(time*100.d0)

  if(mpirank == 0 ) then 
     write(*,&
          '("*** info: Saving data.... time= ",es8.2,1x," saveflags= ",5(i1))')&
          time,iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask,&
          iSaveSolidVelocity
  endif

  if((iSaveVelocity.ne.0).or.(iSaveVorticity.ne.0).or.(iSavePress.ne.0)) then
       !-----------------------------------------------
       !--Calculate ux and uy in physical space
       !-----------------------------------------------
       call cofitxyz(uk(:,:,:,1),u(:,:,:,1))
       call cofitxyz(uk(:,:,:,2),u(:,:,:,2))
       call cofitxyz(uk(:,:,:,3),u(:,:,:,3))
       !-----------------------------------------------
       !-- SaveVelocity
       !----------------------------------------------- 
       if(iSaveVelocity == 1) then
          call Save_scalar_HDF5 ( time, './fields/ux_'//trim(adjustl(name)), u(:,:,:,1), "ux" )
          call Save_scalar_HDF5 ( time, './fields/uy_'//trim(adjustl(name)), u(:,:,:,2), "uy" )
          call Save_scalar_HDF5 ( time, './fields/uz_'//trim(adjustl(name)), u(:,:,:,3), "uz" )
       endif
  
       if((iSaveVorticity.ne.0).or.(iSavePress.ne.0)) then
          !-----------------------------------------------
          !-- compute vorticity
          !-----------------------------------------------
          do iy=ca(3),cb(3)    ! ky : 0..ny/2-1 ,then,-ny/2..-1     
             ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
             do ix=ca(2),cb(2)  ! kx : 0..nx/2
                kx=scalex*dble(ix)                
                do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then,-nz/2..-1           
                   kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
                   nlk(iz,ix,iy,1)=dcmplx(0d0,1d0)*(ky*uk(iz,ix,iy,3) &
                        - kz*uk(iz,ix,iy,2) )
                   nlk(iz,ix,iy,2)=dcmplx(0d0,1d0)*(kz*uk(iz,ix,iy,1) &
                        - kx*uk(iz,ix,iy,3) )
                   nlk(iz,ix,iy,3)=dcmplx(0d0,1d0)*(kx*uk(iz,ix,iy,2) &
                        - ky*uk(iz,ix,iy,1) )
                enddo
             enddo
          enddo
          ! Transform it to physical space
          call cofitxyz(nlk(:,:,:,1),vort(:,:,:,1)) 
          call cofitxyz(nlk(:,:,:,2),vort(:,:,:,2))
          call cofitxyz(nlk(:,:,:,3),vort(:,:,:,3))
          !-----------------------------------------------
          !-- Save Vorticity
          !----------------------------------------------- 
          if(iSaveVorticity == 1) then
             call Save_scalar_HDF5 ( time, './fields/vorx_'//trim(adjustl(name)), vort(:,:,:,1), "vorx" )
             call Save_scalar_HDF5 ( time, './fields/vory_'//trim(adjustl(name)), vort(:,:,:,2), "vory" )
             call Save_scalar_HDF5 ( time, './fields/vorz_'//trim(adjustl(name)), vort(:,:,:,3), "vorz" )

             ! I don't think we'll keep this for very long:
             work=sqrt(vort(:,:,:,1)**2 + vort(:,:,:,2)**2 +vort(:,:,:,3)**2 )
             call Save_scalar_HDF5 ( time, './fields/vorabs_'//trim(adjustl(name)), work, "vorabs" )
          endif
  
          if(iSavePress == 1) then  
             !-------------------------------------------------------------
             !-- Calculate omega x u(cross-product)
             !-- and transform the result into Fourier space 
             !-------------------------------------------------------------
             if((iPenalization == 1).and.(iMoving==0)) then
                work=u(:,:,:,2)*vort(:,:,:,3)&
                     -u(:,:,:,3)*vort(:,:,:,2)&
                     -u(:,:,:,1)*mask
                call coftxyz(work,nlk(:,:,:,1))
                work=u(:,:,:,3)*vort(:,:,:,1)&
                     -u(:,:,:,1)*vort(:,:,:,3)&
                     -u(:,:,:,2)*mask
                call coftxyz(work,nlk(:,:,:,2))
                work=u(:,:,:,1)*vort(:,:,:,2)&
                     -u(:,:,:,2)*vort(:,:,:,1)&
                     -u(:,:,:,3)*mask
                call coftxyz(work,nlk(:,:,:,3))
             elseif((iPenalization==1).and.(iMoving==1)) then
                work=u(:,:,:,2)*vort(:,:,:,3)&
                     -u(:,:,:,3)*vort(:,:,:,2)&
                     -(u(:,:,:,1)-us(:,:,:,1))*mask
                call coftxyz(work,nlk(:,:,:,1))
                work=u(:,:,:,3)*vort(:,:,:,1)&
                     -u(:,:,:,1)*vort(:,:,:,3)&
                     -(u(:,:,:,2)-us(:,:,:,2))*mask
                call coftxyz(work,nlk(:,:,:,2))
                work=u(:,:,:,1)*vort(:,:,:,2)&
                     -u(:,:,:,2)*vort(:,:,:,1)&
                     -(u(:,:,:,3)-us(:,:,:,3))*mask
                call coftxyz(work,nlk(:,:,:,3))
             else
                work=u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2)
                call coftxyz(work,nlk(:,:,:,1))
                work=u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3)
                call coftxyz(work,nlk(:,:,:,2))
                work=u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1)
                call coftxyz(work,nlk(:,:,:,3))  
             endif
             !-------------------------------------------------------------
             !-- add pressure, new version
             !-- p=(i*kx*sxk + i*ky*syk + i*kz*szk) / k**2
             !-- note: we use rotational formulation: p is NOT the
             !physical pressure
             !-------------------------------------------------------------
             do iy=ca(3),cb(3)  ! ky : 0..ny/2-1 ,then, -ny/2..-1     
                ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
                ky2=ky*ky
                do ix=ca(2),cb(2) ! kx : 0..nx/2
                   kx=scalex*dble(ix)                
                   kx2=kx*kx
                   do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1
                      kz     =scalez*dble(modulo(iz+nz/2,nz)-nz/2)
                      kz2    =kz*kz
                      k_abs_2=kx2+ky2+kz2
                      if(abs(k_abs_2) .ne. 0.0) then  
                         nlk(iz,ix,iy,1)=&
                              (kx*nlk(iz,ix,iy,1)&
                              +ky*nlk(iz,ix,iy,2)&
                              +kz*nlk(iz,ix,iy,3)&
                              )/k_abs_2
                      endif
                   enddo
                enddo
             enddo
             call cofitxyz(nlk(:,:,:,1),work)
             ! work contains total pressure, remove kinetic energy to
             ! get "physical" pressure
             work=work - 0.5d0*(&
                   u(:,:,:,1)*u(:,:,:,1)&
                  +u(:,:,:,2)*u(:,:,:,2)&
                  +u(:,:,:,3)*u(:,:,:,3)&
                  )
             call Save_scalar_HDF5(time, './fields/p_'//trim(adjustl(name)), work, "p" )
          endif
       endif
    endif

  !-----------------------------------------------
  !-- Save Mask
  !----------------------------------------------- 

  if((iSaveMask==1).and.(iPenalization==1)) then
    call Save_scalar_HDF5(time,'./fields/mask_'//trim(adjustl(name)),mask, "mask" )
  endif
  
  if((iSaveSolidVelocity==1).and.(iPenalization==1).and.(iMoving==1)) then
    call Save_scalar_HDF5(time,'./fields/usx_'//trim(adjustl(name)),us(:,:,:,1), "usx" )
    call Save_scalar_HDF5(time,'./fields/usy_'//trim(adjustl(name)),us(:,:,:,2), "usy" )
    call Save_scalar_HDF5(time,'./fields/usz_'//trim(adjustl(name)),us(:,:,:,3), "usz" )
  endif


  time_save=time_save + MPI_wtime() - t1
end subroutine save_fields_new

subroutine write_attribute_dble(adims,aname,attribute,dim,dset_id)
  use mpi_header
  use share_vars
  use HDF5
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  real (kind=pr), DIMENSION(dim), intent (in) :: attribute
  character(len=4), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id       ! dataset identifier 
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier

  ! Determine the dataspace identifier ASPACE_ID
  call h5screate_simple_f(arank,adims,aspace_id,error)

  ! set attr_id, ie create an attribute attached to the object DSET_ID
  call h5acreate_f(dset_id,aname,H5T_NATIVE_DOUBLE,aspace_id,attr_id,error)
  
  ! Write the attribute data ATTRIBUTE to the attribute identifier ATTR_ID.
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attribute, adims, error)  

  call h5aclose_f(attr_id, error) ! Close the attribute.
  call h5sclose_f(aspace_id, error) ! Terminate access to the data space.
  
end subroutine write_attribute_dble

subroutine write_attribute_int(adims,aname,attribute,dim,dset_id)
  use mpi_header
  use share_vars
  use HDF5
  implicit none

  integer, intent(in) :: dim
  integer(hsize_t), DIMENSION(dim), intent(in) :: adims  ! Attribute dimension
  integer, DIMENSION(dim), intent (in) :: attribute
  character(len=4), intent(in) :: aname ! attribute name
  integer(hid_t),intent(in) :: dset_id       ! dataset identifier 
  integer, parameter :: arank = 1

  integer :: error  ! error flags
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier

  ! Determine the dataspace identifier aspace_id
  call h5screate_simple_f(arank, adims, aspace_id, error)
  ! set attr_id, ie create an attribute attached to the object dset_id
  call h5acreate_f(dset_id,aname,H5T_NATIVE_INTEGER,aspace_id,attr_id,error)
  ! Write the attribute data ATTRIBUTE to the attribute identifier attr_id.
  call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,attribute,adims,error)  
  call h5aclose_f(attr_id,error) ! Close the attribute.
  call h5sclose_f(aspace_id,error) ! Terminate access to the data space.
  
end subroutine write_attribute_int
