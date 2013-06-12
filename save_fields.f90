subroutine Save_scalar_HDF5 ( time,  filename, field_out, dsetname )
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

  integer(hsize_t), dimension(rank) :: dimensions_file  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  
  integer(hssize_t), dimension(rank) :: offset 
  integer(hsize_t),  dimension(rank) :: stride
  integer(hsize_t),  dimension(rank) :: block
  integer :: error, error_n  ! error flags
  
  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  INTEGER(HID_T) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  character(len=11) :: aname2="domain_size" ! attribute name
  !----------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this first part is to tell HDF5 how our (real!) data is organized
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  dimensions_file = (/nx,ny,nz/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1

  ! stride is spacing between elements, this is one here
  stride(1) = 1 
  stride(2) = 1 
  stride(3) = 1

  ! how many blocks to select from dataspace
  count(1) =  1 
  count(2) =  1 
  count(3) =  1 

  ! the block contains how many data points to store
  block(1) = dimensions_local(1)
  block(2) = dimensions_local(2)
  block(3) = dimensions_local(3)

  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)

  write (*,'("mpirank=",i1," dims = ("i2,",",i2,",",i2,"), offset=(",i2,",",i2,",",i2,")"  )') &
       mpirank, dimensions_local(1),dimensions_local(2),dimensions_local(3),offset(1),offset(2),offset(3)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! now we call the HDF subroutines from the "chunks" example
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

  !
  ! Initialize HDF5 library and Fortran interfaces.
  !
  CALL h5open_f(error) 

  ! -----------------------------------------------------------
  ! Setup file access property list with parallel I/O access.
  ! -----------------------------------------------------------
  ! this sets up a property list ("plist_id") with standard values for FILE_ACCESS
  CALL H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  CALL H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! ----------------------------------
  ! Create the file collectively.
  ! ---------------------------------
  CALL H5Fcreate_f ( trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  ! this closes the property list (we'll re-use it)
  CALL H5Pclose_f(plist_id, error)

  ! -----------------------------------
  ! Create the data space for the  dataset. 
  ! -----------------------------------
  ! dataspace in the file: contains all data from all procs
  CALL H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  CALL H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset.
  CALL H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  CALL H5Pset_chunk_f(plist_id, rank, dimensions_local, error)
  CALL H5Dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, & ! double precision in memory...
       dset_id, error, plist_id)
  CALL H5Sclose_f(filespace, error)


  ! Select hyperslab in the file.
  CALL H5Dget_space_f(dset_id, filespace, error)
  CALL H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error , &
       stride, block)

  ! Create property list for collective dataset write
  CALL H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  CALL H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


  ! Write the dataset collectively. 
  CALL H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, error, & ! but single precision on disk..
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

       
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
  ! now the data is written, we take care of the ATTRIBUTE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ------
  ! time
  ! ------
  ! Create scalar data space for the attribute.
  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error)  
  ! Create dataset attribute.
  aname = "time"
  CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)
  ! Write the attribute data.
  CALL h5awrite_f(attr_id, atype_id, time, adims, error)  ! here the value of the attribute is "time"
  ! Close the attribute.
  CALL h5aclose_f(attr_id, error)
  ! Terminate access to the data space.
  CALL h5sclose_f(aspace_id, error)  
  ! ------
  ! eps
  ! ------
  ! Create scalar data space for the attribute.
  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error) 
  ! Create dataset attribute.
  aname = "epsi"
  CALL h5acreate_f(dset_id, aname, atype_id, aspace_id, attr_id, error)
  ! Write the attribute data.
  CALL h5awrite_f(attr_id, atype_id, eps, adims, error)  ! here the value of the attribute is "time"
  ! Close the attribute.
  CALL h5aclose_f(attr_id, error)
  ! Terminate access to the data space.
  CALL h5sclose_f(aspace_id, error)  
  ! ------
  ! domain size
  ! ------
  adims = (/3/)
  ! Create scalar data space for the attribute.
  CALL h5screate_simple_f(arank, adims, aspace_id, error)
  ! Create datatype for the attribute.
  CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype_id, error) 
  ! Create dataset attribute.
  aname2 = "domain_size"
  CALL h5acreate_f(dset_id, aname2, atype_id, aspace_id, attr_id, error)
  ! Write the attribute data.
  CALL h5awrite_f(attr_id, atype_id, (/xl, yl, zl/), adims, error)  ! here the value of the attribute is "time"
  ! Close the attribute.
  CALL h5aclose_f(attr_id, error)
  ! Terminate access to the data space.
  CALL h5sclose_f(aspace_id, error)   
  
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  ! time to close everything
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  
  

  
  ! Close dataspaces.
  CALL H5Sclose_f(filespace, error)
  CALL H5Sclose_f(memspace, error)
  ! Close the dataset.
  CALL H5Dclose_f(dset_id, error)
  ! Close the property list.
  CALL H5Pclose_f(plist_id, error)
  ! Close the file.
  CALL H5Fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library.
  CALL h5close_f(error)

  if (mpirank==0) then
     call WriteXML ( time, trim(adjustl(filename)) , trim(adjustl(dsetname)) )
  endif

end subroutine Save_scalar_HDF5




subroutine Save_vector_HDF5 ( time,  filename, field_out, dsetname )
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  use HDF5
  implicit none
  real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent(in) :: field_out
  integer, parameter :: rank = 4 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename   ! file name
  character(len=*), intent (in) :: dsetname 

  integer(hid_t) :: file_id       ! file identifier 
  integer(hid_t) :: dset_id       ! dataset identifier 
  integer(hid_t) :: filespace     ! dataspace identifier in file 
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier 

  integer(hsize_t), dimension(rank) :: dimensions_file  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  
  integer(hssize_t), dimension(rank) :: offset 
  integer(hsize_t),  dimension(rank) :: stride
  integer(hsize_t),  dimension(rank) :: block
  integer :: error, error_n  ! error flags

  !----------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this first part is to tell HDF5 how our (real!) data is organized
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  dimensions_file = (/nx,ny,nz,3/)
  dimensions_local(1) = rb(1)-ra(1) +1
  dimensions_local(2) = rb(2)-ra(2) +1
  dimensions_local(3) = rb(3)-ra(3) +1
  dimensions_local(4) = 3

  ! stride is spacing between elements, this is one here
  stride = 1
  ! how many blocks to select from dataspace
  count =  1 
  
  ! the block contains how many data points to store
  block = dimensions_local

  offset(1) = ra(1)
  offset(2) = ra(2)
  offset(3) = ra(3)
  offset(4) = 0


  ! Initialize HDF5 library and Fortran interfaces.
  CALL h5open_f(error) 

  ! -----------------------------------------------------------
  ! Setup file access property list with parallel I/O access.
  ! -----------------------------------------------------------
  ! this sets up a property list ("plist_id") with standard values for FILE_ACCESS
  CALL H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  CALL H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  !
  ! Create the file collectively.
  ! 
  CALL H5Fcreate_f ( trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
  ! this closes the property list (we'll re-use it)
  CALL H5Pclose_f(plist_id, error)
  !
  ! Create the data space for the  dataset. 
  !
  ! dataspace in the file: contains all data from all procs
  CALL H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  CALL H5Screate_simple_f(rank, dimensions_local, memspace, error)

  !
  ! Create chunked dataset.
  !
  CALL H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  CALL H5Pset_chunk_f(plist_id, rank, dimensions_local, error)
  CALL H5Dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, & ! double precision in memory...
       dset_id, error, plist_id)
  CALL H5Sclose_f(filespace, error)



  ! 
  ! Select hyperslab in the file.
  !
  CALL H5Dget_space_f(dset_id, filespace, error)
  CALL H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, error , &
       stride, block)


  !
  ! Create property list for collective dataset write
  !
  CALL H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
  CALL H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)


  ! Write the dataset collectively. 
  CALL H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, error, & ! but single precision on disk..
       file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)


  ! Close dataspaces.
  CALL H5Sclose_f(filespace, error)
  CALL H5Sclose_f(memspace, error)
  ! Close the dataset.
  CALL H5Dclose_f(dset_id, error)
  ! Close the property list.
  CALL H5Pclose_f(plist_id, error)
  ! Close the file.
  CALL H5Fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library.
  CALL h5close_f(error)

!   if (mpirank==0) then
!      call WriteXML ( time, trim(adjustl(filename)) , trim(adjustl(dsetname)) )
!   endif

end subroutine Save_vector_HDF5




subroutine WriteXML ( time, filename, dsetname )
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

end subroutine WriteXML





subroutine Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,workvis)
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
  real(kind=pr) :: t1,tmp,tmp_local
  integer,dimension(MPI_STATUS_SIZE) :: mpistatus
  character(len=17) :: name
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
end subroutine Dump_Runtime_Backup




subroutine Read_Runtime_Backup(time,dt0,dt1,n1,it,uk,nlk,workvis)
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
  character(len=17) :: name
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

end subroutine Read_Runtime_Backup






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
  integer :: filedesc,mpicode
  integer,dimension(MPI_STATUS_SIZE) :: mpistatus
  character(len=17) :: name
  character(len=1) :: name1
  real(kind=pr) :: u_max_w,divu,divumax,t1=0.d0,t2=0.d0,t3=0.d0,t4=0.d0,&
       t5=0.d0,t6=0.d0
  real(kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  real(kind=pr),dimension(3) :: u_max,u_loc
  complex(kind=pr) :: qk

  ! -------------------------------------------------------
  ! - interface for SaveFile subroutine
  ! -------------------------------------------------------
  interface                                                                
     subroutine SaveFile(filename,field_out)
       use mpi_header ! Module incapsulates mpif.
       use share_vars
       implicit none
       integer,parameter :: pr_out=8   ! double precision array for output
       integer,parameter :: mpireal_out=MPI_DOUBLE_PRECISION 
       ! double precision array for output
       character(len=*),intent(in) :: filename
       real(kind=pr_out),dimension(:,:,:),intent(in) :: field_out
     end subroutine SaveFile

     subroutine SaveVTK(fname,u,vort,p) 
       use share_vars
       implicit none
       real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
            intent(in) :: p
       real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3),&
            intent(in) :: vort,u 
       character(len=*),intent(in) :: fname
     end subroutine SaveVTK

  end interface

  t1=MPI_wtime()

  !--Set up file name base 
  write(name,'(i5.5)') floor(time*100.d0)

  if(mpirank == 0 ) then 
     write(*,&
          '("*** info: Saving data.... time= ",es8.2,1x," saveflags= ",5(i1))')&
          time,iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask,&
          iSaveSolidVelocity
  endif

  !   if((iSaveVelocity.ne.0).or.(iSaveVorticity.ne.0).or.(iSavePress.ne.0)) then
  !      !-----------------------------------------------
  !      !--Calculate ux and uy in physical space
  !      !-----------------------------------------------
  !      call cofitxyz(uk(:,:,:,1),u(:,:,:,1))
  !      call cofitxyz(uk(:,:,:,2),u(:,:,:,2))
  !      call cofitxyz(uk(:,:,:,3),u(:,:,:,3))
  !      !-----------------------------------------------
  !      !-- SaveVelocity
  !      !----------------------------------------------- 
  !      if(iSaveVelocity == 1) then
  !         call SaveFile('ux_'//trim(adjustl(name))//'.mpiio' ,u(:,:,:,1) )
  !         call SaveFile('uy_'//trim(adjustl(name))//'.mpiio' ,u(:,:,:,2) )
  !         call SaveFile('uz_'//trim(adjustl(name))//'.mpiio' ,u(:,:,:,3) )
  !      endif
  ! 
  !      if((iSaveVorticity.ne.0).or.(iSavePress.ne.0)) then
  !         !-----------------------------------------------
  !         !-- compute vorticity
  !         !-----------------------------------------------
  !         do iy=ca(3),cb(3)    ! ky : 0..ny/2-1 ,then,-ny/2..-1     
  !            ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
  !            do ix=ca(2),cb(2)  ! kx : 0..nx/2
  !               kx=scalex*dble(ix)                
  !               do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then,-nz/2..-1           
  !                  kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
  !                  nlk(iz,ix,iy,1)=dcmplx(0d0,1d0)*(ky*uk(iz,ix,iy,3) &
  !                       - kz*uk(iz,ix,iy,2) )
  !                  nlk(iz,ix,iy,2)=dcmplx(0d0,1d0)*(kz*uk(iz,ix,iy,1) &
  !                       - kx*uk(iz,ix,iy,3) )
  !                  nlk(iz,ix,iy,3)=dcmplx(0d0,1d0)*(kx*uk(iz,ix,iy,2) &
  !                       - ky*uk(iz,ix,iy,1) )
  !               enddo
  !            enddo
  !         enddo
  !         ! Transform it to physical space
  !         call cofitxyz(nlk(:,:,:,1),vort(:,:,:,1)) 
  !         call cofitxyz(nlk(:,:,:,2),vort(:,:,:,2))
  !         call cofitxyz(nlk(:,:,:,3),vort(:,:,:,3))
  !         !-----------------------------------------------
  !         !-- Save Vorticity
  !         !----------------------------------------------- 
  !         if(iSaveVorticity == 1) then
  !            call SaveFile('vorx_'//trim(adjustl(name))//'.mpiio',vort(:,:,:,1) )
  !            call SaveFile('vory_'//trim(adjustl(name))//'.mpiio',vort(:,:,:,2) )
  !            call SaveFile('vorz_'//trim(adjustl(name))//'.mpiio',vort(:,:,:,3) )
  !            work=sqrt(vort(:,:,:,1)**2 + vort(:,:,:,2)**2 +vort(:,:,:,3)**2 )
  !            call SaveFile('vorabs_'//trim(adjustl(name))//'.mpiio',work )
  !         endif
  ! 
  !         if(iSavePress == 1) then  
  !            !-------------------------------------------------------------
  !            !-- Calculate omega x u(cross-product)
  !            !-- and transform the result into Fourier space 
  !            !-------------------------------------------------------------
  !            if((iPenalization == 1).and.(iMoving==0)) then
  !               work=u(:,:,:,2)*vort(:,:,:,3)&
  !                    -u(:,:,:,3)*vort(:,:,:,2)&
  !                    -u(:,:,:,1)*mask
  !               call coftxyz(work,nlk(:,:,:,1))
  !               work=u(:,:,:,3)*vort(:,:,:,1)&
  !                    -u(:,:,:,1)*vort(:,:,:,3)&
  !                    -u(:,:,:,2)*mask
  !               call coftxyz(work,nlk(:,:,:,2))
  !               work=u(:,:,:,1)*vort(:,:,:,2)&
  !                    -u(:,:,:,2)*vort(:,:,:,1)&
  !                    -u(:,:,:,3)*mask
  !               call coftxyz(work,nlk(:,:,:,3))
  !            elseif((iPenalization==1).and.(iMoving==1)) then
  !               work=u(:,:,:,2)*vort(:,:,:,3)&
  !                    -u(:,:,:,3)*vort(:,:,:,2)&
  !                    -(u(:,:,:,1)-us(:,:,:,1))*mask
  !               call coftxyz(work,nlk(:,:,:,1))
  !               work=u(:,:,:,3)*vort(:,:,:,1)&
  !                    -u(:,:,:,1)*vort(:,:,:,3)&
  !                    -(u(:,:,:,2)-us(:,:,:,2))*mask
  !               call coftxyz(work,nlk(:,:,:,2))
  !               work=u(:,:,:,1)*vort(:,:,:,2)&
  !                    -u(:,:,:,2)*vort(:,:,:,1)&
  !                    -(u(:,:,:,3)-us(:,:,:,3))*mask
  !               call coftxyz(work,nlk(:,:,:,3))
  !            else
  !               work=u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2)
  !               call coftxyz(work,nlk(:,:,:,1))
  !               work=u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3)
  !               call coftxyz(work,nlk(:,:,:,2))
  !               work=u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1)
  !               call coftxyz(work,nlk(:,:,:,3))  
  !            endif
  !            !-------------------------------------------------------------
  !            !-- add pressure, new version
  !            !-- p=(i*kx*sxk + i*ky*syk + i*kz*szk) / k**2
  !            !-- note: we use rotational formulation: p is NOT the
  !            !physical pressure
  !            !-------------------------------------------------------------
  !            do iy=ca(3),cb(3)  ! ky : 0..ny/2-1 ,then, -ny/2..-1     
  !               ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
  !               ky2=ky*ky
  !               do ix=ca(2),cb(2) ! kx : 0..nx/2
  !                  kx=scalex*dble(ix)                
  !                  kx2=kx*kx
  !                  do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1           
  !                     kz     =scalez*dble(modulo(iz+nz/2,nz)-nz/2)
  !                     kz2    =kz*kz
  !                     k_abs_2=kx2+ky2+kz2
  !                     if(abs(k_abs_2) .ne. 0.0) then  
  !                        nlk(iz,ix,iy,1)=&
  !                             (kx*nlk(iz,ix,iy,1)&
  !                             +ky*nlk(iz,ix,iy,2)&
  !                             +kz*nlk(iz,ix,iy,3)&
  !                             )/k_abs_2
  !                     endif
  !                  enddo
  !               enddo
  !            enddo
  !            call cofitxyz(nlk(:,:,:,1),work)
  !            ! work contains total pressure, remove kinetic energy to
  !            ! get "physical" pressure
  !            work=work - 0.5d0*(&
  !                 u(:,:,:,1)*u(:,:,:,1)&
  !                 +u(:,:,:,2)*u(:,:,:,2)&
  !                 +u(:,:,:,3)*u(:,:,:,3)&
  !                 )
  !            call SaveFile('p_'//trim(adjustl(name))//'.mpiio',work(:,:,:) )
  ! 
  !         endif
  !      endif
  !   endif

  !-----------------------------------------------
  !-- Save Mask
  !----------------------------------------------- 

  if((iSaveMask==1).and.(iPenalization==1)) then
     call Save_scalar_HDF5(time,'mask_'//trim(adjustl(name)),mask, "mask" )
     
     call Save_scalar_HDF5(time,'usx_'//trim(adjustl(name)),us(:,:,:,1), "usx" )
     call Save_scalar_HDF5(time,'usy_'//trim(adjustl(name)),us(:,:,:,2), "usy" )
     call Save_scalar_HDF5(time,'usz_'//trim(adjustl(name)),us(:,:,:,3), "usz" )
     
     
     call Save_vector_HDF5(time,'b_'//trim(adjustl(name)),us(:,:,:,:), "b" )
  endif
  !   if((iSaveSolidVelocity==1).and.(iPenalization==1).and.(iMoving==1)) then
  !      call SaveFile('usx_'//trim(adjustl(name))//'.mpiio',us(:,:,:,1) )
  !      call SaveFile('usy_'//trim(adjustl(name))//'.mpiio',us(:,:,:,2) )
  !      call SaveFile('usz_'//trim(adjustl(name))//'.mpiio',us(:,:,:,3) )
  !   endif






  !------------------------------------------------
  ! TEMP::: compute divergence
  !-----------------------------------------------
  ! compute max val of {|div(.)|/|.|} over entire domain
  !   do iz=ca(1),cb(1)
  !     kz=scalez*(modulo(iz+nz/2,nz) -nz/2)
  !     do iy=ca(3),cb(3)
  !  ky=scaley*(modulo(iy+ny/2,ny) -ny/2)
  !  do ix=ca(2),cb(2)
  !    kx=scalex*ix
  !    ! divergence of velocity field
  !    nlk(iz,ix,iy,1)=dcmplx(0.d0,1.d0)*(kx*uk(iz,ix,iy,1)+ky*uk(iz,ix,iy,2)+kz*uk(iz,ix,iy,3))
  !  enddo
  !     enddo
  !   enddo
  !   ! now nlk(:,:,:,1) contains divergence field
  !   call cofitxyz(nlk(:,:,:,1),work)
  ! 
  !   write(name,'(i5.5)') floor(time*100.d0)
  !   call SaveFile ('divu_'//trim(adjustl(name))//'.mpiio' , work )


  time_save=time_save + MPI_wtime() - t1
end subroutine save_fields_new







subroutine SaveFile(filename,field_out)
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  integer,parameter :: pr_out=8   ! double precision array for output
  integer,parameter :: mpireal_out=MPI_DOUBLE_PRECISION 
  ! double precision array for output
  character(len=*),intent(in) :: filename
  real(kind=pr_out),dimension(:,:,:),intent(in) :: field_out
  integer :: filedesc,mpicode
  integer,dimension(MPI_STATUS_SIZE) :: mpistatus

  ! modified: automatically stores in subfolder fields
  call MPI_FILE_DELETE('./fields/'//filename,MPI_INFO_NULL,mpicode)
  call MPI_FILE_OPEN(MPI_COMM_WORLD,'./fields/'//filename,&
       MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,filedesc,mpicode)
  call MPI_FILE_WRITE_ORDERED(filedesc,field_out,product(rs),mpireal_out,&
       mpistatus,mpicode)
  call MPI_FILE_CLOSE(filedesc,mpicode)    

end subroutine SaveFile

