! Wrapper for saving fields routine
subroutine save_fields_new(time,uk,u,vort,nlk,work)
  use vars
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  select case(method(1:3))
     case("fsi") 
        call save_fields_new_fsi(time,uk,u,vort,nlk,work)
     case("mhd") 
        call save_fields_new_mhd(time,uk,u,vort,nlk,work)
     case default
        if (mpirank == 0) write(*,*) "Error! Unkonwn method in save_fields_new"
        stop
  end select
end subroutine save_fields_new


! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files. 
! The latest version calls cal_nlk_fsi to avoid redudant code. 
! note cal_nlk_fsi returns the NL+penal term in phys space in the work
! array "vort" which is why we have to recompute the vorticity
subroutine save_fields_new_fsi(time,uk,u,vort,nlk,work)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  character(len=17) :: name
  integer :: i

  !--Set up file name base    
  if ( save_only_one_period == "yes" ) then
    ! overwrite files from last period to save disk space
    ! i.e. t=1.05 is written to t=0.05, as well as 2.05 and 3.05
    write(name,'(i5.5)') floor( (time-real(floor(time)))*100.d0 )
  else
    ! name is just the time
    write(name,'(i5.5)') floor(time*100.d0) 
  endif
  
  name=trim(adjustl(name))

  if (mpirank == 0 ) then
    write(*,'("Saving data, time= ",es8.2,1x," flags= ",5(i1)," str=",A)') & 
    time, &
    iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask,iSaveSolidVelocity, &
    trim(adjustl(name))
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cal_nlk_fsi (time,0,nlk,uk,u,vort,work) 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !-------------  
  ! Velocity
  !-------------
  if (iSaveVelocity == 1) then
    call Save_Field_HDF5(time,"./ux_"//name,u(:,:,:,1),"ux")
    call Save_Field_HDF5(time,"./uy_"//name,u(:,:,:,2),"uy")
    call Save_Field_HDF5(time,"./uz_"//name,u(:,:,:,3),"uz")
  endif

  !-------------  
  ! Pressure
  !-------------
  if (iSavePress == 1) then  
    call fft3(nlk,vort) ! nlk now NL+penal term in Fourier space
    ! store the total pressure in the work array
    call compute_pressure(nlk(:,:,:,1),nlk)
    ! total pressure in phys-space
    call ifft(work,nlk(:,:,:,1))
    ! get actuall pressure (we're in the rotational formulation)
    work = work - 0.5d0*( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
    call Save_Field_HDF5(time,'./p_'//name,work,"p")
  endif
     
  !-------------  
  ! Vorticity
  !-------------   
  if (iSaveVorticity==1) then
    ! cal_nlk overwrote the vorticity, recompute it
    call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),& 
               uk(:,:,:,1), uk(:,:,:,2), uk(:,:,:,3)) 
    call ifft3(vort,nlk)      
    !-- Save Vorticity
    if (iSaveVorticity == 1) then
      call Save_Field_HDF5(time,"./vorx_"//name,vort(:,:,:,1),"vorx")
      call Save_Field_HDF5(time,"./vory_"//name,vort(:,:,:,2),"vory")
      call Save_Field_HDF5(time,"./vorz_"//name,vort(:,:,:,3),"vorz")
    endif
  endif
      
  !-------------  
  ! Mask
  !-------------
  if (iSaveMask == 1 .and. iPenalization == 1) then
    call Save_Field_HDF5(time,'./mask_'//name,eps*mask,"mask")
  endif
  
  !-------------  
  ! solid velocity
  !-------------
  if (iSaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call Save_Field_HDF5(time,'./usx_'//name,us(:,:,:,1),"usx")
    call Save_Field_HDF5(time,'./usy_'//name,us(:,:,:,2),"usy")
    call Save_Field_HDF5(time,'./usz_'//name,us(:,:,:,3),"usz")
  endif
  
  if (mpirank==0) write(*,*) "Done saving."
end subroutine save_fields_new_fsi


! Write the field field_out to file filename, saving the name of the
! field, dsetname, (as well as time, eps, and resolution) to the
! metadata of the HDF file .
subroutine Save_Field_HDF5(time,filename,field_out,dsetname)
  use mpi_header
  use vars
  use HDF5
  implicit none
  
  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr), intent (in) :: time
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
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store the MPI IO comminucator
  ! information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
  call H5Fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
       file_id, error, access_prp = plist_id)
  ! this closes the property list plist_id (we'll re-use it)
  call H5Pclose_f(plist_id, error)

  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset.
  ! NB: chunking and hyperslab are unrelated
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, chunking_dims, error)
  ! Output files are single-precition
  call H5Dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, filespace, &
       dset_id, error, plist_id)
  call H5Sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error, stride, dimensions_local)

  ! Create property list for collective dataset write
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively, double precision in memory
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
       error, file_space_id = filespace, mem_space_id = memspace,&
       xfer_prp = plist_id)

  !!! Write the attributes to the HDF files.
  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
  adims = (/1/)
  call write_attribute_dble(adims,"time",(/time/),1,dset_id)
  call write_attribute_dble(adims,"epsi",(/eps/),1,dset_id)
  adims = (/3/)
  call write_attribute_dble(adims,"domain_size",(/xl,yl,zl/),3,dset_id)
  call write_attribute_int(adims,"nxyz",(/nx,ny,nz/),3,dset_id)

  !!! Close dataspaces:
  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Dclose_f(dset_id, error) ! Close the dataset.
  call H5Pclose_f(plist_id, error) ! Close the property list.
  call H5Fclose_f(file_id, error) ! Close the file.
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.

  ! write the XMF data for all of the saved fields
  if ((mpirank==0).and.(iSaveXMF==1)) then
     ! the filename contains a leading "./" which we must remove
     call Write_XMF(time,&
          trim(adjustl(filename(3:len(filename)))),&
          trim(adjustl(dsetname))&
          )
  endif

  time_save=time_save + MPI_wtime() - t1 ! performance analysis
end subroutine Save_Field_HDF5


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


! Write the restart file. nlk(...,0) and nlk(...,1) are saved, the
! time steps, and what else? FIXME: document what is saved.
subroutine Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
  use mpi_header
  use vars
  use hdf5
  implicit none

  real(kind=pr),intent(inout) :: time,dt1,dt0
  integer,intent(inout) :: n1,nbackup,it
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(in)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  character(len=18) :: filename

  real(kind=pr) :: t1
  integer :: error  ! error flags
  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: plist_id      ! property list identifier

  t1=MPI_wtime() ! performance diagnostic

  if(mpirank == 0) then
     write(*,'("*** info: time=",es8.2," dumping runtime_backup",i1,".h5 to disk....")') time, nbackup
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".h5")') nbackup

  ! Initialize HDF5 library and Fortran interfaces:
  call h5open_f(error)

  !!! Setup file access property list with parallel I/O access.
  ! Set up a property list ("plist_id") with standard values for
  ! FILE_ACCESS:
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! Modify the property list and store MPI IO comminucator information
  ! in the file access property list:
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

  ! Create the file collectively. (existing files are overwritten)
  call H5Fcreate_f ( filename, H5F_ACC_TRUNC_F, file_id, error, &
       access_prp = plist_id)
  ! Close the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)

  ! Write the fluid backup field:
  call ifft(work,uk(:,:,:,1))
  call Dump_Field_Backup (work,"ux",time,dt0,dt1,n1,it,file_id)
  call ifft(work,uk(:,:,:,2))
  call Dump_Field_Backup (work,"uy",time,dt0,dt1,n1,it,file_id)
  call ifft(work,uk(:,:,:,3))
  call Dump_Field_Backup (work,"uz",time,dt0,dt1,n1,it,file_id)

  if(method == "mhd") then
     ! Write the MHD backup field:
     call ifft(work,uk(:,:,:,4))
     call Dump_Field_Backup (work,"bx",time,dt0,dt1,n1,it,file_id)
     call ifft(work,uk(:,:,:,5))
     call Dump_Field_Backup (work,"by",time,dt0,dt1,n1,it,file_id)
     call ifft(work,uk(:,:,:,6))
     call Dump_Field_Backup (work,"bz",time,dt0,dt1,n1,it,file_id)
  endif

  ! Write the fluid nonlinear term backup:
  call ifft(work,nlk(:,:,:,1,0))
  call Dump_Field_Backup (work,"nlkx0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,2,0))
  call Dump_Field_Backup (work,"nlky0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,3,0))
  call Dump_Field_Backup (work,"nlkz0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,1,1))
  call Dump_Field_Backup (work,"nlkx1",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,2,1))
  call Dump_Field_Backup (work,"nlky1",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,3,1))
  call Dump_Field_Backup (work,"nlkz1",time,dt0,dt1,n1,it,file_id)
  
  if(method == "mhd") then
     ! Write the MHD backup field:
     call ifft(work,nlk(:,:,:,4,0))
     call Dump_Field_Backup (work,"bnlkx0",time,dt0,dt1,n1,it,file_id)
     call ifft(work,nlk(:,:,:,5,0))
     call Dump_Field_Backup (work,"bnlky0",time,dt0,dt1,n1,it,file_id)
     call ifft(work,nlk(:,:,:,6,0))
     call Dump_Field_Backup (work,"bnlkz0",time,dt0,dt1,n1,it,file_id)
     call ifft(work,nlk(:,:,:,4,1))
     call Dump_Field_Backup (work,"bnlkx1",time,dt0,dt1,n1,it,file_id)
     call ifft(work,nlk(:,:,:,5,1))
     call Dump_Field_Backup (work,"bnlky1",time,dt0,dt1,n1,it,file_id)
     call ifft(work,nlk(:,:,:,6,1))
     call Dump_Field_Backup (work,"bnlkz1",time,dt0,dt1,n1,it,file_id)
  endif

  ! Close the file:
  call H5Fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library:
  call h5close_f(error)

  nbackup = 1 - nbackup
  time_bckp=time_bckp + MPI_wtime() -t1 ! Performance diagnostic

  if(mpirank ==0 ) write(*,'("<<< info: done saving backup.")')
end subroutine Dump_Runtime_Backup


! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "file_id". Attributes are stores in one attribute
! "bckp" which contains 8 values
subroutine Dump_Field_Backup (field,dsetname,time,dt0,dt1,n1,it,file_id  )
  use mpi_header
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
     call MPI_REDUCE ( dimensions_local(i), chunking_dims(i),1, &
          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST  ( chunking_dims(i), 1, MPI_INTEGER8, 0, &
          MPI_COMM_WORLD, mpierror )
  enddo

  ! -----------------------------------
  ! Create the data space for the  dataset.
  ! -----------------------------------
  ! Dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)
  ! Create chunked dataset.
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, chunking_dims, error)
  call H5Dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, filespace, &
       dset_id, error, plist_id)
  call H5Sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset write
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  ! Write the dataset collectively.
  call H5Dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field, dimensions_file, error, &
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
  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Dclose_f(dset_id, error)
  call H5Pclose_f(plist_id, error)

  deallocate(attributes)
end subroutine Dump_Field_Backup





! Read in a single file that follows the naming convention
! this is a serial routine (parallel version below)
! note you need to know what dimension the file has,
! call fetch_attributes first
subroutine Read_Single_File_serial ( filename, field )
  use mpi_header
  use vars
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr), intent (out) :: field(0:nx-1,0:ny-1,0:nz-1)
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time,dt1,dt0, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: n1,it
  integer                       :: nx_file,ny_file,nz_file, mpierror, i  

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
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  !integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes  
  
  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )

  if (mpisize>1) then
    write (*,*) "this routine is currently serial only"
    stop
  endif
  
  ! check if file exist
  call check_file_exists( filename )

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
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call H5Dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call H5Dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
  call H5Dclose_f(dset_id, error)
  call H5Fclose_f(file_id,error)
  call H5close_f(error)
  
  
end subroutine Read_Single_File_serial




! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine Read_Single_File ( filename, field )
  use mpi_header
  use vars
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),&
  dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
  intent (out) :: field
  
  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time,dt1,dt0, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: n1,it
  integer                       :: nx_file,ny_file,nz_file, mpierror, i  

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
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension
  integer(hid_t) :: aspace_id     ! Attribute Dataspace identifier
  !integer(hid_t) :: atype_id      ! Attribute Dataspace identifier
  integer(hid_t) :: attr_id       ! Attribute identifier
  character(len=4) :: aname ! attribute name
  real (kind=pr), dimension (:), allocatable :: attributes  
  logical :: exist1

  ! the dataset is named the same way as the file: (this is convention)
  dsetname = filename ( 1:index( filename, '_' )-1 )
  if (mpirank==0) then
    write (*,'("Reading file ",A)') trim(adjustl(filename))
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
    stop
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
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call H5Fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)
  
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
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call H5Dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call H5Dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Pclose_f(plist_id, error) ! note the dataset remains opened

  ! Close dataset
  call H5Dclose_f(dset_id, error)
  call H5Fclose_f(file_id,error)
  call H5close_f(error)
  
  
end subroutine Read_Single_File




! Load backup data from disk to initialize run for restart
subroutine Read_Runtime_Backup(filename,time,dt0,dt1,n1,it,uk,nlk,explin,work)
  use mpi_header
  use vars
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),intent(out) :: time,dt1,dt0
  integer,intent(out) :: n1,it
  complex(kind=pr), intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real(kind=pr),intent(out) :: explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

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
  call H5Pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call H5Pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call H5Fopen_f (filename, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call H5Pclose_f(plist_id, error)

  ! Read fluid backup field:
  call Read_Field_Backup(work,"ux",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,1),work)
  call Read_Field_Backup(work,"uy",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,2),work)
  call Read_Field_Backup(work,"uz",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,3),work)

  if(method == "mhd") then
     ! Read MHD backup field:
     call Read_Field_Backup(work,"bx",time,dt0,dt1,n1,it,file_id)
     call fft(uk(:,:,:,4),work)
     call Read_Field_Backup(work,"by",time,dt0,dt1,n1,it,file_id)
     call fft(uk(:,:,:,5),work)
     call Read_Field_Backup(work,"bz",time,dt0,dt1,n1,it,file_id)
     call fft(uk(:,:,:,6),work)
  endif

  ! Read fluid nonlinear source term backup:
  call Read_Field_Backup(work,"nlkx0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,1,0),work)
  call Read_Field_Backup(work,"nlky0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,2,0),work)
  call Read_Field_Backup(work,"nlkz0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,3,0),work)
  call Read_Field_Backup(work,"nlkx1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,1,1),work)
  call Read_Field_Backup(work,"nlky1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,2,1),work)
  call Read_Field_Backup(work,"nlkz1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,3,1),work)

  if(method == "mhd") then
     ! Read MHD nonlinear source term backup too:
     call Read_Field_Backup(work,"bnlkx0",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,4,0),work)
     call Read_Field_Backup(work,"bnlky0",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,5,0),work)
     call Read_Field_Backup(work,"bnlkz0",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,6,0),work)
     call Read_Field_Backup(work,"bnlkx1",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,4,1),work)
     call Read_Field_Backup(work,"bnlky1",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,5,1),work)
     call Read_Field_Backup(work,"bnlkz1",time,dt0,dt1,n1,it,file_id)
     call fft(nlk(:,:,:,6,1),work)
  endif

  call H5Fclose_f(file_id,error)
  call H5close_f(error)

  ! It is important to have explin, because it won't be initialized
  ! if both time steps dt0 and dt1 match so we compute it here (TOMMY:
  ! are you sure about dt1??? TODO) 
  ! FIXME: only compute if dt0=dt1?
  call cal_vis(dt1,explin)

  if(mpirank == 0) then
     write(*,'("time=",es15.8," dt0=",es15.8)') time, dt0
     write(*,'("!!! DONE READING BACKUP (thats good news!)")')
     write(*,'("---------")')
  endif

end subroutine Read_Runtime_Backup


! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
subroutine Read_Field_Backup(field,dsetname,time,dt0,dt1,n1,it,file_id)
  use mpi_header
  use fsi_vars
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
     call MPI_REDUCE ( dimensions_local(i), chunking_dims(i),1, &
          MPI_INTEGER8, MPI_MAX,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST  ( chunking_dims(i), 1, MPI_INTEGER8, 0, &
          MPI_COMM_WORLD, mpierror )
  enddo

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call H5Screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call H5Screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call H5Pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call H5Dopen_f(file_id, dsetname, dset_id, error)

  ! Select hyperslab in the file.
  call H5Dget_space_f(dset_id, filespace, error)
  call H5Sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
       error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call H5Pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call H5Pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  call H5Dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
       mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  call H5Sclose_f(filespace, error)
  call H5Sclose_f(memspace, error)
  call H5Pclose_f(plist_id, error) ! note the dataset remains opened

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
     stop
  endif

  deallocate (attributes)

  ! Close dataset
  call H5Dclose_f(dset_id, error)
end subroutine Read_Field_Backup


! Write a given attribute with attribute name aname and dimensions
! adims/dims to a given dataset identifier dset_id. Double version.
subroutine write_attribute_dble(adims,aname,attribute,dim,dset_id)
  use mpi_header
  use vars
  use HDF5
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
  use mpi_header
  use vars
  use HDF5
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


! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
subroutine save_fields_new_mhd(time,ubk,ub,wj,nlk,work)
  use mpi_header
  use mhd_vars
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in) :: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  character(len=17) :: name
  integer :: i

  !--Set up file name base
  write(name,'(i5.5)') floor(time*100.d0)
  name=trim(adjustl(name))

  if(mpirank == 0 ) write(*,*) "Saving output fields..."

  ! We need the velocity for saving the velocity and/or vorticity
  if(iSaveVelocity == 1 .or. iSaveVorticity == 1) then
     do i=1,3
        call ifft(ub(:,:,:,i),ubk(:,:,:,i))
     enddo
  endif
    
  ! We need the magnetic fields velocity for saving the magnetic field
  ! and/or current density
  if(iSaveMagneticField == 1  .or. iSaveCurrent == 1) then
     do i=4,6
        call ifft(ub(:,:,:,i),ubk(:,:,:,i))
     enddo
  endif
    
  ! Save the velocity
  if(iSaveVelocity == 1) then
     call Save_Field_HDF5(time,'./ux_'//name,ub(:,:,:,1),"ux")
     call Save_Field_HDF5(time,'./uy_'//name,ub(:,:,:,2),"uy")
     call Save_Field_HDF5(time,'./uz_'//name,ub(:,:,:,3),"uz")
  endif
  
  ! Save the vorticity
  if(iSaveVorticity == 1) then
     ! compute vorticity
     call curl(&
          nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
          ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3)) 
     do i=1,3
        call ifft(wj(:,:,:,i),nlk(:,:,:,i))
     enddo
     call Save_Field_HDF5(time,'./vorx_'//name,wj(:,:,:,1),"vorx")
     call Save_Field_HDF5(time,'./vory_'//name,wj(:,:,:,2),"vory")
     call Save_Field_HDF5(time,'./vorz_'//name,wj(:,:,:,3),"vorz")
  endif
  
  ! Save the magnetic field
  if(iSaveMagneticField == 1) then
     call Save_Field_HDF5(time,'./bx_'//name,ub(:,:,:,4),"bx")
     call Save_Field_HDF5(time,'./by_'//name,ub(:,:,:,5),"by")
     call Save_Field_HDF5(time,'./bz_'//name,ub(:,:,:,6),"bz")
  endif

  ! Save the current density
  if(iSaveCurrent == 1) then
     call curl(&
          nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6),&
          ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6)) 
     do i=4,6
        call ifft(wj(:,:,:,i),nlk(:,:,:,i))
     enddo
     call Save_Field_HDF5(time,'./jx_'//name,wj(:,:,:,4),"jx")
     call Save_Field_HDF5(time,'./jy_'//name,wj(:,:,:,5),"jy")
     call Save_Field_HDF5(time,'./jz_'//name,wj(:,:,:,6),"jz")
  endif
  
  ! Save Mask
  ! FIXME: for stationary masks, this should be done only once
  if((iSaveMask == 1).and.(iPenalization == 1)) then
     call Save_Field_HDF5(time,'./mask_'//name,mask,"mask")
  endif

  if(mpirank == 0 ) write(*,*) "   ...finished saving output fields."
end subroutine save_fields_new_mhd





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
  integer(hid_t) :: atype_id      ! attribute dataspace identifier
  integer(hsize_t), dimension(1) :: adims = (/1/) ! attribute dimension
  integer     ::   arank = 1                      ! attribure rank
  integer(size_t) :: attrlen    ! length of the attribute string

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
subroutine check_file_exists ( fname )
  use vars
  implicit none
  
  character (len=*), intent(in) :: fname
  logical :: exist1

  if (mpirank == 0) then
    inquire ( file=fname, exist=exist1 )
    if ( exist1 .eqv. .false.) then
      write (*,'("ERROR! file: ",A," not found")') trim(fname) 
      stop
    endif  
  endif
  
end subroutine check_file_exists