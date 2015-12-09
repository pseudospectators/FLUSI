! Wrapper for saving fields routine
subroutine save_fields(time,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()


  select case(method)
  case("fsi")
    call save_fields_fsi(time,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  case("mhd")
    call save_fields_mhd(time,uk,u,vort,nlk)
  case default
    if (mpirank == 0) write(*,*) "Error! Unkonwn method in save_fields"
    call abort()
  end select

  time_save=time_save + MPI_wtime() - t1 ! performance analysis
end subroutine save_fields


!-------------------------------------------------------------------------------
! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
! The latest version calls cal_nlk_fsi to avoid redudant code.
!-------------------------------------------------------------------------------
subroutine save_fields_fsi(time,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  use fsi_vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use insect_module
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  real(kind=pr):: volume
  character(len=6) :: name
  character(len=7) :: scalar_name
  integer :: j
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect

  !--Set up file name base
  if ( save_only_one_period == "yes" ) then
    ! overwrite files from last period to save disk space
    ! i.e. t=1.05 is written to t=0.05, as well as 2.05 and 3.05
    write(name,'(i6.6)') nint( (time-dble(floor(time/tsave_period)))*1000.d0 )
  else
    ! name is just the time
    write(name,'(i6.6)') nint(time*1000.d0)
  endif

  if (mpirank == 0 ) then
    write(*,'("Saving data, time= ",es12.4,1x," flags= ",5(i1)," name=",A," ...")',advance='no') &
    time,isaveVelocity,isaveVorticity,isavePress,isaveMask,isaveSolidVelocity,name
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ensure that the mask function is at the right time, if it is not constant
  if (iMoving==1) call create_mask (time, Insect, beams)
  ! if we save the pressure, we must compute the right hand side now:
  if (isavePress==1) then
    call cal_nlk_fsi (time,0,nlk,uk,u,vort,work,workc,Insect)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  ! Velocity (returned in x-space by cal_nlk_fsi)
  !-----------------------------------------------------------------------------
  if (isaveVelocity == 1) then
    if (iSavePress==0) then
      ! the cal_nlk has not been called, and we need to do the IFFT
      call ifft3(ink=uk(:,:,:,1:nd),outx=u(:,:,:,1:nd))
    endif
    call save_field_hdf5(time,"./ux_"//name,u(:,:,:,1))
    call save_field_hdf5(time,"./uy_"//name,u(:,:,:,2))
    call save_field_hdf5(time,"./uz_"//name,u(:,:,:,3))
  endif

  !-----------------------------------------------------------------------------
  ! Pressure
  !-----------------------------------------------------------------------------
  if (isavePress == 1) then
    ! compute pressure (remember NLK is *not* divergence free)
    call pressure( nlk,workc(:,:,:,1) )
    ! total pressure in x-space
    call ifft( ink=workc(:,:,:,1), outx=work(:,:,:,1) )
    ! get actuall pressure (we're in the rotational formulation)
    work(:,:,:,1) = work(:,:,:,1) - 0.5d0*( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
    call save_field_hdf5(time,'./p_'//name,work(:,:,:,1))
  endif

  !-----------------------------------------------------------------------------
  ! Vorticity
  !-----------------------------------------------------------------------------
  if ((isaveVorticity==1).or.(iSaveMagVorticity==1)) then
    !-- compute vorticity:
    call curl( ink=uk, outk=nlk)

    ! nlk=uk
    ! call curl_2nd(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3))
    call ifft3( ink=nlk, outx=vort )
  endif

  if (isaveVorticity==1) then
    !-- save Vorticity
    call save_field_hdf5(time,"./vorx_"//name,vort(:,:,:,1))
    call save_field_hdf5(time,"./vory_"//name,vort(:,:,:,2))
    call save_field_hdf5(time,"./vorz_"//name,vort(:,:,:,3))
  endif

  if (iSaveMagVorticity==1) then
    !-- save vorticity magnitude
    work(:,:,:,1) = dsqrt( vort(:,:,:,1)**2 + vort(:,:,:,1)**2 + vort(:,:,:,1)**2 )
    call save_field_hdf5(time,"./vorabs_"//name,work(:,:,:,1))
  endif


  !-----------------------------------------------------------------------------
  ! Mask
  !-----------------------------------------------------------------------------
  if (isaveMask == 1 .and. iPenalization == 1) then
    mask = mask*eps
    call compute_mask_volume(volume)
    if ((mpirank==0).and.(volume<1.0d-10)) write(*,*) "WARNING: saving empty mask"
    call save_field_hdf5(time,'./mask_'//name,mask)
    ! call save_field_hdf5(time,'./color_'//name,dble(mask_color))
    mask = mask/eps
  endif

  !-----------------------------------------------------------------------------
  ! solid velocity
  !-----------------------------------------------------------------------------
  if (isaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call save_field_hdf5(time,'./usx_'//name,us(:,:,:,1))
    call save_field_hdf5(time,'./usy_'//name,us(:,:,:,2))
    call save_field_hdf5(time,'./usz_'//name,us(:,:,:,3))
  endif

  !-----------------------------------------------------------------------------
  ! passive scalar
  !-----------------------------------------------------------------------------
  if (use_passive_scalar==1) then
    do j=1,n_scalars
      write(scalar_name,'("scalar",i1)') j
      call save_field_hdf5(time,scalar_name//'_'//name,&
           scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j))
    enddo
  endif

  !-----------------------------------------------------------------------------
  ! if time-avg velocity field is computed, save this here
  !-----------------------------------------------------------------------------
  if ((time_avg=="yes").and.(time>=tstart_avg)) then
    if (vel_avg=="yes") then
      ! time averaged velocity is computed in fourier space. therefore, before
      ! saving it here, we have to transform back to x-space
      call ifft3( ink=uk_avg(:,:,:,1:3), outx=vort(:,:,:,1:3) )
      if (save_one_only=="yes") then
        ! we always write into the same file (this would be the usual case)
        call save_field_hdf5(time,"./uavgx_000",vort(:,:,:,1))
        call save_field_hdf5(time,"./uavgy_000",vort(:,:,:,2))
        call save_field_hdf5(time,"./uavgz_000",vort(:,:,:,3))
      else
        ! everytime we save to a new file with the time-code in the name
        call save_field_hdf5(time,"./uavgx_"//name,vort(:,:,:,1))
        call save_field_hdf5(time,"./uavgy_"//name,vort(:,:,:,2))
        call save_field_hdf5(time,"./uavgz_"//name,vort(:,:,:,3))
      endif
    endif
    if (ekin_avg=="yes") then
      if (save_one_only=="yes") then
        call save_field_hdf5(time,"./ekinavg_000",e_avg)
      else
        call save_field_hdf5(time,"./ekinavg_"//name,e_avg)
      endif
    endif
    if (enstrophy_avg=="yes") then
      if (save_one_only=="yes") then
        call save_field_hdf5(time,"./zavg_000",Z_avg)
      else
        call save_field_hdf5(time,"./zavg_"//name,Z_avg)
      endif
    endif
  endif

  if (mpirank==0) write(*,*) " ...DONE!"
end subroutine save_fields_fsi


!-------------------------------------------------------------------------------
! Standart wrapper for the HDF5 library, saves a single array (e.g. one vector component)
! to a single HDF5 file. A bunch of useful attributes (resolution, domain size,
! penalization, viscosity, time) are stored as well.
!-------------------------------------------------------------------------------
subroutine save_field_hdf5(time,filename,field_out)
  use vars
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename

  if (striding<1) call abort("Striding value is bad, exit!")

  if (striding==1) then
    ! save the entire field to disk (no striding)
    call flusi_hdf5_wrapper(time,filename,ra,rb,field_out)
  else
    ! save strided field to disk
    call save_field_hdf5_strided(time,filename,field_out)
  endif
end subroutine save_field_hdf5



!-------------------------------------------------------------------------------
! Same as save_field_hdf5, but the x-axis is not the entire field
! useful for saving single slices (nnx=1) or smaller subsets, or slice collections
! (nnx>nx-1)
!-------------------------------------------------------------------------------
subroutine save_field_hdf5_xvar(time,filename,field_out,nnx)
  use mpi
  use vars
  use hdf5
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(0:nnx-1,ra(2):rb(2),ra(3):rb(3))
  integer, intent(in) :: nnx
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
  integer :: rared(1:3), rbred(1:3)

  ! array bounds on MPI-procs:
  rared = (/0,ra(2),ra(3)/)
  rbred = (/nnx-1,rb(2),rb(3)/)

  call flusi_hdf5_wrapper( time, filename, rared, rbred, field_out)
end subroutine save_field_hdf5_xvar


! Write the restart file. nlk(...,0) and nlk(...,1) are saved, the
! time steps, and what else? FIXME: document what is saved.
subroutine dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,ub,nlk,&
  work,scalars,scalars_rhs,Insect,beams)
  use mpi
  use vars
  use fsi_vars
  use hdf5
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none

  real(kind=pr),intent(inout) :: time,dt1,dt0
  integer,intent(inout) :: n1,nbackup,it
  complex(kind=pr),intent(in) :: ub(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(in)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  type(solid), dimension(1), intent(in) :: beams
  type(diptera), intent(in) :: Insect

  character(len=18) :: filename
  character(len=1) :: scalar_id

  real(kind=pr) :: t1
  integer :: error,j  ! error flags
  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: plist_id      ! property list identifier

  t1=MPI_wtime() ! performance diagnostic

  if(mpirank == 0) then
    write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")',&
    advance='no') nbackup, time
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
  call ifft(work,ub(:,:,:,1))
  call dump_field_backup(work,"ux",time,dt0,dt1,n1,it,file_id)
  call ifft(work,ub(:,:,:,2))
  call dump_field_backup(work,"uy",time,dt0,dt1,n1,it,file_id)
  call ifft(work,ub(:,:,:,3))
  call dump_field_backup(work,"uz",time,dt0,dt1,n1,it,file_id)


  if(method == "mhd") then
    ! Write the MHD backup field:
    call ifft(work,ub(:,:,:,4))
    call dump_field_backup(work,"bx",time,dt0,dt1,n1,it,file_id)
    call ifft(work,ub(:,:,:,5))
    call dump_field_backup(work,"by",time,dt0,dt1,n1,it,file_id)
    call ifft(work,ub(:,:,:,6))
    call dump_field_backup(work,"bz",time,dt0,dt1,n1,it,file_id)
  endif

  ! Write the fluid nonlinear term backup:
  call ifft(work,nlk(:,:,:,1,0))
  call dump_field_backup(work,"nlkx0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,2,0))
  call dump_field_backup(work,"nlky0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,3,0))
  call dump_field_backup(work,"nlkz0",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,1,1))
  call dump_field_backup(work,"nlkx1",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,2,1))
  call dump_field_backup(work,"nlky1",time,dt0,dt1,n1,it,file_id)
  call ifft(work,nlk(:,:,:,3,1))
  call dump_field_backup(work,"nlkz1",time,dt0,dt1,n1,it,file_id)


  if((method=="fsi").and.(use_passive_scalar==1)) then
    do j = 1, n_scalars
      write (scalar_id,'(i1)') j
      call dump_field_backup(scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j),&
           "scalar"//scalar_id,time,dt0,dt1,n1,it,file_id)
      call dump_field_backup(scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,0),&
           "scalar"//scalar_id//"_nlk0",time,dt0,dt1,n1,it,file_id)
      call dump_field_backup(scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,1),&
           "scalar"//scalar_id//"_nlk1",time,dt0,dt1,n1,it,file_id)
    enddo
  endif

  if(method == "mhd") then
    ! Write the MHD backup field:
    call ifft(work,nlk(:,:,:,4,0))
    call dump_field_backup(work,"bnlkx0",time,dt0,dt1,n1,it,file_id)
    call ifft(work,nlk(:,:,:,5,0))
    call dump_field_backup(work,"bnlky0",time,dt0,dt1,n1,it,file_id)
    call ifft(work,nlk(:,:,:,6,0))
    call dump_field_backup(work,"bnlkz0",time,dt0,dt1,n1,it,file_id)
    call ifft(work,nlk(:,:,:,4,1))
    call dump_field_backup(work,"bnlkx1",time,dt0,dt1,n1,it,file_id)
    call ifft(work,nlk(:,:,:,5,1))
    call dump_field_backup(work,"bnlky1",time,dt0,dt1,n1,it,file_id)
    call ifft(work,nlk(:,:,:,6,1))
    call dump_field_backup(work,"bnlkz1",time,dt0,dt1,n1,it,file_id)
  endif

  !-- initialize runnning avg from file
  if((method=="fsi").and.(time_avg=="yes").and.(vel_avg=="yes")) then
    call ifft ( outx=work , ink=uk_avg(:,:,:,1) )
    call dump_field_backup(work,"uavgx",time,dt0,dt1,n1,it,file_id)
    call ifft ( outx=work , ink=uk_avg(:,:,:,2) )
    call dump_field_backup(work,"uavgy",time,dt0,dt1,n1,it,file_id)
    call ifft ( outx=work , ink=uk_avg(:,:,:,3) )
    call dump_field_backup(work,"uavgz",time,dt0,dt1,n1,it,file_id)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(ekin_avg=="yes")) then
    call dump_field_backup(e_avg,"ekinavg",time,dt0,dt1,n1,it,file_id)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(enstrophy_avg=="yes")) then
    call dump_field_backup(Z_avg,"Z_avg",time,dt0,dt1,n1,it,file_id)
  endif
  ! Close the file:
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interfaces and HDF5 library:
  call h5close_f(error)


  !-------------------------------------------------------------------------
  ! backup for the rigid body solver (free-flight insect)
  !-------------------------------------------------------------------------
  if ((method=="fsi").and.(mpirank==0)) then
    if ((iMask=="Insect").and.(Insect%BodyMotion=="free_flight")) then
      write (*,'(A)',advance="no") "insect bckp in "//filename//".rigidsolver"
      open(10, file=filename//".rigidsolver", form='formatted', status='replace')
      write(10, *) time, Insect%STATE, Insect%RHS_old, Insect%RHS_this, Insect%M_body_quaternion
      close(10)
    endif
  endif

  !-------------------------------------------------------------------------
  !-- backup solid solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_solid_model=="yes").and.(method=="fsi")) then
    call dump_solid_backup( time, beams, nbackup )
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
! note you need to know the dimensions and domain decomposition before
! calling it.
subroutine Read_Single_File ( filename, field )
  use mpi
  use vars
  use hdf5
  use helpers
  use basic_operators
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),&
  dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
  intent (out) :: field

  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time, xl_file, yl_file, zl_file, t1
  real (kind=pr)                :: fmax,fmin,favg,viscosity_dummy
  character(len=strlen)             :: dsetname
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
  integer :: error, mpicode  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1
  t1 = MPI_wtime()
  ! the dataset is named the same way as the file: (this is convention)
  dsetname = get_dsetname(filename)
  if (mpirank==0) then
    write(*,'(40("~"))')
    write(*,'("Reading from file ",A)') trim(adjustl(filename))
    write(*,'("dsetname=",A)') trim(adjustl(dsetname))
  endif

  !-----------------------------------------------------------------------------
  ! perform tests
  !-----------------------------------------------------------------------------
  call check_file_exists ( filename )

  ! fetch attributes from file to see if it is a good idea to load it
  call Fetch_attributes( filename,nx_file,ny_file,nz_file,&
  xl_file,yl_file ,zl_file,time, viscosity_dummy )

  if (mpirank==0) then
    write(*,'("nx=",i4," ny=",i4," nz=",i4," time=",g12.4," viscosity=",g12.4)')&
     nx_file,ny_file,nz_file,time,viscosity_dummy
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') xl_file,yl_file ,zl_file
  endif

  ! if the domain size doesn't match, proceed, but yell.
  if ((xl.ne.xl_file).or.(yl.ne.yl_file).or.(zl.ne.zl_file)) then
    if (mpirank == 0) then
      write (*,'(A)') " WARNING! Domain size mismatch."
      write (*,'("in memory:   xl=",es12.4,"yl=",es12.4,"zl=",es12.4)')&
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

  ! check if the size of attributed memory is the size in the file.
  call MPI_ALLREDUCE (dimensions_local(1),nx_file,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (dimensions_local(2),ny_file,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (dimensions_local(3),nz_file,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)

  ! if the resolutions do not match, yell and hang yourself
  if ((nx.ne.nx_file).or.(ny.ne.ny_file).or.(nz.ne.nz_file)) then
    if (mpirank == 0) then
      write (*,'(A)') "ERROR! Resolution mismatch"
      write (*,'(A)') "This happens if ra(:) and rb(:) are not properly initialized."
      write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
      write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nx_file,ny_file,nz_file
      call abort()
    endif
  endif


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

  fmax=fieldmax(field)
  fmin=fieldmin(field)
  favg=fieldmean(field)

  if (mpirank==0) then
    write (*,'("max=",g12.4," min=",g12.4," mean=",g12.4)') fmax,fmin,favg
  endif

  call checknan(field,"recently loaded field")

  if (mpirank==0) then
    write (*,'("Done reading file, Elapsed time=",g12.4,"s")') MPI_wtime() - t1
    write(*,'(40("~"))')
  endif

end subroutine Read_Single_File




! Load backup data from disk to initialize run for restart
subroutine read_runtime_backup(filename,time,dt0,dt1,n1,it,uk,nlk,explin,work,scalars,scalars_rhs)
  use mpi
  use fsi_vars
  use p3dfft_wrapper
  use hdf5
  implicit none

  character(len=*),intent(in) :: filename
  real(kind=pr),intent(out) :: time,dt1,dt0
  integer,intent(out) :: n1,it
  complex(kind=pr), intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(out)::&
  nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
  real(kind=pr),intent(out) :: explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)

  integer :: error  ! Error flag
  integer :: j
  character(len=1) :: scalar_id
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
  call read_field_backup(work,"ux",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,1),work)
  call read_field_backup(work,"uy",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,2),work)
  call read_field_backup(work,"uz",time,dt0,dt1,n1,it,file_id)
  call fft(uk(:,:,:,3),work)


  if(method == "mhd") then
    ! Read MHD backup field:
    call read_field_backup(work,"bx",time,dt0,dt1,n1,it,file_id)
    call fft(uk(:,:,:,4),work)
    call read_field_backup(work,"by",time,dt0,dt1,n1,it,file_id)
    call fft(uk(:,:,:,5),work)
    call read_field_backup(work,"bz",time,dt0,dt1,n1,it,file_id)
    call fft(uk(:,:,:,6),work)
  endif

  ! Read fluid nonlinear source term backup:
  call read_field_backup(work,"nlkx0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,1,0),work)
  call read_field_backup(work,"nlky0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,2,0),work)
  call read_field_backup(work,"nlkz0",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,3,0),work)
  call read_field_backup(work,"nlkx1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,1,1),work)
  call read_field_backup(work,"nlky1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,2,1),work)
  call read_field_backup(work,"nlkz1",time,dt0,dt1,n1,it,file_id)
  call fft(nlk(:,:,:,3,1),work)

  if((method=="fsi").and.(use_passive_scalar==1)) then
    do j = 1, n_scalars
      write (scalar_id,'(i1)') j
      call read_field_backup(scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j),&
           "scalar"//scalar_id,time,dt0,dt1,n1,it,file_id)
      call read_field_backup(scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,0),&
           "scalar"//scalar_id//"_nlk0",time,dt0,dt1,n1,it,file_id)
      call read_field_backup(scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,1),&
           "scalar"//scalar_id//"_nlk1",time,dt0,dt1,n1,it,file_id)
    enddo
  endif

  !-- initialize runnning avg from file
  if((method=="fsi").and.(time_avg=="yes").and.(vel_avg=="yes")) then
    call read_field_backup(work,"uavgx",time,dt0,dt1,n1,it,file_id)
    call fft ( inx=work , outk=uk_avg(:,:,:,1) )
    call read_field_backup(work,"uavgy",time,dt0,dt1,n1,it,file_id)
    call fft ( inx=work , outk=uk_avg(:,:,:,2) )
    call read_field_backup(work,"uavgz",time,dt0,dt1,n1,it,file_id)
    call fft ( inx=work , outk=uk_avg(:,:,:,3) )
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(ekin_avg=="yes")) then
    call read_field_backup(e_avg,"ekinavg",time,dt0,dt1,n1,it,file_id)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(enstrophy_avg=="yes")) then
    call read_field_backup(Z_avg,"Z_avg",time,dt0,dt1,n1,it,file_id)
  endif

  if(method == "mhd") then
    ! Read MHD nonlinear source term backup too:
    call read_field_backup(work,"bnlkx0",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,4,0),work)
    call read_field_backup(work,"bnlky0",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,5,0),work)
    call read_field_backup(work,"bnlkz0",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,6,0),work)
    call read_field_backup(work,"bnlkx1",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,4,1),work)
    call read_field_backup(work,"bnlky1",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,5,1),work)
    call read_field_backup(work,"bnlkz1",time,dt0,dt1,n1,it,file_id)
    call fft(nlk(:,:,:,6,1),work)
  endif

  call h5fclose_f(file_id,error)
  call H5close_f(error)

  ! It is important to have explin, because it won't be initialized
  ! if both time steps dt0 and dt1 match so we compute it here (TOMMY:
  ! are you sure about dt1??? TODO)
  ! FIXME: only compute if dt0=dt1?
  call cal_vis(dt1,explin)

  ! note when we started this run
  tstart = time
  if(mpirank == 0) then
    write(*,'("time=",es15.8," dt0=",es15.8)') time, dt0
    write(*,'("!!! DONE READING BACKUP (thats good news!)")')
    write(*,'("---------")')
  endif

end subroutine read_runtime_backup


! This routine reads a single field "dsetname" from a backup file
! "file_id". the field has the attribute "attributes", which is an 8x1
! array containing scalar backup information
subroutine read_field_backup(field,dsetname,time,dt0,dt1,n1,it,file_id)
  use mpi
  use fsi_vars
  use hdf5
  use basic_operators
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

  if (mpirank==0) then
    write(*,'("Reading ",A," from backup file")') trim(adjustl(dsetname))
  endif
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
  call checknan(field,'recently read backup file!!')
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


! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
subroutine save_fields_mhd(time,ubk,ub,wj,nlk)
  use mpi
  use mhd_vars
  use p3dfft_wrapper
  use basic_operators
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(in) :: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  character(len=17) :: name
  integer :: i

  !--Set up file name base
  write(name,'(i5.5)') floor(time*100.d0)
  name=trim(adjustl(name))

  if(mpirank == 0 ) write(*,*) "Saving output fields..."

  ! We need the velocity for saving the velocity and/or vorticity
  if(isaveVelocity == 1 .or. isaveVorticity == 1) then
    do i=1,3
      call ifft(ub(:,:,:,i),ubk(:,:,:,i))
    enddo
  endif

  ! We need the magnetic fields velocity for saving the magnetic field
  ! and/or current density
  if(isaveMagneticfield == 1  .or. isaveCurrent == 1) then
    do i=4,6
      call ifft(ub(:,:,:,i),ubk(:,:,:,i))
    enddo
  endif

  ! save the velocity
  if(isaveVelocity == 1) then
    call save_field_hdf5(time,'./ux_'//name,ub(:,:,:,1))
    call save_field_hdf5(time,'./uy_'//name,ub(:,:,:,2))
    call save_field_hdf5(time,'./uz_'//name,ub(:,:,:,3))
  endif

  ! save the vorticity
  if(isaveVorticity == 1) then
    ! compute vorticity
    call curl(&
    nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
    ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3))
    do i=1,3
      call ifft(wj(:,:,:,i),nlk(:,:,:,i))
    enddo
    call save_field_hdf5(time,'./vorx_'//name,wj(:,:,:,1))
    call save_field_hdf5(time,'./vory_'//name,wj(:,:,:,2))
    call save_field_hdf5(time,'./vorz_'//name,wj(:,:,:,3))
  endif

  ! save the magnetic field
  if(isaveMagneticfield == 1) then
    call save_field_hdf5(time,'./bx_'//name,ub(:,:,:,4))
    call save_field_hdf5(time,'./by_'//name,ub(:,:,:,5))
    call save_field_hdf5(time,'./bz_'//name,ub(:,:,:,6))
  endif

  ! save the current density
  if(isaveCurrent == 1) then
    call curl(&
    nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6),&
    ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6))
    do i=4,6
      call ifft(wj(:,:,:,i),nlk(:,:,:,i))
    enddo
    call save_field_hdf5(time,'./jx_'//name,wj(:,:,:,4))
    call save_field_hdf5(time,'./jy_'//name,wj(:,:,:,5))
    call save_field_hdf5(time,'./jz_'//name,wj(:,:,:,6))
  endif

  ! save Mask
  ! FIXME: for stationary masks, this should be done only once
  if((isaveMask == 1).and.(iPenalization == 1)) then
    call save_field_hdf5(time,'./mask_'//name,mask)
  endif

  if(mpirank == 0 ) write(*,*) "   ...finished saving output fields."
end subroutine save_fields_mhd



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
subroutine Fetch_attributes( filename, nx, ny, nz, xl, yl ,zl, time, viscosity )
  use hdf5
  use helpers, only : get_dsetname
  use mpi
  implicit none

  integer, parameter :: pr = 8
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time, viscosity

  character(len=*) :: filename  ! file name
  character(len=80) :: dsetname  ! dataset name
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
  integer :: mpirank, mpicode

  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)

  call check_file_exists ( filename )
  ! get dataset name from file name
  dsetname = get_dsetname( filename )

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
  ! open attribute (viscosity)
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  aname = "viscosity"
  CALL h5aopen_f(dset_id, "viscosity", attr_id, error)
  if (error == 0) then
    ! Get dataspace and read
    CALL h5aget_space_f(attr_id, aspace_id, error)
    data_dims(1) = 1
    CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

    viscosity = attr_data
    CALL h5aclose_f(attr_id, error) ! Close the attribute.
    CALL h5sclose_f(aspace_id, error) ! Terminate access to the data space.
  else
    if(mpirank==0) write (*,*) "the file did not contain this attribute!"
    viscosity = 0.d0
  endif

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





!-------------------------------------------------------------------------------
! save a strided field to HDF5. this is certainly not the most elegant way to do
! it, nor the most general, but it works.
!-------------------------------------------------------------------------------
! we figure out what the array bounds of the downsampled array on the local CPU
! are, then copy the data to the smaller field, and then pass both to the HDF5
! wrapper and write it to disk
!-------------------------------------------------------------------------------
subroutine save_field_hdf5_strided(time,filename,field_out)
  use mpi
  use vars
  use hdf5
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
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

  call flusi_hdf5_wrapper( time, filename, rared, rbred, field_red)
end subroutine save_field_hdf5_strided


!-------------------------------------------------------------------------------
! this subroutine writes an array (and just one) to a HDF5 file according to FLUSI
! conventions. Thus, e.g., mask_0000000.h5 is an array at time t=000.000 and con-
! tains a datset "mask" with attributes "nxyz", "domain_size", "epsi", "time" and
! "viscosity"
!-------------------------------------------------------------------------------
! It is a generalization of older routines in the sense that it can save arrays
! of different sizes, which is why the local domain decomposition, i.e. the indices
! of the local CPU's memory, are passed in rared and rbred (opposed to the globals
! ra and rb). This routine does however not work if some CPU have empty arrays,
! so we can only stride and cut the x-direction (since it is contiguous, always!)
!-------------------------------------------------------------------------------
subroutine flusi_hdf5_wrapper( time, filename, rared, rbred, field_out)
  use mpi
  use vars
  use helpers
  use hdf5
  implicit none

  ! The field to be written to disk:
  integer,dimension(1:3), intent(in) :: rared,rbred
  real(kind=pr),intent(in) :: field_out(rared(1):rbred(1),rared(2):rbred(2),rared(3):rbred(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename

  integer, parameter :: rank = 3 ! data dimensionality (2D or 3D)
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
  ! stride is spacing between elements, this is one here. striding is done in the
  ! caller; here, we just write the entire (possibly downsampled) field to disk.
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! HDF attribute variables
  integer, parameter :: arank = 1
  integer(hsize_t), DIMENSION(1) :: adims  ! Attribute dimension

  integer :: mpierror, i, mindim, maxdim, nxred_file, nyred_file, nzred_file, mpicode
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()

  ! ----------------------------------------------------------------------------
  ! Compute the dimension of the complete field (i.e. the union of all CPU's)
  ! which we will write to file.
  ! ----------------------------------------------------------------------------
  call MPI_ALLREDUCE ( rared(1), mindim,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( rbred(1), maxdim,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpicode)
  nxred_file = maxdim-mindim+1

  call MPI_ALLREDUCE ( rared(2), mindim,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( rbred(2), maxdim,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpicode)
  nyred_file = maxdim-mindim+1

  call MPI_ALLREDUCE ( rared(3), mindim,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( rbred(3), maxdim,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpicode)
  nzred_file = maxdim-mindim+1

  ! Tell HDF5 how our  data is organized:
  dimensions_file = (/nxred_file, nyred_file, nzred_file/)
  offset(1) = rared(1)
  offset(2) = rared(2)
  offset(3) = rared(3)
  dimensions_local(1) = rbred(1)-rared(1) + 1
  dimensions_local(2) = rbred(2)-rared(2) + 1
  dimensions_local(3) = rbred(3)-rared(3) + 1

  ! each process knows how much data it has and where to store it.
  ! now, define the dataset chunking. Chunking is largest dimension in
  ! each direction
  do i = 1, 3
    call MPI_REDUCE(dimensions_local(i),chunking_dims(i),1, &
    MPI_INTEGER8,MPI_MIN,0,MPI_COMM_WORLD,mpierror)
    call MPI_BCAST(chunking_dims(i),1,MPI_INTEGER8,0,MPI_COMM_WORLD,mpierror )
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
if ( index(filename,'.h5')==0 ) then
  ! file does not contain *.h5 ending -> add suffix
  call h5fcreate_f('bglockless:'//trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
else
  ! field DOES contain .h5 ending -> just write
  call h5fcreate_f('bglockless:'//trim(adjustl(filename)), H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
endif
#else
if ( index(filename,'.h5')==0 ) then
  ! file does not contain *.h5 ending -> add suffix
  call h5fcreate_f(trim(adjustl(filename))//'.h5', H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
else
  ! field DOES contain .h5 ending -> just write
  call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, &
  file_id, error, access_prp = plist_id)
endif
#endif

  ! this closes the property list plist_id (we'll re-use it)
  call h5pclose_f(plist_id, error)

  !-----------------------------------------------------------------------------
  ! create dataspace "filespace" to write to
  !-----------------------------------------------------------------------------
  ! Create the data space for the  dataset.
  ! Dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)

  ! Create chunked dataset.
  ! NB: chunking and hyperslab are unrelated
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)
  if (field_precision=="double") then
    ! Output files in double precision
    call h5dcreate_f(file_id, get_dsetname(filename), H5T_NATIVE_DOUBLE, filespace, &
    dset_id, error, plist_id)
  else
    ! Output files in single precision
    call h5dcreate_f(file_id, get_dsetname(filename), H5T_NATIVE_REAL, filespace, &
    dset_id, error, plist_id)
  endif
  call h5sclose_f(filespace, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
  error, stride, dimensions_local)

  ! Create property list for collective dataset write
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

  !-----------------------------------------------------------------------------
  ! create dataspace "memspace" to be written
  !-----------------------------------------------------------------------------
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)



  !-----------------------------------------------------------------------------
  ! actual writing of heavy data
  !-----------------------------------------------------------------------------
  ! Write the dataset collectively, double precision in memory
  call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, field_out, dimensions_file, &
  error, file_space_id = filespace, mem_space_id = memspace,&
  xfer_prp = plist_id)

  !!! Write the attributes to the HDF files.
  ! The attributes written are time, penalisation parameter,
  ! computational resolution, and physical domain size.
  adims = (/1/)
  call write_attribute_dble(adims,"time",(/time/),1,dset_id)
  call write_attribute_dble(adims,"viscosity",(/nu/),1,dset_id)
  call write_attribute_dble(adims,"epsi",(/eps/),1,dset_id)
  adims = (/3/)
  call write_attribute_dble(adims,"domain_size",(/xl,yl,zl/),3,dset_id)
  call write_attribute_int(adims,"nxyz",(/nxred_file, nyred_file, nzred_file/),3,dset_id)

  !!! Close dataspaces:
  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5dclose_f(dset_id, error) ! Close the dataset.
  call h5pclose_f(plist_id, error) ! Close the property list.
  call h5fclose_f(file_id, error) ! Close the file.
  call h5close_f(error) ! Close Fortran interfaces and HDF5 library.


  time_hdf5=time_hdf5 + MPI_wtime() - t1 ! performance analysis
end subroutine flusi_hdf5_wrapper
