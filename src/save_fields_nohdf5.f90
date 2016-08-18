! Wrapper for saving fields routine
subroutine save_fields(time,it,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  real(kind=pr),intent(in) :: time
  integer,intent(in) :: it
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
    call save_fields_fsi(time,it,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  case("mhd")
    call save_fields_mhd(time,uk,u,vort,nlk)
  case default
    if (mpirank == 0) write(*,*) "Error! Unkonwn method in save_fields"
    call abort(9)
  end select

  time_save=time_save + MPI_wtime() - t1 ! performance analysis
end subroutine save_fields


!-------------------------------------------------------------------------------
! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
! The latest version calls cal_nlk_fsi to avoid redudant code.
!-------------------------------------------------------------------------------
subroutine save_fields_fsi(time,it,uk,u,vort,nlk,work,workc,scalars,scalars_rhs,Insect,beams)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use insect_module
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in) :: time
  integer,intent(in) :: it
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  real(kind=pr):: volume, t
  character(len=6) :: name
  character(len=7) :: scalar_name
  integer :: j
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect

  ! set up base name for files. we can use either the physical time or the time step
  ! for this. If we use the time, and we have save_only_one_period == "yes", then
  ! we save t/T for only one period
  select case(naming)
  case('timestep')
    ! write time step-based names
    write(name,'(i6.6)') it
  case default
    ! write time-based names by default
    if ( save_only_one_period == "yes" ) then
      ! overwrite files from last period to save disk space
      ! i.e. t=1.05 is written to t=0.05, as well as 2.05 and 3.05
      write(name,'(i6.6)') nint( (time-dble(floor(time/tsave_period)))*1000.d0 )
    else
      ! name is just the time
      write(name,'(i6.6)') nint(time*1000.d0)
    endif
  end select


  if (mpirank == 0 ) then
    write(*,'(80("~"))')
    write(*,'("Saving data, time= ",g12.4,1x," flags= ",5(i1)," name=",A," ...")') &
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
    call save_field_hdf5(time,"ux_"//name,u(:,:,:,1))
    call save_field_hdf5(time,"uy_"//name,u(:,:,:,2))
    call save_field_hdf5(time,"uz_"//name,u(:,:,:,3))
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
    call save_field_hdf5(time,'p_'//name,work(:,:,:,1))
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
    call save_field_hdf5(time,"vorx_"//name,vort(:,:,:,1))
    call save_field_hdf5(time,"vory_"//name,vort(:,:,:,2))
    call save_field_hdf5(time,"vorz_"//name,vort(:,:,:,3))
  endif

  if (iSaveMagVorticity==1) then
    !-- save vorticity magnitude
    work(:,:,:,1) = dsqrt( vort(:,:,:,1)**2 + vort(:,:,:,2)**2 + vort(:,:,:,3)**2 )
    call save_field_hdf5(time,"vorabs_"//name,work(:,:,:,1))
  endif


  !-----------------------------------------------------------------------------
  ! Mask
  !-----------------------------------------------------------------------------
  if (isaveMask == 1 .and. iPenalization == 1) then
    mask = mask*eps
    call compute_mask_volume(volume)
    if ((mpirank==0).and.(volume<1.0d-10)) write(*,*) "WARNING: saving empty mask"
    call save_field_hdf5(time,'mask_'//name,mask)
    ! call save_field_hdf5(time,'color_'//name,dble(mask_color))
    mask = mask/eps
  endif

  !-----------------------------------------------------------------------------
  ! solid velocity
  !-----------------------------------------------------------------------------
  if (isaveSolidVelocity == 1 .and. iPenalization == 1 .and. iMoving == 1) then
    call save_field_hdf5(time,'usx_'//name,us(:,:,:,1))
    call save_field_hdf5(time,'usy_'//name,us(:,:,:,2))
    call save_field_hdf5(time,'usz_'//name,us(:,:,:,3))
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
        call save_field_hdf5(time,"uavgx_000",vort(:,:,:,1))
        call save_field_hdf5(time,"uavgy_000",vort(:,:,:,2))
        call save_field_hdf5(time,"uavgz_000",vort(:,:,:,3))
      else
        ! everytime we save to a new file with the time-code in the name
        call save_field_hdf5(time,"uavgx_"//name,vort(:,:,:,1))
        call save_field_hdf5(time,"uavgy_"//name,vort(:,:,:,2))
        call save_field_hdf5(time,"uavgz_"//name,vort(:,:,:,3))
      endif
    endif
    if (ekin_avg=="yes") then
      if (save_one_only=="yes") then
        call save_field_hdf5(time,"ekinavg_000",e_avg)
      else
        call save_field_hdf5(time,"ekinavg_"//name,e_avg)
      endif
    endif
    if (enstrophy_avg=="yes") then
      if (save_one_only=="yes") then
        call save_field_hdf5(time,"zavg_000",Z_avg)
      else
        call save_field_hdf5(time,"zavg_"//name,Z_avg)
      endif
    endif
  endif

  if (root) then
    write(*,*) " ...done saving!"
    write(*,'(80("~"))')
  endif
end subroutine save_fields_fsi



! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
subroutine save_fields_mhd(time,ubk,ub,wj,nlk)
  use mpi
  use vars
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
    call save_field_hdf5(time,'ux_'//name,ub(:,:,:,1))
    call save_field_hdf5(time,'uy_'//name,ub(:,:,:,2))
    call save_field_hdf5(time,'uz_'//name,ub(:,:,:,3))
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
    call save_field_hdf5(time,'vorx_'//name,wj(:,:,:,1))
    call save_field_hdf5(time,'vory_'//name,wj(:,:,:,2))
    call save_field_hdf5(time,'vorz_'//name,wj(:,:,:,3))
  endif

  ! save the magnetic field
  if(isaveMagneticfield == 1) then
    call save_field_hdf5(time,'bx_'//name,ub(:,:,:,4))
    call save_field_hdf5(time,'by_'//name,ub(:,:,:,5))
    call save_field_hdf5(time,'bz_'//name,ub(:,:,:,6))
  endif

  ! save the current density
  if(isaveCurrent == 1) then
    call curl(&
    nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6),&
    ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6))
    do i=4,6
      call ifft(wj(:,:,:,i),nlk(:,:,:,i))
    enddo
    call save_field_hdf5(time,'jx_'//name,wj(:,:,:,4))
    call save_field_hdf5(time,'jy_'//name,wj(:,:,:,5))
    call save_field_hdf5(time,'jz_'//name,wj(:,:,:,6))
  endif

  ! save Mask
  ! FIXME: for stationary masks, this should be done only once
  if((isaveMask == 1).and.(iPenalization == 1)) then
    call save_field_hdf5(time,'mask_'//name,mask)
  endif

  if(mpirank == 0 ) write(*,*) "   ...finished saving output fields."
end subroutine save_fields_mhd



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

  real(kind=pr) :: t1 ! diagnostic used for performance analysis.
  t1 = MPI_wtime()
  
  ! Array bounds and sizes
  sz_out(1) = rb(1)-ra(1) +1
  sz_out(2) = rb(2)-ra(2) +1
  sz_out(3) = rb(3)-ra(3) +1

  ! Write Fortran binary file
  open(10, file = trim(adjustl(filename)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3))
  write (10,rec=1) field_out
  close (10)

  time_hdf5=time_hdf5 + MPI_wtime() - t1 ! performance analysis
end subroutine save_field_hdf5






! Write the restart file. nlk(...,0) and nlk(...,1) are saved, the
! time steps, and what else? FIXME: document what is saved.
subroutine dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,ub,nlk,&
  work,scalars,scalars_rhs,Insect,beams)
  use mpi
  use vars
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

  t1=MPI_wtime() ! performance diagnostic

  !!! NOT YET IMPLEMENTED !!!

  time_bckp=time_bckp + MPI_wtime() -t1 ! Performance diagnostic

  if(mpirank == 0) write(*,'(A)') "...DONE!"
end subroutine dump_runtime_backup



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
  open(10, file = trim(adjustl(filename)), form='unformatted', access='direct', recl=sz_out(1)*sz_out(2)*sz_out(3))
  read (10,rec=1) field
  close (10)
  
  if (mpirank==0) then
    write (*,'("...DONE! ")',advance='yes')
  endif
  
end subroutine Read_Single_File




!-------------------------------------------------------------------------------
! Load backup data from disk to initialize run for restart
!-------------------------------------------------------------------------------
subroutine read_runtime_backup(filename,time,dt0,dt1,n1,it,uk,nlk,explin,work,scalars,scalars_rhs)
  use vars
  use p3dfft_wrapper
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
  integer :: j,nx_file,ny_file,nz_file
  character(len=1) :: scalar_id
  real(kind=pr), dimension(1:8) :: attributes


  !!! NOT YET IMPLEMENTED

end subroutine read_runtime_backup


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
subroutine Fetch_attributes( filename, nx, ny, nz, xl, yl ,zl, time, viscosity )
  use helpers, only : get_dsetname
  use vars, only : pr
  use mpi
  implicit none

  character(len=*), intent(in) :: filename  ! file name
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time, viscosity

  real(kind=pr),dimension(1) :: attr_data1, attr_data0
  real(kind=pr),dimension(1:3) :: attr_data2
  integer,dimension(1:3) :: attr_data3

  call check_file_exists ( filename )
  !call read_attribute( filename, get_dsetname(filename), "time", attr_data0)
  !call read_attribute( filename, get_dsetname(filename), "viscosity", attr_data1)
  !call read_attribute( filename, get_dsetname(filename), "domain_size", attr_data2)
  !call read_attribute( filename, get_dsetname(filename), "nxyz", attr_data3)

  time = 0 
  viscosity = 0 
  xl = 0 
  yl = 0 
  zl = 0 
  nx = 0 
  ny = 0 
  nz = 0 
end subroutine Fetch_attributes

