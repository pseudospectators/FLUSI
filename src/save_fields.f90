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


! Write the restart file. nlk(...,0) and nlk(...,1) are saved, the
! time steps, and what else? FIXME: document what is saved.
subroutine dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,ub,nlk,&
  work,scalars,scalars_rhs,Insect,beams)
  use mpi
  use vars
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

  t1=MPI_wtime() ! performance diagnostic

  if(mpirank == 0) then
    write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")',&
    advance='no') nbackup, time
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1,".h5")') nbackup

  if (root) then
    call init_empty_file(filename)
  endif

  ! Write the fluid backup field:
  call ifft(work,ub(:,:,:,1))
  call dump_field_backup(filename,work,"ux",time,dt0,dt1,n1,it)
  call ifft(work,ub(:,:,:,2))
  call dump_field_backup(filename,work,"uy",time,dt0,dt1,n1,it)
  call ifft(work,ub(:,:,:,3))
  call dump_field_backup(filename,work,"uz",time,dt0,dt1,n1,it)
  ! Write the fluid nonlinear term backup:
  call ifft(work,nlk(:,:,:,1,0))
  call dump_field_backup(filename,work,"nlkx0",time,dt0,dt1,n1,it)
  call ifft(work,nlk(:,:,:,2,0))
  call dump_field_backup(filename,work,"nlky0",time,dt0,dt1,n1,it)
  call ifft(work,nlk(:,:,:,3,0))
  call dump_field_backup(filename,work,"nlkz0",time,dt0,dt1,n1,it)
  call ifft(work,nlk(:,:,:,1,1))
  call dump_field_backup(filename,work,"nlkx1",time,dt0,dt1,n1,it)
  call ifft(work,nlk(:,:,:,2,1))
  call dump_field_backup(filename,work,"nlky1",time,dt0,dt1,n1,it)
  call ifft(work,nlk(:,:,:,3,1))
  call dump_field_backup(filename,work,"nlkz1",time,dt0,dt1,n1,it)

  if(method == "mhd") then
    ! Write the MHD backup field:
    call ifft(work,ub(:,:,:,4))
    call dump_field_backup(filename,work,"bx",time,dt0,dt1,n1,it)
    call ifft(work,ub(:,:,:,5))
    call dump_field_backup(filename,work,"by",time,dt0,dt1,n1,it)
    call ifft(work,ub(:,:,:,6))
    call dump_field_backup(filename,work,"bz",time,dt0,dt1,n1,it)
    ! Write the MHD backup field:
    call ifft(work,nlk(:,:,:,4,0))
    call dump_field_backup(filename,work,"bnlkx0",time,dt0,dt1,n1,it)
    call ifft(work,nlk(:,:,:,5,0))
    call dump_field_backup(filename,work,"bnlky0",time,dt0,dt1,n1,it)
    call ifft(work,nlk(:,:,:,6,0))
    call dump_field_backup(filename,work,"bnlkz0",time,dt0,dt1,n1,it)
    call ifft(work,nlk(:,:,:,4,1))
    call dump_field_backup(filename,work,"bnlkx1",time,dt0,dt1,n1,it)
    call ifft(work,nlk(:,:,:,5,1))
    call dump_field_backup(filename,work,"bnlky1",time,dt0,dt1,n1,it)
    call ifft(work,nlk(:,:,:,6,1))
    call dump_field_backup(filename,work,"bnlkz1",time,dt0,dt1,n1,it)
  endif

  if((method=="fsi").and.(use_passive_scalar==1)) then
    do j = 1, n_scalars
      write (scalar_id,'(i1)') j
      call dump_field_backup(filename,scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j),&
           "scalar"//scalar_id,time,dt0,dt1,n1,it)
      call dump_field_backup(filename,scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,0),&
           "scalar"//scalar_id//"_nlk0",time,dt0,dt1,n1,it)
      call dump_field_backup(filename,scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,1),&
           "scalar"//scalar_id//"_nlk1",time,dt0,dt1,n1,it)
    enddo
  endif

  !-- initialize runnning avg from file
  if((method=="fsi").and.(time_avg=="yes").and.(vel_avg=="yes")) then
    call ifft ( outx=work , ink=uk_avg(:,:,:,1) )
    call dump_field_backup(filename,work,"uavgx",time,dt0,dt1,n1,it)
    call ifft ( outx=work , ink=uk_avg(:,:,:,2) )
    call dump_field_backup(filename,work,"uavgy",time,dt0,dt1,n1,it)
    call ifft ( outx=work , ink=uk_avg(:,:,:,3) )
    call dump_field_backup(filename,work,"uavgz",time,dt0,dt1,n1,it)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(ekin_avg=="yes")) then
    call dump_field_backup(filename,e_avg,"ekinavg",time,dt0,dt1,n1,it)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(enstrophy_avg=="yes")) then
    call dump_field_backup(filename,Z_avg,"Z_avg",time,dt0,dt1,n1,it)
  endif

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



!-------------------------------------------------------------------------------
! This routine dumps a single field "field" as a dataset "dsetname" to
! a backup file "filename". Attributes are stores in one attribute
! "bckp" which contains 8 values
!-------------------------------------------------------------------------------
subroutine dump_field_backup(filename,field,dsetname,time,dt0,dt1,n1,it)
  use vars
  use hdf5_wrapper
  implicit none

  integer,intent(in) :: n1,it
  real(kind=pr), intent (in) :: time,dt1,dt0
  real(kind=pr),intent(in) :: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*), intent (in) :: dsetname, filename

  call write_field_hdf5( filename,dsetname, ra, rb, field, overwrite=.false.)
  call write_attribute( filename,dsetname,"bckp",(/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/) )
end subroutine dump_field_backup



!-------------------------------------------------------------------------------
! Read in a single file that follows the naming convention
! note you need to know the dimensions and domain decomposition before
! calling it.
!-------------------------------------------------------------------------------
subroutine Read_Single_File ( filename, field )
  use vars
  use hdf5_wrapper
  use basic_operators, only : fieldmax, fieldmin, fieldmean, checknan
  use helpers, only : get_dsetname
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

  if (mpirank==0) then
    write(*,'(40("~"))')
    write(*,'("Reading from file ",A)') trim(adjustl(filename))
    write(*,'("dsetname=",A)') trim(adjustl(get_dsetname(filename)))
    write(*,'("nx=",i4," ny=",i4," nz=",i4," time=",g12.4," viscosity=",g16.4)') nxyz,ttime(1),viscosity_dummy(1)
    write(*,'("xl=",g12.4," yl=",g12.4," zl=",g12.4)') domain

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
      call abort(125)
    endif
  endif

  ! actual reading of file
  call read_field_hdf5 ( filename, get_dsetname(filename), ra,rb, field)
  ! check if field contains NaN
  call checknan(field,"recently loaded field")

  fmax = fieldmax(field)
  fmin = fieldmin(field)
  favg = fieldmean(field)
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1.0d+6
  t2 = MPI_wtime() - t1

  if (mpirank==0) then
    write (*,'("read ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') mbyte,t2,mbyte/t2
    write (*,'("max=",g12.4," min=",g12.4," mean=",g12.4)') fmax,fmin,favg
    write (*,'("Done reading file, Elapsed time=",g12.4,"s")') MPI_wtime() - t1
    write(*,'(40("~"))')
  endif


end subroutine Read_Single_File



!-------------------------------------------------------------------------------
! Load backup data from disk to initialize run for restart
!-------------------------------------------------------------------------------
subroutine read_runtime_backup(filename,time,dt0,dt1,n1,it,uk,nlk,explin,work,scalars,scalars_rhs)
  use vars
  use p3dfft_wrapper
  use hdf5_wrapper
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


  if(mpirank == 0) then
    write(*,'("---------")')
    write(*,'(A)') "!!! I'm trying to resume a backup file: "//trim(adjustl(filename))
  endif

  call check_file_exists ( filename )

  ! read the attribute
  call read_attribute( filename, "ux", "bckp", attributes )
  ! and extract the values
  time    = attributes(1)
  dt1     = attributes(2)
  dt0     = attributes(3)
  n1      = int(attributes(4))
  it      = int(attributes(5))
  nx_file = int(attributes(6))
  ny_file = int(attributes(7))
  nz_file = int(attributes(8))

  if ((nx/=nx_file).or.(ny/=ny_file).or.(nz/=nz_file)) then
    write (*,'(A)') "ERROR! Resolution mismatch"
    write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
    write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nx_file,ny_file,nz_file
    call abort(77776)
  endif

  ! Read fluid backup field:
  call read_field_backup(filename,"ux",work)
  call fft(uk(:,:,:,1),work)
  call read_field_backup(filename,"uy",work)
  call fft(uk(:,:,:,2),work)
  call read_field_backup(filename,"uz",work)
  call fft(uk(:,:,:,3),work)
  ! Read fluid nonlinear source term backup:
  call read_field_backup(filename,"nlkx0",work)
  call fft(nlk(:,:,:,1,0),work)
  call read_field_backup(filename,"nlky0",work)
  call fft(nlk(:,:,:,2,0),work)
  call read_field_backup(filename,"nlkz0",work)
  call fft(nlk(:,:,:,3,0),work)
  call read_field_backup(filename,"nlkx1",work)
  call fft(nlk(:,:,:,1,1),work)
  call read_field_backup(filename,"nlky1",work)
  call fft(nlk(:,:,:,2,1),work)
  call read_field_backup(filename,"nlkz1",work)
  call fft(nlk(:,:,:,3,1),work)

  if(method == "mhd") then
    ! Read MHD backup field:
    call read_field_backup(filename,"bx",work)
    call fft(uk(:,:,:,4),work)
    call read_field_backup(filename,"by",work)
    call fft(uk(:,:,:,5),work)
    call read_field_backup(filename,"bz",work)
    call fft(uk(:,:,:,6),work)
    ! Read MHD nonlinear source term backup too:
    call read_field_backup(filename,"bnlkx0",work)
    call fft(nlk(:,:,:,4,0),work)
    call read_field_backup(filename,"bnlky0",work)
    call fft(nlk(:,:,:,5,0),work)
    call read_field_backup(filename,"bnlkz0",work)
    call fft(nlk(:,:,:,6,0),work)
    call read_field_backup(filename,"bnlkx1",work)
    call fft(nlk(:,:,:,4,1),work)
    call read_field_backup(filename,"bnlky1",work)
    call fft(nlk(:,:,:,5,1),work)
    call read_field_backup(filename,"bnlkz1",work)
    call fft(nlk(:,:,:,6,1),work)
  endif


  if((method=="fsi").and.(use_passive_scalar==1)) then
    do j = 1, n_scalars
      write (scalar_id,'(i1)') j
      call read_field_backup(filename,"scalar"//scalar_id, scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j))
      call read_field_backup(filename,"scalar"//scalar_id//"_nlk0",scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,0))
      call read_field_backup(filename,"scalar"//scalar_id//"_nlk1",scalars_rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j,1))
    enddo
  endif

  !-- initialize runnning avg from file
  if((method=="fsi").and.(time_avg=="yes").and.(vel_avg=="yes")) then
    call read_field_backup(filename,"uavgx",work)
    call fft ( inx=work , outk=uk_avg(:,:,:,1) )
    call read_field_backup(filename,"uavgy",work)
    call fft ( inx=work , outk=uk_avg(:,:,:,2) )
    call read_field_backup(filename,"uavgz",work)
    call fft ( inx=work , outk=uk_avg(:,:,:,3) )
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(ekin_avg=="yes")) then
    call read_field_backup(filename,"ekinavg",e_avg)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(enstrophy_avg=="yes")) then
    call read_field_backup(filename,"Z_avg",Z_avg)
  endif

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

  call read_field_hdf5( filename, dsetname, ra, rb, field )
  call checknan(field,'recently read backup file!!')

  mbyte = dble(nx)*dble(ny)*dble(nz)*8.d0/1.0d+6
  t1 = MPI_wtime() -t1
  if (root) write(*,'(".. read ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1
end subroutine read_field_backup


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
  use hdf5_wrapper
  use helpers, only : get_dsetname
  use mpi
  implicit none

  character(len=*), intent(in) :: filename  ! file name
  integer, intent (out) :: nx, ny, nz
  real (kind=pr), intent(out) :: xl,yl,zl, time, viscosity

  real(kind=pr),dimension(1) :: attr_data1, attr_data0
  real(kind=pr),dimension(1:3) :: attr_data2
  integer,dimension(1:3) :: attr_data3

  call check_file_exists ( filename )
  call read_attribute( filename, get_dsetname(filename), "time", attr_data0)
  call read_attribute( filename, get_dsetname(filename), "viscosity", attr_data1)
  call read_attribute( filename, get_dsetname(filename), "domain_size", attr_data2)
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
! Standart wrapper for the HDF5 library, saves a single array (e.g. one vector component)
! to a single HDF5 file. A bunch of useful attributes (resolution, domain size,
! penalization, viscosity, time) are stored as well.
!-------------------------------------------------------------------------------
subroutine save_field_hdf5(time,filename,field_out)
  use vars
  use hdf5_wrapper
  use helpers
  implicit none

  ! The field to be written to disk:
  real(kind=pr),intent(in) :: field_out(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent (in) :: time
  character(len=*), intent (in) :: filename
  character(len=strlen) :: fname,fname2
  integer, dimension(1:3) :: tmp
  real(kind=pr) :: t1, mbyte

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

  if (striding<1) call abort("Striding value is bad, exit!")

  if (striding==1) then
    ! save the entire field to disk (no striding)
    ! write actual field to the file
    call write_field_hdf5( fname, get_dsetname(fname), ra, rb, field_out)
    ! append some useful attributes to the field in the file
    call write_attribute( fname, get_dsetname(fname), "time",(/time/))
    call write_attribute( fname, get_dsetname(fname), "viscosity",(/nu/))
    call write_attribute( fname, get_dsetname(fname), "epsi",(/eps/))
    call write_attribute( fname, get_dsetname(fname), "domain_size",(/xl,yl,zl/))
    call write_attribute( fname, get_dsetname(fname), "nxyz",(/nx,ny,nz/))
  else
    ! save strided field to disk
    call save_field_hdf5_strided(time,fname,field_out)
  endif

  ! footer
  t1 = MPI_wtime() - t1
  mbyte = dble(nx)*dble(ny)*dble(nz)*4.d0/1.0d+6
  if (root) write(*,'(".. wrote ",f7.2," MB in ",f7.2," s (",f7.2,"MB/s)")') &
  mbyte, t1, mbyte/t1
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
  use helpers
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
  call write_attribute( fname, get_dsetname(fname), "nxyz",(/nx/2,ny/2,nz/2/))

  deallocate (field_red)
end subroutine save_field_hdf5_strided
