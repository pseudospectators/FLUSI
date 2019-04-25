! Wrapper for saving fields routine
subroutine save_fields(time,it,uk,u,vort,nlk,work,workc,press,scalars,scalars_rhs,Insect,beams,wings)
  use vars
  use solid_model
  use flexible_model
  use module_insects
  implicit none

  real(kind=pr),intent(in) :: time
  integer,intent(in) :: it
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  type(flexible_wing),dimension(1:nWings), intent(inout) :: Wings
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: t1 ! diagnostic used for performance analysis.

  t1 = MPI_wtime()

  select case(method)
  case("fsi")
    call save_fields_fsi(time,it,uk,u,vort,nlk,work,workc,press,scalars,scalars_rhs,Insect,beams,wings)
  case("mhd")
    call save_fields_mhd(time,uk,u,vort,nlk)
  case default
    call abort(9,"Error! Unknown method in save_fields")
  end select

  call toc("IO (save_fields to HDF5)", MPI_wtime() - t1)
end subroutine save_fields


!-------------------------------------------------------------------------------
! Main save routine for fields for fsi. it computes missing values
! (such as p and vorticity) and stores the fields in several HDF5
! files.
! The latest version calls cal_nlk_fsi to avoid redudant code.
!-------------------------------------------------------------------------------
subroutine save_fields_fsi(time,it,uk,u,vort,nlk,work,workc,press,scalars,scalars_rhs,Insect,beams,wings)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use flexible_model
  use module_insects
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in) :: time
  integer,intent(in) :: it
  complex(kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  real(kind=pr):: volume, t
  character(len=6) :: name
  character(len=7) :: scalar_name
  integer :: j
  type(flexible_wing),dimension(1:nWings), intent(inout) :: Wings
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
  if (iMoving==1) call create_mask (time, Insect, beams, wings)
  ! if we save the pressure, we must compute the right hand side now:
  if (isavePress==1 .and. equation/='artificial-compressibility') then
    call cal_nlk_fsi (time,0,nlk,uk,u,vort,work,workc,Insect)
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !-----------------------------------------------------------------------------
  ! Velocity (returned in x-space by cal_nlk_fsi)
  !-----------------------------------------------------------------------------
  if (isaveVelocity == 1) then
    if (iSavePress==0) then
      ! the cal_nlk has not been called, and we need to do the IFFT
      call ifft3(ink=uk(:,:,:,1:nd), outx=u(:,:,:,1:nd))
    endif
    call save_field_hdf5(time, "ux_"//name,u(:,:,:,1))
    call save_field_hdf5(time, "uy_"//name,u(:,:,:,2))
    call save_field_hdf5(time, "uz_"//name,u(:,:,:,3))
  endif

  !-----------------------------------------------------------------------------
  ! Pressure
  !-----------------------------------------------------------------------------
  if (isavePress == 1 .and. equation/="artificial-compressibility") then
      ! compute pressure (remember NLK is *not* divergence free)
      call pressure( nlk,workc(:,:,:,1) )
      ! total pressure in x-space
      call ifft( ink=workc(:,:,:,1), outx=work(:,:,:,1) )
      ! get actuall pressure (we're in the rotational formulation)
      work(:,:,:,1) = work(:,:,:,1) - 0.5d0*( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
      call save_field_hdf5(time,'p_'//name,work(:,:,:,1))

  elseif (isavePress == 1 .and. equation=="artificial-compressibility") then
    call ifft( ink=uk(:,:,:,4), outx=work(:,:,:,1) )
    call save_field_hdf5(time, "p_"//name, work(:,:,:,1))

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
    call compute_mask_volume(volume)
    if ((mpirank==0).and.(volume<1.0d-10)) write(*,*) "WARNING: saving empty mask"
    call save_field_hdf5(time, 'mask_'//name, mask)
    ! call save_field_hdf5(time,'color_'//name,dble(mask_color))
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
      call ifft(ubk(:,:,:,i), ub(:,:,:,i))
    enddo
  endif

  ! We need the magnetic fields velocity for saving the magnetic field
  ! and/or current density
  if(isaveMagneticfield == 1  .or. isaveCurrent == 1) then
    do i=4,6
      call ifft(ubk(:,:,:,i), ub(:,:,:,i))
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
      call ifft(nlk(:,:,:,i), wj(:,:,:,i))
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
      call ifft(nlk(:,:,:,i), wj(:,:,:,i))
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
