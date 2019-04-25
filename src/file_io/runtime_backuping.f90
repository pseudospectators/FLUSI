! Write the restart file. nlk(...,0) and nlk(...,1) are saved, the
! time steps, and what else? FIXME: document what is saved.
subroutine dump_runtime_backup(time,dt0,dt1,n1,it,nbackup,ub,nlk,&
  work,scalars,scalars_rhs,Insect,beams,wings)
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use flexible_model
  use module_insects
  use module_helpers, only : check_file_exists
  implicit none

  real(kind=pr),intent(inout) :: time,dt1,dt0
  integer,intent(inout) :: n1,nbackup,it
  complex(kind=pr),intent(in) :: ub(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(in)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  type(solid), dimension(1), intent(in) :: beams
  type(flexible_wing), dimension (1:nWings), intent (in) :: Wings
  type(diptera), intent(in) :: Insect
  character(len=1) :: scalar_id
  real(kind=pr) :: t1
  integer :: error,j  ! error flags
  character(len=15) :: filename15
#ifndef NOHDF5
  character(len=15) :: filename
#else
  character(len=256) :: filename
  character(len=5) :: suffix
#endif

  t1=MPI_wtime() ! performance diagnostic

#ifndef NOHDF5
  ! Write backup file with parallel HDF5 support
  if(mpirank == 0) then
    write(*,'(80("~"))')
    write(*,'("Dumping runtime_backup",i1,".h5 (time=",es12.4,") to disk....")') nbackup, time
    write(*,*) "Backup type is:", backup_type
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1)') nbackup

  if (root .and. backup_type=="one-file-backup") then
    ! if all fields are written to one file, we need to initialize the empty file first
    call init_empty_file(filename//".h5")
  endif
#else
  ! Write backup files without parallel HDF5 support
  if(mpirank == 0) then
    write(*,'(80("~"))')
    write(*,'("Dumping runtime_backup",i1,".np* (time=",es12.4,") to disk....")') nbackup, time
  endif

  ! Create current filename:
  write(filename,'("runtime_backup",i1)') nbackup
  write(suffix,'(i5.5)') mpirank
  suffix = trim(adjustl(suffix))
  filename = trim(adjustl(filename))//'.np'//suffix

  ! Open file for output
  open(11, file = trim(adjustl(filename)), form='unformatted', access='sequential')

  ! Write attributes
  write(11) time, dt1, dt0, n1, it, nx, ny, nz
#endif

  ! Write the fluid backup field:
  call ifft(ub(:,:,:,1), work)
  call dump_field_backup(filename,work,"ux",time,dt0,dt1,n1,it)
  call ifft(ub(:,:,:,2), work)
  call dump_field_backup(filename,work,"uy",time,dt0,dt1,n1,it)
  call ifft(ub(:,:,:,3), work)
  call dump_field_backup(filename,work,"uz",time,dt0,dt1,n1,it)
  ! Write the fluid nonlinear term backup:
  call ifft(nlk(:,:,:,1,0), work)
  call dump_field_backup(filename,work,"nlkx0",time,dt0,dt1,n1,it)
  call ifft(nlk(:,:,:,2,0), work)
  call dump_field_backup(filename,work,"nlky0",time,dt0,dt1,n1,it)
  call ifft(nlk(:,:,:,3,0), work)
  call dump_field_backup(filename,work,"nlkz0",time,dt0,dt1,n1,it)
  call ifft(nlk(:,:,:,1,1), work)
  call dump_field_backup(filename,work,"nlkx1",time,dt0,dt1,n1,it)
  call ifft(nlk(:,:,:,2,1), work)
  call dump_field_backup(filename,work,"nlky1",time,dt0,dt1,n1,it)
  call ifft(nlk(:,:,:,3,1), work)
  call dump_field_backup(filename,work,"nlkz1",time,dt0,dt1,n1,it)

  if(method == "mhd") then
    ! Write the MHD backup field:
    call ifft(ub(:,:,:,4), work)
    call dump_field_backup(filename,work,"bx",time,dt0,dt1,n1,it)
    call ifft(ub(:,:,:,5), work)
    call dump_field_backup(filename,work,"by",time,dt0,dt1,n1,it)
    call ifft(ub(:,:,:,6), work)
    call dump_field_backup(filename,work,"bz",time,dt0,dt1,n1,it)
    ! Write the MHD backup field:
    call ifft(nlk(:,:,:,4,0), work)
    call dump_field_backup(filename,work,"bnlkx0",time,dt0,dt1,n1,it)
    call ifft(nlk(:,:,:,5,0), work)
    call dump_field_backup(filename,work,"bnlky0",time,dt0,dt1,n1,it)
    call ifft(nlk(:,:,:,6,0), work)
    call dump_field_backup(filename,work,"bnlkz0",time,dt0,dt1,n1,it)
    call ifft(nlk(:,:,:,4,1), work)
    call dump_field_backup(filename,work,"bnlkx1",time,dt0,dt1,n1,it)
    call ifft(nlk(:,:,:,5,1), work)
    call dump_field_backup(filename,work,"bnlky1",time,dt0,dt1,n1,it)
    call ifft(nlk(:,:,:,6,1), work)
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

  ! in the ACM case, pressure is a state variable and needs to be stored / read
  ! when performing runtime backups
  if ((method == "fsi") .and. (equation == "artificial-compressibility")) then
      call ifft(ub(:,:,:,4), work)
      call dump_field_backup(filename,work,"p",time,dt0,dt1,n1,it)
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

#ifdef NOHDF5
  ! close backup file
  close(11)
#endif

  !-------------------------------------------------------------------------
  ! backup for the rigid body solver (free-flight insect)
  !-------------------------------------------------------------------------
  if ((method=="fsi").and.(mpirank==0)) then
    if ((iMask=="Insect").and.(Insect%BodyMotion=="free_flight")) then
      write(filename15,'("runtime_backup",i1)') nbackup
      write (*,'(A)',advance="no") "insect bckp in "//filename15//".rigidsolver"
      open(10, file=filename15//".rigidsolver", form='formatted', status='replace')
      write(10, *) time, Insect%STATE, Insect%RHS_old, Insect%RHS_this
      close(10)
    endif
  endif

  !-------------------------------------------------------------------------
  !-- backup solid solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_solid_model=="yes").and.(method=="fsi")) then
    call dump_solid_backup( time, beams, nbackup )
  endif

  !-------------------------------------------------------------------------
  !-- backup flexible wing solver, if doing active FSI
  !-------------------------------------------------------------------------
  if((use_flexible_wing_model=="yes").and.(method=="fsi")) then
    call dump_flexible_wing_backup( time, wings, nbackup )
  endif

  nbackup = 1 - nbackup
  call toc("IO (dump_runtime_backup)", MPI_wtime() - t1)

  if(mpirank == 0) then
    write(*,'(A)') "done writing backup."
    write(*,'(80("~"))')
  endif
end subroutine dump_runtime_backup



!-------------------------------------------------------------------------------
! Load backup data from disk to initialize run for restart
!-------------------------------------------------------------------------------
subroutine read_runtime_backup(filename2,time,dt0,dt1,n1,it,uk,nlk,explin,work,scalars,scalars_rhs)
  use vars
  use mpi
  use p3dfft_wrapper
  use module_helpers, only : check_file_exists
#ifndef NOHDF5
  use hdf5_wrapper
#endif
  implicit none

  character(len=*),intent(in) :: filename2
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
#ifndef NOHDF5
  character(len=15) :: filename
#else
  character(len=256) :: filename
  character(len=5) :: suffix
#endif

#ifndef NOHDF5
  ! Read backup file with parallel HDF5 support
  filename = filename2(1:15) ! "runtime_backup0" or "runtime_backup1"

  if(mpirank == 0) then
    write(*,'("---------")')
    write(*,'(A)') "!!! I'm trying to resume a backup file: "//trim(adjustl(filename2))
    write(*,'(A)') filename
  endif

  if (backup_type == "one-file-backup") then
    call check_file_exists ( filename2 )
    ! read the attribute
    call read_attribute( filename//".h5", "ux", "bckp", attributes )
  else
    ! read the attribute
    call read_attribute( filename//"_ux.h5", "ux", "bckp", attributes )
  endif

  ! and extract the values
  time    = attributes(1)
  dt1     = attributes(2)
  dt0     = attributes(3)
  n1      = int(attributes(4))
  it      = int(attributes(5))
  nx_file = int(attributes(6))
  ny_file = int(attributes(7))
  nz_file = int(attributes(8))
#else
  ! Read backup files without parallel HDF5 support
  if(mpirank == 0) then
    write(*,'("---------")')
    write(*,'(A)') "!!! I'm trying to resume a backup file: "//filename2(1:15)
  endif
  ! Create current filename:
  write(suffix,'(i5.5)') mpirank
  suffix = trim(adjustl(suffix))
  filename = filename2(1:15)//'.np'//suffix
  ! Open file for output
  call check_file_exists ( trim(adjustl(filename)) )
  open(11, file = trim(adjustl(filename)), form='unformatted', access='sequential', status='old')
  ! read the attributes
  read(11) time, dt1, dt0, n1, it, nx_file, ny_file, nz_file
#endif
  if ((nx/=nx_file).or.(ny/=ny_file).or.(nz/=nz_file)) then
    write (*,'(A)') "ERROR! Resolution mismatch"
    write (*,'("in memory:   nx=",i4," ny=",i4," nz=",i4)') nx,ny,nz
    write (*,'("but in file: nx=",i4," ny=",i4," nz=",i4)') nx_file,ny_file,nz_file
    call abort(77776, "ERROR! read_runtime_backup(): Resolution mismatch")
  endif
  ! Read fluid backup field:
  call read_field_backup(filename,"ux",work)
  call fft(work, uk(:,:,:,1))
  call read_field_backup(filename,"uy",work)
  call fft(work, uk(:,:,:,2))
  call read_field_backup(filename,"uz",work)
  call fft(work, uk(:,:,:,3))
  ! Read fluid nonlinear source term backup:
  call read_field_backup(filename,"nlkx0",work)
  call fft(work, nlk(:,:,:,1,0))
  call read_field_backup(filename,"nlky0",work)
  call fft(work, nlk(:,:,:,2,0))
  call read_field_backup(filename,"nlkz0",work)
  call fft(work, nlk(:,:,:,3,0))
  call read_field_backup(filename,"nlkx1",work)
  call fft(work, nlk(:,:,:,1,1))
  call read_field_backup(filename,"nlky1",work)
  call fft(work, nlk(:,:,:,2,1))
  call read_field_backup(filename,"nlkz1",work)
  call fft(work, nlk(:,:,:,3,1))
  if(method == "mhd") then
    ! Read MHD backup field:
    call read_field_backup(filename,"bx",work)
    call fft(work, uk(:,:,:,4))
    call read_field_backup(filename,"by",work)
    call fft(work, uk(:,:,:,5))
    call read_field_backup(filename,"bz",work)
    call fft(work, uk(:,:,:,6))
    ! Read MHD nonlinear source term backup too:
    call read_field_backup(filename,"bnlkx0",work)
    call fft(work, nlk(:,:,:,4,0))
    call read_field_backup(filename,"bnlky0",work)
    call fft(work, nlk(:,:,:,5,0))
    call read_field_backup(filename,"bnlkz0",work)
    call fft(work, nlk(:,:,:,6,0))
    call read_field_backup(filename,"bnlkx1",work)
    call fft(work, nlk(:,:,:,4,1))
    call read_field_backup(filename,"bnlky1",work)
    call fft(work, nlk(:,:,:,5,1))
    call read_field_backup(filename,"bnlkz1",work)
    call fft(work, nlk(:,:,:,6,1))
  endif

  ! in the ACM case, pressure is a state variable and needs to be stored / read
  ! when performing runtime backups
  if ((method == "fsi") .and. (equation == "artificial-compressibility")) then
      call read_field_backup(filename, "p", work)
      call fft(work, uk(:,:,:,4))
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
    call fft ( work, uk_avg(:,:,:,1) )
    call read_field_backup(filename,"uavgy",work)
    call fft ( work, uk_avg(:,:,:,2) )
    call read_field_backup(filename,"uavgz",work)
    call fft ( work, uk_avg(:,:,:,3) )
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(ekin_avg=="yes")) then
    call read_field_backup(filename,"ekinavg",e_avg)
  endif

  if((method=="fsi").and.(time_avg=="yes").and.(enstrophy_avg=="yes")) then
    call read_field_backup(filename,"Z_avg",Z_avg)
  endif
#ifdef NOHDF5
  ! close backup file
  close(11)
#endif
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
