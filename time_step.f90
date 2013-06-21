subroutine time_step(u,uk,nlk,vort,work,workvis)
  use mpi_header
  use fsi_vars
  implicit none
  
  integer :: inter,it
  integer :: n0=0,n1=1
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  integer :: it_start
  real(kind=pr) :: time,dt0,dt1,t1,t2
  
  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1)
  real (kind=pr),intent (inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent (inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent (inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real (kind=pr),intent (inout) :: workvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  time=0.0

  dt0=1.0d0 ! useful to trigger cal_vis
  dt1=2.d0  ! just add a comment to test branching...

  ! Initialize vorticity or read values from a backup file
  call init_fields(n1,time,it,dt0,dt1,uk,nlk,vort,workvis)
  n0=1 - n1
  it_start=it

  ! Loop over time steps
  t1=MPI_wtime()
  do while((time<=tmax) .and. (it<=nt))
     dt0=dt1

     ! If the mask is time-dependend,we create it here
     if(iMoving == 1) call Create_Mask(time)

     ! Do a fluid time step
     call FluidTimeStep(time,dt0,dt1,n0,n1,u,uk,nlk,vort,work,workvis,it)

     ! Switch time levels
     inter=n1 ; n1=n0 ; n0=inter
     ! Advance in time so that uk contains the evolved field at time
     ! 'time+dt1'
     time=time + dt1
     it=it + 1

     ! Output of integrals after tdrag
     ! FIXME: what does tdrag mean?
     if(mpirank == 0 .and. modulo(it,itdrag) == 0) call write_integrals(time)
     
     ! Output how much time remains
     if(mpirank == 0) call are_we_there_yet(it,it_start,time,t2,t1,dt1)

     ! Output(after tsave)
     if(modulo(time - tstart,tsave) <= dt1) then
        ! Note: we can safely delete nlk(:,:,:,1:nf,n0). for RK2 it
        ! never matters,and for AB2 this is the one to be overwritten
        ! in the next step.  This frees 3 complex arrays.
        ! FIXME: why is the above comment important?
        call save_fields_new(time,uk,u,vort,nlk(:,:,:,:,n0),work)
        
        ! Backup if that's specified in the PARAMS.ini file
        if(iDoBackup == 1) then
          call Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
        endif
     endif
  end do

  if(mpirank==0) write(*,'("Finished time step; did it=",i5," time steps")') it
end subroutine time_step


! Given the velocity in Fourier space and a work array vortk, compute
! the vorticity in phsycial space.  Arrays are 4-dimensional.
subroutine compute_vorticity(vort,vortk,uk)
  use mpi_header
  use fsi_vars
  implicit none

  ! input: velocity field in Fourier space
  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  ! work: vortk, (at output: vorticity in Fourier space)
  complex(kind=pr),intent(inout)::vortk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  ! output: vorticity
  real(kind=pr),intent(out) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  ! imaginary unit
  complex(kind=pr) :: imag
  imag = dcmplx(0.d0,1.d0)
  
  ! comput vorticity in Fourier space:
  do iy=ca(3),cb(3)    ! ky : 0..ny/2-1 ,then,-ny/2..-1
     ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
     do ix=ca(2),cb(2)  ! kx : 0..nx/2
        kx=scalex*dble(ix)
        do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then,-nz/2..-1
           kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           vortk(iz,ix,iy,1)=imag*(ky*uk(iz,ix,iy,3)-kz*uk(iz,ix,iy,2))
           vortk(iz,ix,iy,2)=imag*(kz*uk(iz,ix,iy,1)-kx*uk(iz,ix,iy,3))
           vortk(iz,ix,iy,3)=imag*(kx*uk(iz,ix,iy,2)-ky*uk(iz,ix,iy,1))
        enddo
     enddo
  enddo

  ! Transform to physical space
  call cofitxyz(vortk(:,:,:,1),vort(:,:,:,1))
  call cofitxyz(vortk(:,:,:,2),vort(:,:,:,2))
  call cofitxyz(vortk(:,:,:,3),vort(:,:,:,3))
end subroutine compute_vorticity


! Output how much time remains in the simulation.
subroutine are_we_there_yet(it,it_start,time,t2,t1,dt1)
  use vars
  implicit none

  real(kind=pr) :: time_left,time,t2,t1,dt1
  integer :: it,it_start

  if(modulo(it,300) == 0) then     
     t2= MPI_wtime() - t1
     time_left=(((tmax-time)/dt1)*(t2/dble(it-it_start)))
     write(*,'("time left: ",i3,"d ",i2,"h ",i2,"m ",i2,"s dt=",es7.1)') &
          floor(time_left/(24.d0*3600.d0))   ,&
          floor(mod(time_left,24.*3600.d0)/3600.d0),&
          floor(mod(time_left,3600.d0)/60.d0),&
          floor(mod(mod(time_left,3600.d0),60.d0)),dt1
  endif
end subroutine are_we_there_yet


! Wrapper for writing integral quantities to file
subroutine write_integrals(time)
  use vars
  implicit none

  real(kind=pr) :: time

  select case(method(1:3))
  case("fsi") 
     call write_integrals_fsi(time)
  case("mhd") 
     call write_integrals_mhd(time)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in write_integrals"
     call abort
  end select
end subroutine write_integrals


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time)
  use fsi_vars
  implicit none

  real(kind=pr) :: time
  
  open(14,file='drag_data',status='unknown',position='append')
  write(14,'(7(es12.4,1x))')  time,GlobIntegrals%E_kin,&
       GlobIntegrals%Dissip, GlobIntegrals%Force(1),&
       GlobIntegrals%Force(2),GlobIntegrals%Force(3),&
       GlobIntegrals%Volume
  close(14)
end subroutine write_integrals_fsi


! mhd version of writing integral quantities to disk
subroutine write_integrals_mhd(time)
  use mhd_vars
  implicit none

  real(kind=pr) :: time
 
  !FIXME: do things here?
end subroutine write_integrals_mhd
