subroutine time_step 
  use mpi_header
  use fsi_vars
  implicit none

  integer :: inter,it
  integer :: n0=0,n1=1
  integer :: nbackup=0  ! 0 - backup to file runtime_backup0,1 - to
  ! runtime_backup1,2 - no backup
  integer :: it_start
  real(kind=pr) :: time,dt0,dt1,t1,t2,time_left
  real(kind=pr),dimension(:,:,:),allocatable :: workvis  
  real(kind=pr),dimension(:,:,:,:),allocatable :: u,vort
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  complex(kind=pr),dimension(:,:,:,:,:),allocatable :: nlk  
  real(kind=pr),dimension(:,:,:),allocatable :: work

  time=0.0

  dt0=1.0d0 ! useful to trigger cal_vis
  dt1=2.d0  ! just add a comment to test branching...

  ! Allocate memory
  allocate(workvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  ! velocity in fourier space
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3,0:1))
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))   
  ! velocity in phy space
  allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))   
  ! vorticity in phy space
  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  ! Create obstacle mask
  if(iPenalization>0) then 
     ! you need the mask field only if you want to actually do penalization
     allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
     if(iMoving == 1) then
        ! if your obstacle moves,you'll need this field for its velocity field
        allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
     endif
     ! create mask(this one call can be redundant,but who cares.)
     call Create_Mask(time)
  endif

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

     !--Switch time levels
     inter=n1 ; n1=n0 ; n0=inter
     !--Advance in time: at this point,uk contains the velocity field
     !--at time 'time+dt1'
     time=time + dt1
     it=it + 1

     ! Output(after tdrag)
     if(modulo(it,itdrag) == 0) then
        if(mpirank == 0) then
           ! NB: the following line does not compile on turing:
           ! open (14,file='',status='unknown',access='append')
           ! so I replaced it with:
           open(14,file='drag_data',status='unknown',position='append')
           write(14,'(7(es12.4,1x))')  time,GlobIntegrals%E_kin,&
                GlobIntegrals%Dissip, GlobIntegrals%Force(1),&
                GlobIntegrals%Force(2),GlobIntegrals%Force(3),&
                GlobIntegrals%Volume
           close(14)
        endif
     endif

     if(modulo(it,300) == 0) then     
        if(mpirank == 0) then 
           t2= MPI_wtime() - t1
           time_left=(((tmax-time)/dt1)*(t2/dble(it-it_start)))
           write(*,'("time left: ",i3,"d ",i2,"h ",i2,"m ",i2,"s dt=",es7.1)') &
                floor(time_left/(24.d0*3600.d0))   ,&
                floor(mod(time_left,24.*3600.d0)/3600.d0),&
                floor(mod(time_left,3600.d0)/60.d0),&
                floor(mod(mod(time_left,3600.d0),60.d0)),dt1
        endif
     endif

     ! Output(after tsave)
     if(modulo(time - tstart,tsave) <= dt1) then
        ! note: we can safely delete nlk(:,:,:,1:3,n0). for RK2 it
        ! never matters,and for AB2 this is the one to be overwritten
        ! in the next step.  this frees 3 complex arrays
        call save_fields_new(time,uk,u,vort,nlk(:,:,:,:,n0),work)
        ! Backup if that's specified in the PARAMS.ini file
        if(iDoBackup == 1) then
          call Dump_Runtime_Backup(time,dt0,dt1,n1,it,nbackup,uk,nlk,work)
        endif
        if(mpirank ==0 ) then 
           write(*,'("<<< info: done saving, returning to time loop")')
        endif
     endif
  end do

  if(mpirank==0) then
     write(*,*) "control values for debugging:"
     write(*,'("Ekin=",es15.8," Dissip=",es15.8," F1=",es15.8," F2=",es15.8," F3=",es15.8," Vol=",es15.8)')&
          GlobIntegrals%E_kin,&
          GlobIntegrals%Dissip, GlobIntegrals%Force(1),&
          GlobIntegrals%Force(2),GlobIntegrals%Force(3),&
          GlobIntegrals%Volume
     write(*,'("did it=",i5," time steps")') it
  endif

  !-- Deallocate memory
  deallocate(workvis)
  deallocate(vort,work)
  deallocate(u,uk,nlk)

  if(iPenalization == 1) then
     deallocate(mask)
     if(iMoving == 1) deallocate(us)
  endif
  ! End time-sepping loop.
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
