program FLUSI
  use mpi_header
  use fsi_vars
  implicit none

  integer                :: mpicode
  real (kind=pr)         :: t1,t2
  character (len=80)     :: infile

  ! Arrays needed for simulation
  real(kind=pr),dimension(:,:,:,:),allocatable :: explin  
  real(kind=pr),dimension(:,:,:,:),allocatable :: u,vort
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  complex(kind=pr),dimension(:,:,:,:,:),allocatable :: nlk  
  real(kind=pr),dimension(:,:,:),allocatable :: work

  ! Set method information in vars module.
  method="fsi" ! We are doing fluid-structure interactions
  nf=1 ! We are evolving one field.
  nd=3*nf ! The one field has three components.
  
  ! Initialize MPI, get size and rank
  call MPI_INIT(mpicode)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpicode)

  time_fft=0.0; time_ifft=0.0; time_vis=0.0; time_mask=0.0;
  time_vor=0.0; time_curl=0.0; time_p=0.0; time_nlk=0.0; time_fluid=0.0;
  time_bckp=0.0; time_save=0.0; time_total=0.0; time_u=0.0;

  if (mpirank == 0) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '  FLUSI'
     write(*,'(A)') '--------------------------------------'
     write(*,'("Running on ",i3," CPUs")') mpisize
  endif

  ! Read input parameters
  allocate(lin(nf)) ! Set up the linear term
  if (mpirank == 0) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(1,infile)
  ! read all parameters from that file
  call get_params(infile)
  ! Initialize FFT
  call fft_initialize 

  !  Set up output directory
  call system('mkdir -p fields') ! NB: this does not work on turing.
  
  ! Overwrite drag_data file? only if we're not resuming a backup!
  if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then 
    open  (14, file = 'drag_data', status = 'replace')
    close (14)
  endif

  ! Print domain decomposition
  if (mpirank == 0) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '*** Domain decomposition:'
     write(*,'(A)') '--------------------------------------'
  endif
  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("mpirank=",i3," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
       &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
       mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)

  ! Allocate memory:
  allocate(explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))
  ! velocity in Fourier space
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd))
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1))
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  ! velocity in physical space
  allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))   
  ! vorticity in physical space
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
     call Create_Mask(0.d0)
  endif

  ! Step forward in time
  call MPI_barrier (MPI_COMM_world, mpicode)
  t1 = MPI_wtime()
  call time_step(u,uk,nlk,vort,work,explin) ! Actual time-stepping function
  t2 = MPI_wtime() - t1
  if (mpirank ==0) then
     write(*,'("$$$ info: total elapsed time time_step=",es12.4, " on ",i2," CPUs")') t2, mpisize
  endif

  if(mpirank==0) then
     write(*,*) "control values for debugging:"
     write(*,'("Ekin=",es15.8," Dissip=",es15.8," F1=",es15.8," F2=",es15.8," F3=",es15.8," Vol=",es15.8)')&
          GlobalIntegrals%Ekin,&
          GlobalIntegrals%Dissip, GlobalIntegrals%Force(1),&
          GlobalIntegrals%Force(2),GlobalIntegrals%Force(3),&
          GlobalIntegrals%Volume
  endif

  ! Deallocate memory
  deallocate(lin)
  deallocate(explin)
  deallocate(vort,work)
  deallocate(u,uk,nlk)
  
  if(iPenalization == 1) then
     deallocate(mask)
     if(iMoving == 1) deallocate(us)
  endif

  ! Show the breakdown of timing information
  if (mpirank == 0) call show_timings(t2)

  call fft_free 
  call MPI_FINALIZE(mpicode)
end program FLUSI


! Output information on where the algorithm spent the most time.
subroutine show_timings(t2)
  use fsi_vars
  implicit none

  real (kind=pr) :: t2,tmp

  write(*,'(A)') '--------------------------------------'
  write(*,'(A)') '*** Timings'
  write(*,'(A)') '--------------------------------------'
  write(*,'("of the total time ",es12.4,", FLUSI spend ",es12.4," (",f5.1,"%) on FFTS")') &
       t2, time_fft+time_ifft,100.0*(time_fft+time_ifft)/t2 
  write(*,'(A)') '--------------------------------------'
  write(*,'("Time Stepping contributions:")')
  write(*,'("Fluid      : ",es12.4," (",f5.1,"%)")') time_fluid, 100.0*time_fluid/t2
  write(*,'("Mask       : ",es12.4," (",f5.1,"%)")') time_mask, 100.0*time_mask/t2
  write(*,'("Save Fields: ",es12.4," (",f5.1,"%)")') time_save, 100.0*time_save/t2
  write(*,'("Backuping  : ",es12.4," (",f5.1,"%)")') time_bckp, 100.0*time_bckp/t2
  tmp = t2 - (time_fluid + time_mask + time_save + time_bckp)
  write(*,'("Misc       : ",es12.4," (",f5.1,"%)")') tmp, 100.0*tmp/t2
  write(*,'(A)') '--------------------------------------'
  write(*,'(A)') "The time spend for the fluid decomposes into:"
  write(*,'("cal_nlk: ",es12.4," (",f5.1,"%)")') time_nlk, 100.0*time_nlk/time_fluid
  write(*,'("cal_vis: ",es12.4," (",f5.1,"%)")') time_vis, 100.0*time_vis/time_fluid
  tmp = time_fluid - time_nlk - time_vis
  write(*,'("explin: ",es12.4," (",f5.1,"%)")') tmp, 100.0*tmp/time_fluid
  write(*,'(A)') '--------------------------------------'
  write(*,'(A)') "cal_nlk decomposes into:"
  write(*,'("ifft(uk)       : ",es12.4," (",f5.1,"%)")') time_u, 100.0*time_u/time_nlk
  write(*,'("curl(uk)       : ",es12.4," (",f5.1,"%)")') time_vor, 100.0*time_vor/time_nlk
  write(*,'("vor x u - chi*u: ",es12.4," (",f5.1,"%)")') time_curl, 100.0*time_curl/time_nlk
  write(*,'("projection     : ",es12.4," (",f5.1,"%)")') time_p, 100.0*time_p/time_nlk
  tmp = time_nlk - time_u - time_vor - time_curl - time_p
  write(*,'("Misc           : ",es12.4," (",f5.1,"%)")') tmp, 100.0*tmp/time_nlk
  write (*,'(A)') "cal_nlk: FFTs and local operations:"
  write(*,'("FFTs           : ",es12.4," (",f5.1,"%)")') time_nlk_fft, 100.0*time_nlk_fft/time_nlk
  tmp = time_nlk-time_nlk_fft
  write(*,'("local          : ",es12.4," (",f5.1,"%)")') tmp, 100.0*tmp/time_nlk

  write(*,'(A)') '--------------------------------------'
  write(*,'(A)') 'Finalizing computation....'
  write(*,'(A)') '--------------------------------------'
end subroutine show_timings
