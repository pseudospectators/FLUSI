
program FLUSI
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  integer                :: mpicode
  real (kind=pr)         :: t1,t2, tmp
  character (len=80)     :: infile

  !---- Initialize MPI, get size and rank, create mpi data types
  call MPI_INIT (mpicode)
  call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
  call MPI_COMM_RANK (MPI_COMM_WORLD,mpirank,mpicode)
  
  time_fft=0.0; time_ifft=0.0; time_vis=0.0; time_mask=0.0;
  time_vor=0.0; time_curl=0.0; time_p=0.0; time_nlk=0.0; time_fluid=0.0;
  time_bckp=0.0; time_save=0.0; time_total=0.0; time_u=0.0;

  if (mpirank == 0) then
    write(*,'(A)') '--------------------------------------'
    write(*,'(A)') '		FLUSI'
    write(*,'(A)') '--------------------------------------'
    write(*,'("Running on ",i3," CPUs")') mpisize
    call system('mkdir fields')
  endif

  mpiinteger 	= MPI_INTEGER
  mpireal 	= MPI_DOUBLE_PRECISION
  mpicomplex 	= MPI_DOUBLE_COMPLEX

  
  
  !---- Read input parameters
  if (mpirank == 0) then
    write(*,'(A)') '*** info: Reading input data...'
  endif

  ! get filename of PARAMS file from command line
  call get_command_argument (1, infile)
  ! read all parameters from that file
  call get_params (infile)

  !---- Step forward in time
  if (mpirank == 0) then
    write(*,'(A)') '--------------------------------------'
    write(*,'(A)') '*** info: Starting time iterations...'
    write(*,'(A)') '--------------------------------------'
  endif
  
  ! overwrite drag_data file
  if ((inicond .ne. 99).and.(mpirank==0)) then
  open  (14, file = 'drag_data', status = 'replace')
  close (14)
  endif
  
  
  t1 = MPI_wtime()
  call time_step ()
  t2 = MPI_wtime() - t1
  
  if (mpirank ==0) then
    write(*,'("$$$ info: total elapsed time time_step=",es12.4, " on ",i2," CPUs")') t2, mpisize
  endif

  !---- Deallocate memory
  if (mpirank == 0) then
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
    write(*,'("workvis: ",es12.4," (",f5.1,"%)")') tmp, 100.0*tmp/time_fluid
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
    write(*,'(A)') '!!! warning: Finalizing computation...'
    write(*,'(A)') '--------------------------------------'
  endif
  
  call fft_free 
  
  !---- Finalize MPI
  call MPI_FINALIZE (mpicode)

end program FLUSI

