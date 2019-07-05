subroutine save_lagrangian_backup(time,fpart,npart_loc,nbackup)
  use mpi
  use vars
  use hdf5
  use interpolation
  !use lagrangian_vars
  implicit none

  integer, intent(in) :: nbackup, npart_loc
  real(kind=pr),intent(in):: time
  real(kind=pr),intent(in):: fpart(npart_proc,n_lagr)

  integer:: fip, i
  integer(8):: record, sizerec
  character(len=18) :: filename

  real(kind=pr) :: t1
  
  t1=MPI_wtime() 
  !  
  fip=1115
  !sizerec=16*(n_lagr+1) ! ifort
  !sizerec=32*(n_lagr+1) ! gfortran
  sizerec=reclc*(n_lagr+1) ! gfortran

  write(filename,'("PARTICLE_RESTART",i1)') 1-nbackup
  if (mpirank==0) print*,'Write ', filename
  OPEN(unit=fip,file=filename,action="write",status="unknown",form="unformatted", access='direct', recl=sizerec)
    do i=1,npart_loc
       record=int(fpart(i,1),8) ! define the record line (time 0: fpart(1:npart), time 1:fpart(1:npart), ...)    
       write(fip,rec=record) time,fpart(i,:) ! write the date on the record line
   enddo
  CLOSE(fip)

! timing module: sum the time for all time steps
call toc("Particles (backup)", MPI_wtime()-t1)
  
end subroutine save_lagrangian_backup

subroutine read_lagrangian_backup(filename,index_part,xpart)
  use mpi
  use vars
  use hdf5
  use interpolation
  !use lagrangian_vars
  implicit none

  character(len=*),intent(in) :: filename  
  real(kind=pr),intent(out):: index_part(npart)
  real(kind=pr),intent(out):: xpart(npart,1:3)

  integer:: fip, i, npos
  integer(8):: record, sizerec
  real(kind=pr):: time


  if (nx==1) then
     npos=2
  else
     npos=3
  endif

  if(mpirank == 0) then
     write(*,'("---------")')
     write(*,'(A)') "!!! I'm trying to resume a lagrangian backup file: "//filename
     if(inicond(1:8) /= "backup::") then
         write(*,*) 'ATTENTION put inicond=backup -> here inicond(1:8)='//inicond(1:8) 
         stop
     endif
  endif
  
  fip=1115
  !sizerec=16*(n_lagr+1) ! ifort
  !sizerec=32*(n_lagr+1) ! gfortran
  sizerec=reclc*(n_lagr+1) ! gfortran
  
  if (npos==2) then
      xpart(:,1)=xl
      OPEN(unit=fip,file=filename,action="read",status="old",form="unformatted", access='direct', recl=sizerec)
      do i=1,npart  
           read(fip,rec=i) time,index_part(i),xpart(i,2:3) ! write the date on the record line
       enddo
      CLOSE(fip)
      xpart(:,1)=xl
  else
       OPEN(unit=fip,file=filename,action="read",status="old",form="unformatted", access='direct', recl=sizerec)
      do i=1,npart   
           read(fip,rec=i) time,index_part(i),xpart(i,1:3) ! write the date on the record line
       enddo
      CLOSE(fip) 
  endif 
  
  if(mpirank == 0) then
      write(*,'(A)') "!!! Reading lagrangian backup file: "//trim(filename)//" --> Success"
      write(*,'("---------")')
  endif
  
  !if (mpirank==0) then
  !   do i=1,npart 
  !       print*,'time=',time,' xpart=',index_part(i),xpart(i,:)
  !   enddo
  !endif

  
end subroutine read_lagrangian_backup


