subroutine lagrangian(it,time,dt,velk,vel,fpart,npart_loc,it_lagr) 
 ! Benjamin 03 July 2019
 ! This subroutine is computing the lagrangian tracers
 ! The time-advancement of the particles is based on RK2 scheme
 ! The migration of the particles between procs is also include
  use mpi
  use vars
  use p3dfft_wrapper
  use basic_operators
  use penalization ! mask array etc
  use interpolation
  use ghosts
 implicit none
 
 integer,intent(in):: it
 real(kind=pr),intent(in):: time, dt
 real(kind=pr),intent(in):: vel(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
 complex(kind=pr),intent(in):: velk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) 

 real(kind=pr) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
 complex(kind=pr) :: worck(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) 

 integer,intent(inout):: npart_loc
 integer,intent(inout) :: it_lagr  ! counter for saving Particles file 
 real(kind=pr),intent(inout):: fpart(npart_proc,n_lagr)
 
 real(kind=pr),dimension(:,:,:,:),allocatable :: data_interp  ! Data to interpolate with ghost points
 real(kind=pr),dimension(:,:,:,:),allocatable :: acc, data_lagr 
 real(kind=pr),dimension(1:3)     :: x_interp,  x_interp2
 real(kind=pr),dimension(1:n_interp,0:1) :: value_interp  
 
 integer:: ix, iy, iz, i_interp, i, j, id0, interp_lag, ilagr
 integer:: fip, ndata, npos
 integer(8):: record, record_temp, sizerec
 real(kind=pr):: l(1:3), dl(1:3)
 real(kind=pr):: t0, t1, t2, t3, t4
 character(len=17) :: filename, name

 t0 = MPI_wtime()
 !
 l(1)=xl; l(2)=yl; l(3)=zl;
 dl(1)=dx; dl(2)=dy; dl(3)=dz;
 !

 if (nx==1) then
    npos=2
    x_interp(1)=xl
    x_interp2(1)=xl
 else
   npos=3
 endif

 ! Interpolation: Create the ghost points
 ! Define the data to interpolate
 interp_lag=0
 allocate(data_interp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_interp))
 data_interp=1.d16
     
 ilagr=0
 if (nx==1) then
    id0=1
    ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vel(:,:,:,2)
    ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vel(:,:,:,3)
 else
    id0=0
    ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vel(:,:,:,1)
    ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vel(:,:,:,2)
    ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vel(:,:,:,3)
 endif

 ! Synchronise the ghost points with the different procs   
 call synchronize_ghosts ( data_interp(:,:,:,1:npos), npos)

 !!!!!!!!!  Additionnal Lagrangian statistics !!!!!!!!!
 ! Compute Lagrangian acceleration
 if ((ilagr_acc==1).and.((modulo(it,floor(tsave_part/dt))==0).or.(time==tmax))) then
    interp_lag=1
    allocate(acc(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
    call lagrangian_acc(velk,vel,acc)         
    if (nx==1) then
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=acc(:,:,:,2)
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=acc(:,:,:,3)
    else
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=acc(:,:,:,1)
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=acc(:,:,:,2)
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=acc(:,:,:,3)
   endif
 endif

 ! Compute Lagrangian Vorticity
 if ((ilagr_vort==1).and.((modulo(it,floor(tsave_part/dt))==0).or.(time==tmax))) then
    interp_lag=1
    !-- compute vorticity:
    call curl( ink=velk, outk=worck)
    call ifft3( ink=worck, outx=vort )
      
    if (nx==1) then
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vort(:,:,:,1)
    else
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vort(:,:,:,1)
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vort(:,:,:,2)
       ilagr=ilagr+1; data_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),ilagr)=vort(:,:,:,3)
    endif
 endif       

 if (interp_lag==1) then
  ! Synchronise the ghost points with the different procs
  call synchronize_ghosts ( data_interp(:,:,:,npos+1:n_interp), n_interp-npos)
  if (ilagr_acc==1) deallocate(acc)     
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !--------------------------------------------------------------------------------------------------
 ! Advancement of particles using RK2 time scheme
 t2 = MPI_wtime()
 do i=1,npart_loc  
    do j=1,npos
       x_interp(j+id0)=fpart(i,j+1)
    enddo
    value_interp=0.d0
    ! Interpolation for the current time
    do j=1,npos
       value_interp(j,0) = trilinear_interp((/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/), &
                                            (/dx,dy,dz/), data_interp(:,:,:,j), x_interp, .false. )

       x_interp2(j+id0)=x_interp(j+id0)+dt*value_interp(j,0)
       fpart(i,j+npos+1)=value_interp(j,0)
    enddo
 
    t3 = MPI_wtime()
    if (interp_lag==1) then
       do j=npos+1,n_interp
          value_interp(j,0) = trilinear_interp((/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/), &
                                               (/dx,dy,dz/), data_interp(:,:,:,j), x_interp, .false. )
          fpart(i,j+npos+1)=value_interp(j,0)
       enddo    
    endif
    call toc("Particles (interpolation)", MPI_wtime()-t3)

    ! Interpolation for the time+dt       
    do j=1,npos
       value_interp(j,1) = trilinear_interp((/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/), &
                                            (/dx,dy,dz/), data_interp(:,:,:,j), x_interp2, .false. )
       fpart(i,j+1)=x_interp(j+id0)+dt/2*(value_interp(j,0)+value_interp(j,1)) ! RK2 scheme
    enddo
 enddo
 call toc("Particles (time scheme)", MPI_wtime()-t2)
 !--------------------------------------------------------------------------------------------------

 call migration_part(it,npart_loc, fpart ) ! Migration of the particles
 do i=1,npart_loc
    do j=1,npos
       fpart(i,j+1)=per_x(fpart(i,j+1),l(j+id0),dl(j+id0)) ! periodization of the particle positions
    enddo
 enddo

 ! Check if the particles are well spread into the different procs
 !write(filename,'("PART.",i1)') mpirank
 !if (it==1)  then
 !    open(100+mpirank,file=filename,status='unknown') ! Ascii file
 !    close(100+mpirank,status='delete')
 !endif
 

 ! SAVING DATA IN PARTICLE.evol files
 ! Saving could be improved (compressed data) by using hdf5 library --> future work
 if ((modulo(it,floor(tsave_part/dt))==0).or.(time==tmax)) then
    t1=MPI_wtime()      
    !
    ! Define and write the binary file fp   
    fip=111111
    sizerec=16*(n_lagr+1)
! ACCESS DIRECT
!         OPEN(unit=fip,file="PARTICLE.evol",action="write",status="unknown",form="unformatted", access='direct', recl=sizerec)
!         do i=1,npart_loc
!           record=int(fpart(i,1),8)+it_lagr*npart ! define the record line (time 0: fpart(1:npart), time 1:fpart(1:npart), ...)    
!           write(fip,rec=record) time,fpart(i,:) ! write the data on the record line
!
! ACCESS STREAM
!         OPEN(unit=fip,file="PARTICLE.evol",action="write",status="unknown",form="unformatted", access='stream', &
!         convert='LITTLE_ENDIAN')
!         do i=1,npart_loc
!            record_temp=int(fpart(i,1),8)+it_lagr*npart ! define the record line (time 0: fpart(1:npart), time 1:fpart(1:npart), ...)  
!            record= record_temp*((n_lagr+1)*8 + 1)
!            write(fip,pos=record) time,fpart(i,:) ! write the data on the record line
!ACCESS ONE FILE BY PROC
   write(fip+mpirank) it,time,npart_loc
   write(fip*10+mpirank,*) it,time,npart_loc

   do i=1,npart_loc
      !print*,'floor(fpart(i,1))=',floor(fpart(i,1))
      write(fip+mpirank) floor(fpart(i,1)),fpart(i,2:n_lagr) ! write the data on the record line in file PARTICLE_proc_%i
!
!
!
      if ((fpart(i,1)==1d0).and.(itraj_test==1)) then
         t4=MPI_wtime() 
         if (it_lagr==0)  then
            open(99,file='traj_test',status='unknown') ! Ascii file
            close(99,status='delete')
         endif
         open(99,file='traj_test',status='unknown',position='append') ! Ascii file
             write(99,34) time,fpart(i,:)
         close(99)
         call toc( "Particles (save::traj_test)", MPI_wtime() - t4)
      endif

!      TEST ALL THE PARTICLES IN FORT.1*** files
!      do j=1,npart
!         if (fpart(i,1)==j) then
!            if (it_lagr==0)  then
!               open(1000+j,status='unknown') ! Ascii file
!               close(1000+j,status='delete')
!            endif
!            open(1000+j,status='unknown',position='append') ! Ascii file
!                write(1000+j,34) time,fpart(i,:)
!            close(1000+j)
!         endif
!      enddo

   enddo

 ! CLOSE(fip)
   it_lagr=it_lagr+1

   if (((modulo(time,tsave)<dt).and.(time>=tsave_first)).or.(time==tmax)) then
      write(filename,'(es12.6)') time
      do i=1,npart_loc
         open(2,file="cal_traj_time_"//trim(filename),form="formatted",status='unknown',position='append')
              write(2,34) time,fpart(i,:)         
         close(2)   
      enddo
  endif
  call toc("Particles (save)", MPI_wtime()-t1)

!     endif
 
  !! Check if the particles are well spread into the different procs
!  write(filename,'("PART.",i1)') mpirank
!  if (it==1)  open(100+mpirank,file=filename,status='replace') ! Ascii file
!  write(filename,'("PART.",i1)') mpirank         
!  open(100+mpirank,file=filename,position="append") ! Ascii file
!      do i=1,npart_loc
!         write(100+mpirank,31) time,fpart(i,:)
!      enddo               
!  close(100+mpirank)
 endif
 deallocate(data_interp)
 
 ! save timing
 call toc( "Particles (lagrangian MAIN)", MPI_wtime() - t0)

24 FORMAT(4(1X,E20.10)) 
25 FORMAT(5(1X,E20.10))
26 FORMAT(6(1X,E20.10))
27 FORMAT(7(1X,E20.10))   
28 FORMAT(8(1X,E20.10)) 
29 FORMAT(9(1X,E20.10))
30 FORMAT(10(1X,E20.10))
31 FORMAT(11(1X,E20.10))     
33 FORMAT(13(1X,E20.10))    
34 FORMAT(14(1X,E20.10)) 
end subroutine lagrangian
