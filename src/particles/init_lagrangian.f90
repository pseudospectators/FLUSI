subroutine init_lagrangian(fpart,npart_loc)
 ! Benjamin 30 September 2015
 ! This subroutine is computing the lagrangian tracers
 ! The time-advancement of the particles is based on RK2 scheme
 ! The migration of the particles between procs is also include
 use mpi
 use vars
 use p3dfft_wrapper
 use basic_operators
 use penalization ! mask array etc
 use interpolation
 use qsort_mod
 !use particles 
 use penalization ! mask array etc
 implicit none
 
 integer,intent(out):: npart_loc
 real(kind=pr),intent(out):: fpart(npart_proc,n_lagr)

 integer:: ix, iy, iz, ipart, i, index, index2, npart_all_proc, nskip, npart_max, itest, npos
 integer :: mpicode, status(MPI_status_size)
 integer(8):: record, sizerec
 integer,dimension(8) :: values
 real (kind=pr) :: rand, xd(3), tmp, r_blob
 character(len=17) :: filename
 
 real(kind=pr), dimension(:,:), allocatable :: xpart
 real(kind=pr), dimension(:), allocatable :: index_part, xm, ym, zm

 allocate(xpart(npart,1:3),index_part(npart))
 allocate(xm(0:nx-1),ym(0:ny-1),zm(0:nz-1))

 ! Define the mesh
 xm=0.d0; ym=0.d0; zm=0.d0
 do ix = 0, nx-1
    xm(ix) = dble(ix)*dx
 enddo
 do iy = 0, ny-1
    ym(iy) = dble(iy)*dy
 enddo
 do iz = 0, nz-1
    zm(iz) = dble(iz)*dz
 enddo


 if (mpirank==0) print*,'Initialize the particles'
 if (nx==1) then
    npos=2
    if (mpirank==0) print*,"2D case"
 else
   npos=3
   if (mpirank==0) print*,"3D case"
 endif
 !
 index=0
 index2=1
 xpart=0.d0
 xd(1)=(rb(1)-ra(1)+1)*dx
 xd(2)=(rb(2)-ra(2)+1)*dy
 xd(3)=(rb(3)-ra(3)+1)*dz
 
 if (iPenalization==1) then
    if (trim(inicond_lagrangian)=='random')       inicond_lagrangian='random_mask'
 endif
 
 select case(inicond_lagrangian)
    
 case("random")
     index2=-1
     if (mpisize>1) stop "Only for one proc"
     !--------- RANDOM SEED ----------
     !call random_seed()       
     !call date_and_time(VALUES=values)
     !do i=1,values(8)
     !   call random_number(rand)
     !enddo
     !---------------------------------   
     do ipart=1,npart     
        call random_number(rand)
        xpart(ipart,1)=rand*xl
        if (nx==1) xpart(ipart,1)=xl
        call random_number(rand)
        xpart(ipart,2)=rand*yl
        call random_number(rand)
        xpart(ipart,3)=rand*zl  
     enddo
    
 case("random_mask")
     index2=-1
     if (mpisize>1) stop "Only for one proc"
     !--------- RANDOM SEED ----------
     !call random_seed()       
     !call date_and_time(VALUES=values)
     !do i=1,values(8)
     !   call random_number(rand)
     !enddo
     !---------------------------------   
     do ipart=1,npart
        itest=0
        if (nx==1) then
            do while (itest==0)            
                       call random_number(rand)
                       xpart(ipart,1)=rand*xl
                       if (nx==1) xpart(ipart,1)=xl
                       call random_number(rand)
                       xpart(ipart,2)=rand*yl
                       call random_number(rand)
                       xpart(ipart,3)=rand*zl
            
                       if (xpart(ipart,2)>(yl-dy)) xpart(ipart,2)=(yl-dy)
                       if (xpart(ipart,3)>(zl-dz)) xpart(ipart,3)=(zl-dz)                   
                       ix=0
                       iy=int(xpart(ipart,2)/dy)
                       iz=int(xpart(ipart,3)/dz)
                   
                       if (iy>=ny) iy=ny-1
                       if (iz>=nz) iz=nz-1
                       
                       if ((mask(ix,iy,iz)>0)  .or. (mask(ix,iy+1,iz)>0)  .or. &
                          (mask(ix,iy,iz+1)>0).or. (mask(ix,iy+1,iz+1)>0)  )  then
                           itest=0
                       else
                           itest=1
                       endif
            enddo
        else
             do while (itest==0)            
                       call random_number(rand)
                       xpart(ipart,1)=rand*xl
                       call random_number(rand)
                       xpart(ipart,2)=rand*yl
                       call random_number(rand)
                       xpart(ipart,3)=rand*zl
            
                       if (xpart(ipart,1)>(xl-dx)) xpart(ipart,1)=(xl-dx)
                       if (xpart(ipart,2)>(yl-dy)) xpart(ipart,2)=(yl-dy)
                       if (xpart(ipart,3)>(zl-dz)) xpart(ipart,3)=(zl-dz)                   
                       ix=int(xpart(ipart,1)/dx)
                       iy=int(xpart(ipart,2)/dy)
                       iz=int(xpart(ipart,3)/dz)
                   
                       if (ix>=nx) ix=nx-1
                       if (iy>=ny) iy=ny-1
                       if (iz>=nz) iz=nz-1
                       
                       if ( (mask(ix,iy,iz)>0)     .or. (mask(ix+1,iy,iz)>0)   .or. &
                            (mask(ix,iy+1,iz)>0)   .or. (mask(ix+1,iy+1,iz)>0) .or. &
                            (mask(ix,iy,iz+1)>0)   .or. (mask(ix+1,iy,iz+1)>0) .or. &
                            (mask(ix,iy+1,iz+1)>0) .or. (mask(ix+1,iy+1,iz)>0) ) then
                           itest=0
                       else
                           itest=1
                       endif
            enddo           
        endif
     enddo
        
  case("blob") 
     index2=-1
     if (mpisize>1) stop "Only for one proc"
     r_blob=1.87/2.
     do ipart=1,npart
        itest=0
        if (nx==1) then
            do while (itest==0)            
                       call random_number(rand)
                       xpart(ipart,1)=rand*xl
                       if (nx==1) xpart(ipart,1)=xl
                       call random_number(rand)
                       xpart(ipart,2)=rand*yl
                       call random_number(rand)
                       xpart(ipart,3)=rand*zl
            
                       if (xpart(ipart,2)>(yl-dy)) xpart(ipart,2)=(yl-dy)
                       if (xpart(ipart,3)>(zl-dz)) xpart(ipart,3)=(zl-dz)                   
                       ix=1
                       iy=int(xpart(ipart,2)/dy)
                       iz=int(xpart(ipart,3)/dz)
                   
                       if (iy>=ny) iy=ny-1
                       if (iz>=nz) iz=nz-1
                       
                       if ( sqrt( (xpart(ipart,2)-yl/2)**2+ (xpart(ipart,3)-zl/2-2.)**2 ) >r_blob )  then
                           itest=0
                       else
                           itest=1
                       endif
             enddo
         endif
     enddo

 case("test")
     index2=-1
     if (mpisize>1) stop "Only for one proc"
     !npart=5
     !npart_proc=10
     !deallocate(xpart,index_part)
     !allocate(xpart(npart,1:3),index_part(npart))
     ! circle
     xpart(1,1)=xl;   xpart(1,2)=3*yl/4;     xpart(1,3)=zl-2*dz;
     xpart(2,1)=xl;   xpart(2,2)=3*yl/4;     xpart(2,3)=3*zl/4;
     xpart(3,1)=xl;   xpart(3,2)=yl/2+10*dy; xpart(3,3)=zl/2+10*dy;
     xpart(4,1)=xl;   xpart(4,2)=yl/2+dy;    xpart(4,3)=zl/2+dz;
     xpart(5,1)=xl;   xpart(5,2)=yl/2+5*dy;  xpart(5,3)=zl/2+5*dz;
       
 case("infile")
     open(101,file='PARTICLE_READ',status='OLD') ! Ascii file
         do i=1,npart
            read(101,24) tmp,xpart(i,:)
         enddo
     close(101) 
     
 case default
    if(inicond_lagrangian(1:8) == "backup::") then
       call read_lagrangian_backup(inicond_lagrangian(9:len(inicond_lagrangian)),index_part,xpart)
       index=1     
    else
      if (mpirank==0) then
         !write(*,*) "init_lagrangian:: unknown inicond_lagrangian="//inicond_lagrangian
         call abort(123789,"init_lagrangian:: unknown inicond_lagrangian="//inicond_lagrangian)
      endif
    endif
 end select  

 if (index==1) then
     ! Define the particles into the corresponding proc
     npart_loc=0
     ipart=0
     do i=1,npart
        !iy = floor(xpart(i,2)/dy)
        !iz = floor(xpart(i,3)/dz)
        !!if (((ix>=ra(1)).and.(ix<=rb(1)).or.(nx==1)).and.&
        !if ((iy>=ra(2)).and.(iy<=rb(2)).and.&
        !    (iz>=ra(3)).and.(iz<=rb(3))) then

         if ( (xpart(i,2)>=ym(ra(2))) .and. (xpart(i,2)<=ym(rb(2)))  .and. &
              (xpart(i,3)>=zm(ra(3))) .and. (xpart(i,3)<=zm(rb(3))) ) then
            ipart=ipart+1
            npart_loc=npart_loc+1
            fpart(ipart,1)=index_part(i)
            if (nx==1) then
                fpart(ipart,2:3)=xpart(i,2:3)
            else
                fpart(ipart,2:4)=xpart(i,1:3)
            endif
            !print*,'npart_loc=',npart_loc,', mpirank=',mpirank,',xpart=',xpart(i,:)
            !print*,'xpart=',xpart(i,:)
            !print*,'fpart=',fpart(ipart,:)
        endif   
     enddo
     !print*,''
     !do i=1,npart_loc
     !     print*,'fpart=',fpart(i,:)
     !enddo
 elseif (index==0) then
     fpart=1.d16
     ! Define the particles into the corresponding proc
     npart_loc=0
     ipart=0
     do i=1,npart
        ix = floor(xpart(i,1)/dx)
        if (nx==1) ix=0
        iy = floor(xpart(i,2)/dy)
        iz = floor(xpart(i,3)/dz) 
        !if (( (ix>=ra(1)).and.(ix<=rb(1)).or.(nx==1)).and.&
         if( (iy>=ra(2)).and.(iy<=rb(2)).and.&
             (iz>=ra(3)).and.(iz<=rb(3))) then

         !if ( (xpart(i,2)>=ym(ra(2))) .and. (xpart(i,2)<=ym(rb(2)))  .and. &
         !     (xpart(i,3)>=zm(ra(3))) .and. (xpart(i,3)<=zm(rb(3))) ) then
            
            ipart=ipart+1
            npart_loc=npart_loc+1
            fpart(ipart,1)=i
            if (nx==1) then
                fpart(ipart,2:3)=xpart(i,2:3)
            else
                fpart(ipart,2:4)=xpart(i,1:3)
            endif
            !print*,'npart_loc=',npart_loc,', mpirank=',mpirank,',xpart=',xpart(i,:)
        endif   
     enddo
     
     call MPI_ALLREDUCE (npart_loc, npart_max, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpicode)
     if (npart_max>npart_proc) then
         print*,"npart_max>npart_proc:",npart_max,">",npart_proc
         stop
     endif
 else
     !
 endif
 
     ! print*,''
     !do i=1,npart_loc
     !     print*,'fpart=',fpart(i,:)
     !enddo
     !pause
 !print*,''
 !print*,'npart_loc=',npart_loc,', mpirank=',mpirank
 
 call MPI_ALLREDUCE (npart_loc, npart_all_proc, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpicode)
 !print*,'npart_all_proc=',npart_all_proc
 if (npart_all_proc.ne.npart) then
   call MPI_FINALIZE(mpicode)
   print*,'mpirank',mpirank,npart_all_proc
   stop "npart_all_proc \= npart"
   call exit(0)
 endif
     
 
 if (index<0) then 
    sizerec=16*5
    if (nd==2) sizerec=16*4
    open(99,file='ORIG',status='unknown',form="unformatted", access='direct', recl=sizerec) ! Ascii file
        do i=1,npart_loc
           record=int(fpart(i,1))
           write(99,rec=record) fpart(i,1),fpart(i,2:npos+1)
        enddo
    close(99)
 endif
 
 call MPI_barrier (MPI_COMM_world, mpicode)  
    
 if (index<0) then    
    if(mpirank==0) then
        sizerec=16*5
        if (nd==2) sizerec=16*4        
        open(100,file='ORIG',status='old',form="unformatted", access='direct', recl=sizerec) ! Ascii file
           do i=1,npart
              if (nd==2) then
                  read(100,rec=i) index_part(i),xpart(i,2:3)
              else
                  read(100,rec=i) index_part(i),xpart(i,:)
              endif
           enddo
        close(100)
        
        if (nx==1) xpart(:,1)=xl
        open(101,file='PARTICLE_ORIG',status='replace') ! Ascii file
           do i=1,npart
              !write(101,24) dble(i),xpart(i,:)
              write(101,24) fpart(i,1),fpart(i,2:npos+1)
           enddo
       close(101)         
    endif
    
 else
    if(mpirank==0) then
       open(99,file='PARTICLE_ORIG',status='replace') ! Ascii file
           do i=1,npart
              write(99,24) dble(i),xpart(i,:)
           enddo
       close(99)   
    endif
 endif

if (index2<0) then
   call MPI_FINALIZE(mpicode)
   stop "Reload Simulations with option: infile"
   call exit(0)
endif
  
 !! Check if the particles are well spread into the different procs
 !write(filename,'("ORIG.",i1)') mpirank
 !open(100+mpirank,file=filename,status='replace') ! Ascii file
 !  do i=1,npart_loc
 !     !print*,'fpart=', fpart(i,:)
 !     if (nd==2) write(100+mpirank,23) fpart(i,1:3)
 !     if (nd==3) write(100+mpirank,24) fpart(i,1:4)      
 !  enddo  
 !close(100+mpirank)
 !
23 FORMAT(3(1X,E20.10)) 
24 FORMAT(4(1X,E20.10))  
end subroutine init_lagrangian
