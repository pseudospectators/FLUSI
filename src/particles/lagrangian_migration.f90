subroutine migration_part(it,npart_loc,fpart)
 ! Benjamin 03 July 2019   
 ! This subroutine transfer the particles to the corresponding particles
 ! input: the index of time = it
 ! input: the number of particles allowed by proc = npart_proc
 ! input and output: the number of particles present in the current proc
 ! input and output: the files containing the index, the position and the data interpolated (vel, acc, ...)
 !                   of the particles present in the current proc
 use vars
 use mpi
 use qsort_mod
 !use particles 
 implicit none
 integer, intent(in)    :: it
 integer, intent(inout) :: npart_loc
 real(kind=pr), intent(inout) :: fpart(npart_proc,n_lagr)
 
 integer :: mpicode,idestination, iorigin, statu(MPI_status_size)
 integer :: disp,dir,source,dest
 integer :: i, imigr, iy, iz, nd_mig, ipart
 integer :: nindex, nskip, id_migr, id(2), coord(2), cordout(2)
 real(kind=pr):: l(1:3), dl(1:3), t_migr

 integer, dimension(:), allocatable :: order, n_migr_loc 
 integer, dimension(:,:), allocatable :: n_migr, ind_fpart
 real(kind=pr), dimension(:,:), allocatable :: temp_sort
 real(kind=pr), dimension(:,:,:), allocatable :: value_temp
 type (group), dimension(:), allocatable :: A
 
 
 t_migr=MPI_wtime() 

 l(1)=xl; l(2)=yl; l(3)=zl;
 dl(1)=dx; dl(2)=dy; dl(3)=dz;
  
 if (mpidims(2)==1) then
   nd_mig=2          
 elseif (mpidims(2)>1) then
   nd_mig=8
 endif

 allocate(ind_fpart(nd_mig,npart_proc),n_migr_loc(nd_mig))
 allocate(n_migr(nd_mig,0:mpisize-1))
 n_migr_loc=0
 ind_fpart=0
 n_migr=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!---Compute the number of migrated migrated to the top or the bottom
!---Dectect the corresponding index
!............................................................................
 if (mpidims(2)==1) then ! decomposition=="1D"
    do ipart=1,npart_loc
        iy=0; iz=0
        if (nx==1) then
           iy = floor(fpart(ipart,2)/dy)
           iz = floor(fpart(ipart,3)/dz)
        elseif (nx>1) then
           iy = floor(fpart(ipart,3)/dy)
           iz = floor(fpart(ipart,4)/dz)
        endif

        if (iz>rb(3)) then     ! To the top
           n_migr_loc(1)=n_migr_loc(1)+1
           ind_fpart(1,n_migr_loc(1))=ipart

        elseif (iz<ra(3)) then ! To the bottom
           n_migr_loc(2)=n_migr_loc(2)+1
           ind_fpart(2,n_migr_loc(2))=ipart
        endif

    enddo

 elseif (mpidims(2)>1) then ! decomposition=="2D"
    do ipart=1,npart_loc
        iy = floor(fpart(ipart,3)/dy)
        iz = floor(fpart(ipart,4)/dz)

        if ( (iz>rb(3)).and.(iy.ge.ra(2)).and.(iy.le.rb(2)) ) then     ! Top
           n_migr_loc(1)=n_migr_loc(1)+1
           ind_fpart(1,n_migr_loc(1))=ipart
 
        elseif ( (iz<ra(3)).and.(iy.ge.ra(2)).and.(iy.le.rb(2)) ) then ! Bottom
           n_migr_loc(2)=n_migr_loc(2)+1
           ind_fpart(2,n_migr_loc(2))=ipart
        
        elseif ( (iy>rb(2)).and.(iz.ge.ra(3)).and.(iz.le.rb(3)) ) then ! Right
           n_migr_loc(3)=n_migr_loc(3)+1
           ind_fpart(3,n_migr_loc(3))=ipart
           
        elseif ( (iy<ra(2)).and.(iz.ge.ra(3)).and.(iz.le.rb(3)) ) then ! Left
           n_migr_loc(4)=n_migr_loc(4)+1
           ind_fpart(4,n_migr_loc(4))=ipart
        
        ! Corners
        elseif ( (iz>rb(3)).and.(iy>rb(2)) ) then ! Top Right
           n_migr_loc(5)=n_migr_loc(5)+1
           ind_fpart(5,n_migr_loc(5))=ipart
          
        elseif ( (iz>rb(3)).and.(iy<ra(2)) ) then ! Top Left
           n_migr_loc(6)=n_migr_loc(6)+1
           ind_fpart(6,n_migr_loc(6))=ipart
        
        elseif ( (iz<ra(3)).and.(iy>rb(2)) ) then ! Bottom Right
           n_migr_loc(7)=n_migr_loc(7)+1
           ind_fpart(7,n_migr_loc(7))=ipart

        elseif ( (iz<ra(3)).and.(iy<ra(2)) ) then ! Bottom Left
           n_migr_loc(8)=n_migr_loc(8)+1
           ind_fpart(8,n_migr_loc(8))=ipart
        endif

    enddo
 endif

 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Share the number of migrated particles between all the procs    
  call MPI_ALLGATHER (n_migr_loc, nd_mig , MPI_INTEGER, &
                      n_migr,nd_mig,  MPI_INTEGER, &
                      MPI_COMM_WORLD, mpicode)

 if ((maxval(n_migr).ne.0)) then
    if (mpirank==1000000) then
          print*,'n_migr',n_migr(:,mpirank)
          do i=1,npart_proc
             print*,'START ipart',i,'fpart start',fpart(i,1:4)
          enddo
    endif 
 endif
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Sort the data to send: the migrated particles are put on the end of fpart
 ! fpart(npart_proc,npart_proc-1,...,npart_proc-n_migr(1),npart_proc-n_migr(1)-1,...,npart_proc-n_migr(1)- n_migr(2))
 if ((mpisize.ne.1).and.(maxval(n_migr).ne.0)) then
    nskip=0
    allocate(value_temp(nd_mig,maxval(n_migr(:,mpirank)),1:n_lagr))
    do imigr=1,nd_mig
       !if ((n_migr(nd_mig,mpirank))<0) cycle
       if (imigr==1) then
          nskip=0
       else
          nskip=nskip+n_migr(imigr-1,mpirank)
       endif

       do i=n_migr(imigr,mpirank),1,-1
          fpart(npart_proc-i+1-nskip,:)=fpart(ind_fpart(imigr,i),:)
          value_temp(imigr,i,:)=fpart(ind_fpart(imigr,i),:)
          fpart(ind_fpart(imigr,i),:)=1.d16
       enddo
    enddo
    !
    ! Sort fpart for increasing index
    allocate(A(npart_loc),order(npart_loc),temp_sort(npart_loc,1:n_lagr))
    temp_sort=fpart(1:npart_loc,:)
    do i=1,npart_loc
        A(i)%value = fpart(i,1)
        A(i)%order = i
    enddo
    call QSort(A,npart_loc)
    order=A%order     
    do i=1,npart_loc
        fpart(i,:)=temp_sort(order(i),:)
    enddo
    deallocate(A,order,temp_sort)
 
    nskip=0
    do imigr=1,nd_mig
       if (imigr==1) then
          nskip=0
       else
          nskip=nskip+n_migr(imigr-1,mpirank)
       endif    
       do i=n_migr(imigr,mpirank),1,-1                
          fpart(npart_proc-i+1-nskip,:)=value_temp(imigr,i,:)
       enddo
        
      npart_loc=npart_loc-n_migr_loc(imigr) 

    enddo
   deallocate(value_temp)

 endif

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 coord=0
 call MPI_Cart_coords(mpicommcart,mpirank,2,coord,mpicode)

 if ((mpisize>1).and.(maxval(n_migr).ne.0)) then
    nskip=0
    iorigin=mpirank
    id_migr=0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Send the migrated particles to the corresponding procs
    do imigr=1,nd_mig     
       if (imigr==1) then
          nskip=0
       else
          nskip=nskip+n_migr(imigr-1,mpirank)
       endif

       id=0
       if (imigr==1) then
          id=(/+1, 0/) ! top
          id_migr=2
       elseif (imigr==2) then
          id=(/-1, 0/) ! Bottom
          id_migr=1
       elseif (imigr==3) then 
          id=(/0, +1/) ! Right
          id_migr=4
       elseif (imigr==4) then
          id=(/0, -1/) ! Left
          id_migr=3
       elseif (imigr==5) then
          id=(/+1, +1/) ! Top right
          id_migr=8
       elseif (imigr==6) then
          id=(/+1, -1/) ! Top left
          id_migr=7
       elseif (imigr==7) then
          id=(/-1, +1/) ! Bottom right
          id_migr=6
       elseif (imigr==8) then
          id=(/-1, -1/) ! Bottom Left
          id_migr=5
       endif

       cordout(1)=coord(1)+id(1)
       cordout(2)=coord(2)+id(2)  
       cordout(1)=per(cordout(1),mpidims(1))
       cordout(2)=per(cordout(2),mpidims(2))
       call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)

       if (n_migr(imigr,iorigin)>0) call MPI_send( &
           fpart(npart_proc+1-n_migr(imigr,iorigin)-nskip:npart_proc-nskip,:),&! buffer
           n_lagr*n_migr(imigr,iorigin),&                 ! buffer size
           MPI_DOUBLE_PRECISION,&                         ! datatype
           idestination,&                                 ! idestination
           idestination,&                                 ! tag
           MPI_COMM_WORLD,mpicode) !MPI_COMM_WORLD,statu,mpicode)
    enddo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Receive the migrated particles to the corresponding procs     
    nskip=0
    id_migr=0
    do imigr=1,nd_mig     
       if (imigr==1) then
          nskip=0
       else
          nskip=nskip+n_migr(imigr-1,mpirank)
       endif

       id=0
       if (imigr==1) then
          id=(/+1, 0/) ! top
          id_migr=2
       elseif (imigr==2) then
          id=(/-1, 0/) ! Bottom
          id_migr=1
       elseif (imigr==3) then 
          id=(/0, +1/) ! Right
          id_migr=4
       elseif (imigr==4) then
          id=(/0, -1/) ! Left
          id_migr=3
       elseif (imigr==5) then
          id=(/+1, +1/) ! Top right
          id_migr=8
       elseif (imigr==6) then
          id=(/+1, -1/) ! Top left
          id_migr=7
       elseif (imigr==7) then
          id=(/-1, +1/) ! Bottom right
          id_migr=6
       elseif (imigr==8) then
          id=(/-1, -1/) ! Bottom Left
          id_migr=5
       endif

       cordout(1)=coord(1)+id(1)
       cordout(2)=coord(2)+id(2)  
       cordout(1)=per(cordout(1),mpidims(1))
       cordout(2)=per(cordout(2),mpidims(2))
       call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
       if (n_migr(id_migr,idestination)>0) then 
           call MPI_recv( &
           fpart(npart_loc+1:npart_loc+n_migr(id_migr,idestination),:),&! buffer         
           n_lagr*n_migr(id_migr,idestination),&                ! buffer size
           MPI_DOUBLE_PRECISION,&                               ! datatype
           idestination,&                                       ! source
           iorigin,&                                            ! tag
           MPI_COMM_WORLD,statu,mpicode)           
       ! 
       fpart(npart_loc+n_migr(id_migr,idestination)+1:npart_proc,:)=1.d16
       npart_loc=npart_loc+n_migr(id_migr,idestination) 
       endif
       !
       ! Sort fpart for increasing index
       allocate(A(npart_loc),order(npart_loc),temp_sort(npart_loc,1:n_lagr))
       temp_sort=fpart(1:npart_loc,:)
       do i=1,npart_loc
          A(i)%value = fpart(i,1)
          A(i)%order = i
       enddo
       call QSort(A,npart_loc)
       order=A%order     
       do i=1,npart_loc
          fpart(i,:)=temp_sort(order(i),:)
       enddo
       deallocate(A,order,temp_sort)
    enddo

 endif

!
!............................................................................  
call toc("Particles (migration)", MPI_wtime()-t_migr)

end subroutine migration_part
