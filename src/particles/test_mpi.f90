subroutine test_mpi()
 ! Benjamin 01 July 2019   
 use vars
 use mpi
 use qsort_mod

 use p3dfft
 use fftw3_descriptors

 implicit none

 integer :: mpicode,idestination, iorigin, statu(MPI_status_size)
 integer :: disp,dir,source,dest
 integer :: coord(2), cordout(2), nd_mig, imigr, id(2)
 integer :: nmpidims=2

  ng=2

  ! default case, no decomposition
  decomposition="none"

  if (mpirank==0) then
    write(*,'(A)') "-----------------------------------p3dfft init---------------------------"
    write(*,'("Initializing P3DFFT, n=(",3(i4,1x),")")') nx,ny,nz
    write(*,'("P3DFFT reserves about ",i8,"MB (",i4,"GB) for internal work arrays")') &
    nint(1.6d-5*dble(nx)*dble(ny)*dble(nz)), nint(1.6d-5*dble(nx)*dble(ny)*dble(nz)/1000.d0)
    write(*,'("P3DFFT domain size:     ",3(g15.8,1x))') xl,yl,zl
    write(*,'("P3DFFT scale factors:   ",3(g15.8,1x))') scalex, scaley, scalez
    write(*,'("P3DFFT lattice spacing: ",3(g15.8,1x))') dx,dy,dz


    print*,' '
    print*,'comp',dble(nz/mpisize),1+ng
  endif

  !-- Set up dimensions. It is very important that mpidims(2) > mpidims(1)
  ! because P3Dfft crashes otherwise. This means a 1D decomposition is always
  ! along the z direction in real space.
  if (dble(nz)/dble(mpisize)>1+ng) then
     mpidims(1) = 1             ! due to p3dfft, 1D decomposition is always the
     mpidims(2) = mpisize       ! 3rd index in real space.
     decomposition="1D"
  else
     mpidims(1) = 0
     mpidims(2) = 0
     call MPI_Dims_create(mpisize,nmpidims,mpidims,mpicode)
     if(mpidims(1) > mpidims(2)) then
        mpidims(1) = mpidims(2)
        mpidims(2) = mpisize / mpidims(1)
     endif

     if (dble(nz)/dble(mpidims(2))<1+ng)then
        mpidims(2) = nz/3
        mpidims(1) = mpisize / mpidims(2)
     endif


     !mpidims(1) = nz/3
     !mpidims(2) = mpisize / mpidims(1) 

     decomposition="2D"
  endif

  if (root) write(*,'("mpidims= ",i3,1x,i3)') mpidims
  if (root) write(*,'("ng= ",i3)') ng
  if (root) write(*,'("Using ",A," decomposition!")') trim(adjustl(decomposition))

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
    call abort(33340,"wrong mpidims: change mpisize")
  endif

  !-- Check dimensions order
  if(mpidims(1)>mpidims(2)) then
    call abort(33341,"mpidims(1)>mpidims(2): max mpisize")
  endif


  !-- Initialize P3DFFT
  call p3dfft_setup(mpidims,nx,ny,nz,MPI_COMM_WORLD, overwrite=.false.)

  !-- Get Cartesian topology info
  call p3dfft_get_mpi_info(mpitaskid,mpitasks,mpicommcart)

  !-- Get local sizes
  call p3dfft_get_dims(ra,rb,rs,1)  ! real blocks
  call p3dfft_get_dims(ca,cb,cs,2)  ! complex blocks
  ra(:) = ra(:) - 1
  rb(:) = rb(:) - 1
  ca(:) = ca(:) - 1
  cb(:) = cb(:) - 1

  !-- extends of real arrays that have ghost points. We add ghosts in all
  !-- directions, including the periodic ones.
  ga=ra-ng
  gb=rb+ng

  if (nx==1) then
    ga(1)=0
    gb(1)=0
  endif

  if ( rb(2)-ra(2)+1<(ng+1) .or. rb(3)-ra(3)+1<(ng+1) ) then
    if (mpirank==0) write(*,*) "Too many CPUs: the ghosts span more than one CPU"
    if (mpirank==0) write(*,*) "y", rb(2)-ra(2)+1, "z", rb(3)-ra(3)+1
    ! Benjamin Comment the following lines for particles part
    call abort(33341, "Too many CPUs: the ghosts span more than one CPU")
 call MPI_FINALIZE(mpicode)
 call exit(0)
  endif

  call print_domain_decomposition_test()


 !call print_domain_decomposition_test()
! Attention mpicommcart should be periodic

 if (mpirank==0) print*,'START TEST MPI'
 if (mpirank==0) print*,'mpidims',mpidims
 !print*,'mpirank',mpirank,'y',ra(2),rb(2),'z',ra(3),rb(3)
 


 if ((mpirank==-10) .and. (mpidims(2)==1)) then
    ! mpexex -np 4
    call MPI_Cart_coords(mpicommcart,mpirank,2,coord,mpicode)
    print*,'At mpirank=',mpirank,' my coords are ',coord

    ! TOP
    cordout(1)=coord(1)+1
    cordout(2)=coord(2)  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Top   : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! BOTTOM
    cordout(1)=coord(1)-1
    cordout(2)=coord(2)  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Bottom: At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

endif

if ((mpirank==-10) .and. (mpidims(2)>1)) then
    ! mpexex -np 6
    call MPI_Cart_coords(mpicommcart,mpirank,2,coord,mpicode)
    print*,'At mpirank=',mpirank,' my coords are ',coord

    ! TOP
    cordout(1)=coord(1)+1
    cordout(2)=coord(2)  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Top           : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! BOTTOM
    cordout(1)=coord(1)-1
    cordout(2)=coord(2)  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Bottom        : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! RIGHT
    cordout(1)=coord(1)
    cordout(2)=coord(2)+1  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Right         : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! LEFT
    cordout(1)=coord(1)
    cordout(2)=coord(2)-1  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Left          : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! TOP RIGHT
    cordout(1)=coord(1)+1
    cordout(2)=coord(2)+1
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Top Right     : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! TOP LEFT
    cordout(1)=coord(1)+1
    cordout(2)=coord(2)-1 
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Top Left      : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! BOTTOM RIGHT
    cordout(1)=coord(1)-1
    cordout(2)=coord(2)+1  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Bottom right  : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination

    ! BOTTOM LEFT
    cordout(1)=coord(1)-1
    cordout(2)=coord(2)-1  
    cordout(1)=per(cordout(1),mpidims(1))
    cordout(2)=per(cordout(2),mpidims(2))
    call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
    print*,'Bottom Left   : At mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination
endif

 

if (mpirank==mpisize/2) then
    call MPI_Cart_coords(mpicommcart,mpirank,2,coord,mpicode)
    nd_mig=2
    if (decomposition=="2D") nd_mig=8
    print*,'At mpirank=',mpirank,' my coords are ',coord,' nd_mig=',nd_mig

    do imigr=1,nd_mig
       id=0
       if (imigr==1) id=(/+1, 0/) ! top
       if (imigr==2) id=(/-1, 0/) ! Bottom
       if (imigr==3) id=(/0, +1/) ! Right
       if (imigr==4) id=(/0, -1/) ! Left

       if (imigr==5) id=(/+1, +1/) ! Top right
       if (imigr==6) id=(/+1, -1/) ! Top left
       if (imigr==7) id=(/-1, +1/) ! Bottom right
       if (imigr==8) id=(/-1, -1/) ! Bottom Left
        
    
       cordout(1)=coord(1)+id(1)
       cordout(2)=coord(2)+id(2)  
       cordout(1)=per(cordout(1),mpidims(2))
       cordout(2)=per(cordout(2),mpidims(1))
       call MPI_Cart_Rank (mpicommcart,cordout,idestination,mpicode)
       print*,'imigr=',imigr,' at mpirank=',mpirank, ' coords=',cordout,' the rank is ',idestination
    enddo
endif



 call MPI_FINALIZE(mpicode)
 call exit(0)
end subroutine test_mpi

subroutine print_domain_decomposition_test()
  use vars
  use mpi
  implicit none
  integer :: mpicode
  
  
  if (root) then
     write(*,'(A)') '--------------------------------------'
     write(*,'(A)') '*** Domain decomposition:'
     write(*,'(A)') '--------------------------------------'
  endif
  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("mpirank=",i5," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
       &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
       mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)
  if (root) then
     write(*,'(A)') '--------------------------------------'
  endif
end subroutine 

