! Spherical obstacle
subroutine draw_fractal_tree()
  use mpi
  use vars
  use penalization ! mask array etc
  use insect_module
  implicit none

  integer :: io_error, nlines, mpicode,i
  real (kind=pr) :: R, N_smooth, safety, x1(1:3),x2(1:3)
  character(len=2048) :: dummy
  character(len=strlen) :: file
  type(diptera) :: dummyinsect
  real(kind=pr),allocatable,dimension(:,:) :: treedata

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  N_smooth = 1.5d0
  ! thickness of smoothing layer and safety distance
  smoothing = N_smooth*max(dx,dy,dz)
  safety = 2.d0*N_smooth*max(dx,dy,dz)
  dummyinsect%safety = safety

  file = 'tree_data.in'

  call check_file_exists( file )
  !*****************************************************************************
  ! phase one: read the number of lines (which is the number of rigid cylinders in the tree)
  !*****************************************************************************
  if (root) then
    io_error = 0
    nlines = 0
    open(unit=14,file=file,action='read',status='old')
    do while (io_error==0)
      read (14,'(A)',iostat=io_error) dummy
      if (io_error==0) then
        nlines = nlines+1
      endif
    enddo
    close (14)
    write(*,'("Building a fractal tree with ",i5," rigid cylinders")') nlines
  endif

  call MPI_BCAST(nlines,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode)

  !*****************************************************************************
  ! phase two: read all cylinders from the file into the array and bcast them
  !*****************************************************************************
  allocate( treedata(1:nlines, 1:7) )

  if (root) then
    open(unit=14,file=file,action='read',status='old')
    do i=1,nlines
      read (14,'(A)',iostat=io_error) dummy
      if (io_error==0) then
        read (dummy,*) treedata(i,:)
        write(*,'("read cylinder ",7(es12.4,1x))') treedata(i,:)
      endif
    enddo
    close (14)
  endif

  call MPI_BCAST(treedata,nlines*7,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)

  !*****************************************************************************
  ! phase 3: all ranks draw the individual cylinders...
  ! please note we assume the root point to be 0 0 0 in the file
  !*****************************************************************************
  do i=1, nlines
    x1 = treedata(i,1:3) + (/x0,y0,z0/)
    x2 = treedata(i,4:6) + (/x0,y0,z0/)
    ! the file containes the radius of the cylinder
    R = treedata(i,7)
    call draw_cylinder_new( x1, x2, R, mask, mask_color, us, dummyinsect, int(1,kind=2))
  enddo

  deallocate(treedata)
end subroutine draw_fractal_tree
