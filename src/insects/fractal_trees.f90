subroutine draw_fractal_tree(Insect, mask, mask_color, us)
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: io_error, nlines, mpicode,i
  real (kind=pr) :: R, N_smooth, safety, x1(1:3),x2(1:3), t1
  character(len=2048) :: dummy
  character(len=strlen) :: file
  real(kind=pr),allocatable,dimension(:,:) :: treedata

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  if (smoothing == 0.d0) then
    N_smooth = 1.5d0
    ! thickness of smoothing layer and safety distance
    smoothing = N_smooth*max(dx,dy,dz)
    safety = 2.d0*N_smooth*max(dx,dy,dz)
    Insect%safety = safety
  endif

  file = 'tree_data.in'
  call check_file_exists( file )

  ! error checks
  if (iMoving == 1) then
    call abort(4417,"fractal trees do not move -- set iMoving=0 (avoid creating mask in every iteration)")
  endif

  !*****************************************************************************
  ! phase one: read the number of lines (which is the number of rigid cylinders in the tree)
  !*****************************************************************************
  call count_lines_in_ascii_file_mpi(file, nlines, 0)
  if (root) then
    write(*,'("Building a fractal tree with ",i5," rigid cylinders")') nlines
  endif

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
  t1 = MPI_Wtime()
 do i=1, nlines
    x1 = treedata(i,1:3) + (/x0,y0,z0/)
    x2 = treedata(i,4:6) + (/x0,y0,z0/)
    ! the file containes the radius of the cylinder
    R = treedata(i,7)
    call draw_cylinder_new( x1, x2, R, mask, mask_color, us, Insect, int(1,kind=2))
  end do

  R = MPI_wtime() - t1
  R = mpisum(R)
  if (root) then
    write(*,'(80("-"))')
    write(*,'("done creating fractal tree mask of ",i4," branches")') nlines
    write(*,'("wtime on master process is ",es12.4," secs")') MPI_wtime()-t1
    write(*,'("sum wtime on all cpu ",es12.4," secs (avg=",es12.4,")")') R, R/dble(mpisize)
    write(*,'(80("-"))')
  end if

  deallocate(treedata)
end subroutine draw_fractal_tree
