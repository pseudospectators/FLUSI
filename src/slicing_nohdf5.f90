module slicing
  use vars
  implicit none

  ! module global variables
  integer, parameter :: nslices = 4
  integer, save :: icache = 0, iflush = 0, ncache = 200 !(set in params file)
  real(kind=pr),allocatable, save :: slice_cache(:,:,:,:,:), times_cache(:)



  contains


! save slices to cache, and flush when the cache is full. The module allows to
! save up to nslices yz-cuts with x=const. They are first stored in a module-global
! cache array, where the x-coordinate is the time. When the cache is full, the
! field is flushed into a single HDF5 file. I/O is in parallel, there is no
! communication cost. Since the HDF5 file does not contain the times at which
! the slices are, we save that separately to an ascii file. The cache is also
! flushed when slice_free is called. If the cache is half full and flush_slices
! is called, then we only save half of the cache (no redundancy)
subroutine save_slices( time, u )
  implicit none
  real(kind=pr),intent(in) :: time
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr) :: t1
  integer :: slice_slot

  t1 = MPI_wtime()

  ! loop over slices (we save up to four different ones)
  do slice_slot = 1, nslices
    ! save only if the index is positive (use negative ones to disable)
    if (slices_to_save(slice_slot) >= 0) then
      ! save slice in cache
      slice_cache( icache, :, :, :, slice_slot) = u( slices_to_save(slice_slot) ,:,:,: )
      times_cache( icache ) = time
    endif
  enddo


  ! if the cache is full, then we flush it to disk and start over
  if (icache == ncache) then
      ! the cache is full
      call flush_slices
  else
      ! the cache is not full
      icache = icache + 1
  endif

  tslices = tslices + MPI_wtime() - t1

end subroutine save_slices



! flush slices. Saves the content of the slice cache (module-global) to an HDF5
! file. The x-dimension of that file is the time
! note we save only until icache, so only the
! updated values in the cache, and not the entire array, since entries
! greater then icache are from earlier times
! Thus: If the cache is half full and flush_slices
! is called, then we only save half of the cache (no redundancy)
subroutine flush_slices
  implicit none

  !!! TODO: NOT IMPLEMENTED !!!
  if (mpirank==0) write(*,'(A)',advance='no') "WARNING(flush_slices): Function not implemented without HDF5 support.."

end subroutine flush_slices





! initialize the slicing module.
subroutine slice_init(time)
  implicit none
  real(kind=pr),intent(inout)::time
  integer :: i

  if (mpirank==0) write(*,'(A)',advance='no') "Initializing slicing module.."

  ! overwrite the time-log file (note the time is not included in the HDF5 file)
  if ((mpirank==0).and.(inicond(1:8).ne."backup::")) then
    call init_empty_file('slice_times.t')
  endif

  iflush = nint( time * 1000.d0 )

  if (mpirank==0) then
    do i=1,nslices
      if (slices_to_save(i)>=nx) then
        write(*,*) "Error: slicing index", slices_to_save(i), "out of bounds. stop."
        call abort( 9000, "Error slicing index out of bounds.")
      endif
    enddo
  endif

  ! copy value from global parameters vars
  ncache = ncache_slices

  ! ensure the slice cache is allocated
  if (.not.allocated(slice_cache) ) then
    allocate ( slice_cache(0:ncache-1, ra(2):rb(2), ra(3):rb(3), 1:nd, 1:nslices) )
  endif

  if (.not.allocated(times_cache) ) then
    allocate ( times_cache(0:ncache-1) )
  endif
  if (mpirank==0) write(*,*) "..DONE!"
end subroutine slice_init





! free the memory occupied by slicing and flush the cache to disk, in case
! this has not yet been done.
subroutine slice_free
  implicit none

  call flush_slices

  if (allocated(slice_cache)) deallocate( slice_cache )
  if (allocated(times_cache)) deallocate( times_cache )
end subroutine slice_free



! collect a yz-slice of data from all procs on the root rank.
subroutine gather_slice_yz( u, slice, ixslice )
  implicit none
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout) :: slice(0:,0:)
  integer, intent(in) :: ixslice

  integer :: mpierror,newtype,extent2
  integer :: resizedtype, ierr, dims, iy
  integer(8) :: begin, extent1
  integer, dimension(2) :: sizes, subsizes, starts
  integer, allocatable :: displs(:), counts(:)
  real(kind=pr), allocatable ::local(:,:)

  !-----------------------------------------------------------------------------
  ! allocate storage for the local slice. This way we are sure about it's memory
  ! layout: it is FORTRAN column-major ordering, and the second index is slowest.
  !-----------------------------------------------------------------------------
  allocate( local(ra(2):rb(2),ra(3):rb(3)) )
  local = u(ixslice,:,:)


  !-----------------------------------------------------------------------------
  ! create receiving datatype, relevant only on root. This datatype in non
  ! contiguous in memory. A very good description of the problem is found in
  ! http://stackoverflow.com/questions/17508647/sending-2d-arrays-in-fortran-with-mpi-gather
  !-----------------------------------------------------------------------------
  ! size of global array
  sizes    = (/ny, nz/)
  ! size of local subarrays (we take the biggest one, in case there's different)
  subsizes = (/maxval((rb_table(2,:)-ra_table(2,:)+1)),maxval((rb_table(3,:)-ra_table(3,:)+1))/)
  ! no offset
  starts   = (/0, 0/)
  dims     = 2 ! 2D array
  call MPI_Type_create_subarray(dims, sizes, subsizes, starts,            &
                                MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION,  &
                                newtype, mpierror)
  ! get the sice of a double
  call MPI_Type_size(MPI_DOUBLE_PRECISION, extent2, ierr)
  extent1 = extent2
  begin    = 0
  ! resize the datatype to the size of one double
  call MPI_Type_create_resized(newtype, begin, extent1, resizedtype, ierr)
  call MPI_Type_commit(resizedtype, ierr)



  !-----------------------------------------------------------------------------
  ! root has to know how much data it is going to receive (counts) and where to
  ! put it in a 1D array in memory (displs). For each proc, we allow a different
  ! count (and displacement of course)
  !-----------------------------------------------------------------------------
  allocate ( counts(0:mpisize-1),displs(0:mpisize-1) )
  do iy = 0, mpisize-1
    ! number of elements in local memory
    counts(iy) = (rb_table(2,iy)-ra_table(2,iy)+1)*(rb_table(3,iy)-ra_table(3,iy)+1)
    ! displacement on the 1D line on root
    displs(iy) = ra_table(2,iy) + (nz)*(ra_table(3,iy))
  enddo



  !-----------------------------------------------------------------------------
  ! Note different datatypes for sender and receiver. The sender can simple use
  ! MPI_DOUBLE_PRECISION, since it is sending the entire slice, which is locally
  ! contiguous in memory.
  ! The situation on the receiver, root, is very different: the 2D data is not
  ! at all contiguous, but it has holes in it, which is why we created the data-
  ! type "resizedtype"
  !-----------------------------------------------------------------------------
  call MPI_Gatherv( local, (rb(2)-ra(2)+1)*(rb(3)-ra(3)+1), MPI_DOUBLE_PRECISION, &
                    slice, counts, displs, resizedtype, &
                    0, MPI_COMM_WORLD, ierr)

  ! cleanup
  deallocate( counts, displs, local )
end subroutine gather_slice_yz








subroutine test_slices( u )
  implicit none
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),allocatable :: slice(:,:)
  integer :: ix,iy,iz,mpierror
  real(kind=pr) ::x,y,z


  u = -17.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x = dble(ix)*dx
        y = dble(iy)*dy
        z = dble(iz)*dz
        u(ix,iy,iz,1) = y+z
      enddo
    enddo
  enddo


  if (mpirank==0) then
    allocate( slice(0:ny-1,0:nz-1) )
  else
    allocate( slice(0:0,0:0) )
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  call gather_slice_yz( u(:,:,:,1), slice, 1 )
  !!!!!!!!!!!!!!!!!!!!!!!!!!!


  if (mpirank==0) then
    do iz = 0, nz-1
      write (*,'(18(f4.1,1x))') slice(:,iz)
    enddo

    do iz=0,nz-1
    do iy=0,ny-1
        y = dble(iy)*dy
        z = dble(iz)*dz
        if (abs(slice(iy,iz)-(y+z)) > 1.0e-11 ) write(*,*) "ERROR"
    enddo
    enddo

  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!


  call MPI_barrier (MPI_COMM_world, mpierror)
  call abort(0, "ByeBye as planned!")


end subroutine




end module slicing
