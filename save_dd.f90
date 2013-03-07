subroutine save_dd (time)
!====================================================================
!     Save domain decomposition into a file
!====================================================================
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  real (kind=pr), intent (in) :: time
  integer :: mpicode, irank, j
  integer, dimension (18) :: blocksizes
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus
  character (len=17) :: name

  !--Set up file name base
  write (name, '(es10.4)') time  

  if ( mpirank == 0 ) then

     open (15, file = 'domain_dec_'//name, status = 'replace')
     write (15, '(19(1x,i6))') 0, (ra(j), j=1,3), (rb(j), j=1,3), (rs(j), j=1,3), &
                                 (ca(j), j=1,3), (cb(j), j=1,3), (cs(j), j=1,3)
     do irank = 1,mpisize-1
        call MPI_RECV (blocksizes, 18, mpiinteger, irank, 101, MPI_COMM_WORLD, mpistatus, mpicode)
        write (15, '(19(1x,i6))') irank, (blocksizes(j), j=1,18)
     enddo
     close (15)

  else

     blocksizes(1:3) = ra(1:3)
     blocksizes(4:6) = rb(1:3)
     blocksizes(7:9) = rs(1:3)
     blocksizes(10:12) = ca(1:3)
     blocksizes(13:15) = cb(1:3)
     blocksizes(16:18) = cs(1:3)
     call MPI_SEND (blocksizes, 18, mpiinteger, 0, 101, MPI_COMM_WORLD, mpicode)

  endif

end subroutine save_dd

