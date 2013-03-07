! ---------------------------------------------------------------
! MPI three-dimensional array transpose
! ---------------------------------------------------------------
!
! idir  - direction
!    idir = 1 - permute 1 and 2 indices
!    idir = 2 - permute 1 and 3 indices 
! n  - sizes
! ka, kb, ks  - input extents and sizes
!    ka(1) = 0
!    kb(1) = n(1)-1
!    ka(2) = mpicoords(2)*n(2)/mpidims(2)
!    kb(2) = (mpicoords(2)+1)*n(2)/mpidims(2)-1
!    ka(3) = mpicoords(1)*n(3)/mpidims(1)
!    kb(3) = (mpicoords(1)+1)*n(3)/mpidims(1)-1
!    ks(:) = kb(:)-ka(:)+1
! kat, kbt, kst  - output extents and sizes
!    kat(1) = 0
!    kbt(1) = n(idir+1)-1
!    kat(1+idir) = mpicoords(3-idir)*n(1)/mpidims(3-idir)
!    kbt(1+idir) = (mpicoords(3-idir)+1)*n(1)/mpidims(3-idir)-1
!    kat(4-idir) = ka(4-idir)
!    kbt(4-idir) = kb(4-idir)
!    kst(:) = kbt(:)-kat(:)+1
! mpidims  - sizes of MPI Cartesian topology
! mpicoords  - coordimates in the MPI Cartesian topology
! mpicommslab  - communicators  of the MPI Cartesian sub-topologies
! mpireal  - MPI real type, corresponding to datain and dataout
! datain  - input array
! dataout  - output array
!
subroutine subtr ( idir, n, ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, mpireal, datain, dataout )
  use mpi_header ! Module incapsulates mpif.
  implicit none

  integer, intent (in) :: idir, mpireal
  integer, dimension (1:3), intent (in) :: n, ka, kb, ks, kat, kbt, kst
  integer, dimension (2), intent (in) :: mpidims, mpicoords, mpicommslab
  real (kind=kind(0.0)), dimension (ka(1):kb(1),ka(2):kb(2),ka(3):kb(3)), intent (in) :: datain
  real (kind=kind(0.0)), dimension (kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)), intent (out) :: dataout
  integer :: jj, kk
  integer :: mpicode, mpirealsize, mpiblock
  integer, dimension(:), allocatable :: types, counts
  integer, dimension(:), allocatable :: displs
  logical :: reorganization
  real (kind=kind(0.0)), dimension (:,:,:), allocatable :: datatemp

  !--Check direction
  if ( (idir/=1) .and. (idir/=2) ) then
     print *, 'tr2dplus: illegal idir'
     stop
  endif
  
  write (*,*) "is this used?"

  !--Allocate temporary array
  allocate ( datatemp(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )

  !--Transpose blocks locally
  if ( idir == 1 ) then
     do kk = ka(3), kb(3)
        do jj = 0, mpidims(2)-1
           dataout((jj*ks(2)):((jj+1)*ks(2)-1),:,kk) = transpose (datain((jj*kst(2)):((jj+1)*kst(2)-1),:,kk))
        enddo
     enddo
  else   ! if idir == 2
     do kk = ka(2), kb(2)
        do jj = 0, mpidims(1)-1
           dataout((jj*ks(3)):((jj+1)*ks(3)-1),kk,:) = transpose (datain((jj*kst(3)):((jj+1)*kst(3)-1),kk,:))
        enddo
     enddo
  endif

  !--Communications between processes are required when transpose plans are non-local 
  if ( mpidims(3-idir) > 1 ) then

     !--Copy locally transposed blocks 
     datatemp(:,:,:) = dataout(:,:,:)  

     !--Create type for communication
     call MPI_TYPE_VECTOR (kst(2)*kst(3),ks(1+idir),kst(1),mpireal,mpiblock,mpicode)
     call MPI_TYPE_COMMIT (mpiblock,mpicode)

     !--Compute sizes and displacements
     allocate (counts(mpidims(3-idir)),displs(mpidims(3-idir)),types(mpidims(3-idir)))
     call MPI_TYPE_SIZE (mpireal,mpirealsize,mpicode)
     counts(:) = 1
     displs(1) = 0
     do jj = 2, mpidims(3-idir)
        displs(jj) = displs(jj-1) + ks(idir+1)*mpirealsize
     enddo
     types(:) = mpiblock

     !--Communicate
     call MPI_ALLTOALLW (datatemp,counts,displs,types,dataout,counts,displs,types,mpicommslab(idir),mpicode)

     !--Deallocate types, sizes and displacements
     call MPI_TYPE_FREE(mpiblock,mpicode)
     deallocate (counts,displs,types)

  endif

  !--Deallocate temporary array
  deallocate (datatemp)

end subroutine subtr

