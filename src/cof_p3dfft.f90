!====================================================================
!====================================================================
!
!     Fourier transform subroutines using P3DFFT and FFTW3
!
!====================================================================
!====================================================================

module fftw3_descriptors
  use vars
  implicit none
  integer*8,dimension(1:3),save :: Desc_Handle_1D_f,Desc_Handle_1D_b

end module fftw3_descriptors

!---------------------------------------------------------------------
! module that contains wrappers for P3DFFT
!---------------------------------------------------------------------

module p3dfft_wrapper
  use vars
  implicit none
  
  contains

! Compute the FFT of the real-valued 3D array inx and save the output
! in the complex-valued 3D array outk.
subroutine fft(outk,inx)
    use mpi
    use vars ! For precision specficiation and array sizes    
    real(kind=pr),intent(in)::inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

    call coftxyz(inx,outk)
end subroutine fft


! Compute the inverse FFT of the complex-valued 3D array ink and save the
! output in the real-valued 3D array outx.
subroutine ifft(outx,ink)
    use mpi
    use vars ! For precision specficiation and array sizes    
    complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
    real(kind=pr),intent(out)::outx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

    call cofitxyz(ink,outx)
end subroutine ifft


! Compute the FFT of the real-valued 3D array inx and save the output
! in the complex-valued 3D array outk.
subroutine fft3(outk,inx)
    use mpi
    use vars ! For precision specficiation and array sizes    
    real(kind=pr),intent(in)::inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
    complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
    call coftxyz(inx(:,:,:,1),outk(:,:,:,1))
    call coftxyz(inx(:,:,:,2),outk(:,:,:,2))
    call coftxyz(inx(:,:,:,3),outk(:,:,:,3))
end subroutine fft3


! Compute the inverse FFT of the complex-valued 3D array ink and save the
! output in the real-valued 3D array outx.
subroutine ifft3(outx,ink)
    use mpi
    use vars ! For precision specficiation and array sizes    
    complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
    real(kind=pr),intent(out)::outx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
    call cofitxyz(ink(:,:,:,1),outx(:,:,:,1))
    call cofitxyz(ink(:,:,:,2),outx(:,:,:,2))
    call cofitxyz(ink(:,:,:,3),outx(:,:,:,3))    
end subroutine ifft3

subroutine fft_initialize
  !-----------------------------------------------------------------------------
  !     Allocate memory and initialize FFT
  !-----------------------------------------------------------------------------
  ! This subroutine performs two major tasks:
  !     * Initialize P3DFFT and provide the memory management (how big are the 
  !       chunks of data on each process (coftxyz and cofitxyz)
  !     * Initialize 1D transforms that rely on FFTW (cofts and cofits)
  !       the auxiliary subroutine textents is used here
  !-----------------------------------------------------------------------------
  use mpi ! Module incapsulates mpif.
  use vars
  use p3dfft
  use fftw3_descriptors
  implicit none
  include 'fftw3.f'

  integer,parameter :: nmpidims = 2
  integer :: mpicode,idir,L,n
  integer,dimension(1:3) :: ka,kb,ks,kat,kbt,kst
  logical,dimension(2) :: subcart
  real(kind=pr),dimension(:,:),allocatable :: f,ft
  integer, dimension(:,:), allocatable :: yz_plane_local

  !-----------------------------------------------------------------------------
  !------ Three-dimensional FFT                                           ------
  !-----------------------------------------------------------------------------
  ! default case, no decomposition
  decomposition="none"
  
  !-- Set up dimensions. It is very important that mpidims(2) > mpidims(1)
  ! because P3Dfft crashes otherwise. This means a 1D decomposition is always
  ! along the z direction in real space. 
  if (nz>mpisize.and.nx==1) then
     mpidims(1) = 1             ! due to p3dfft, 1D decomposition is always the
     mpidims(2) = mpisize       ! 3rd index in real space.
     decomposition="1D"
  else
     ! unfortunately, 2D data decomposition does not work with 2D code (nx==1)
     mpidims = 0
     call MPI_Dims_create(mpisize,nmpidims,mpidims,mpicode)
     if(mpidims(1) > mpidims(2)) then
        mpidims(1) = mpidims(2)
        mpidims(2) = mpisize / mpidims(1)
     endif
     decomposition="2D"
  endif
  
  if (root) write(*,'("mpidims= ",i3,1x,i3)') mpidims
  if (root) write(*,'("Using ",A," decomposition!")') trim(adjustl(decomposition))

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
     print *, 'wrong mpidims: change mpisize'
     call abort()
  endif

  !-- Initialize P3DFFT
  call p3dfft_setup(mpidims,nx,ny,nz,MPI_COMM_WORLD, overwrite=.false.)

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
  
  if ( rb(2)-ra(2)+1<2*ng .or. rb(3)-ra(3)+1<2*ng ) then
    if (mpirank==0) write(*,*) "Too many CPUs: the ghosts span more than one CPU"
    if (mpirank==0) write(*,*) "y", rb(2)-ra(2)+1, "z", rb(3)-ra(3)+1
    call abort()
  endif
  
  
  !-- Allocate domain partitioning tables and gather sizes from all processes 
  !-- (only for real arrays)
  allocate ( ra_table(1:3,0:mpisize-1), rb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER (ra, 3, MPI_INTEGER, ra_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  call MPI_ALLGATHER (rb, 3, MPI_INTEGER, rb_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)

  !-- create a 2D array of size ny*nz that holds the mpirank of every point
  !-- in the y-z plane
  allocate (yz_plane_ranks(0:ny-1,0:nz-1) )
  allocate (yz_plane_local(0:ny-1,0:nz-1) )
  
  yz_plane_local = 0
  yz_plane_local( ra(2):rb(2), ra(3):rb(3) ) = mpirank
  
  call MPI_ALLREDUCE ( yz_plane_local,yz_plane_ranks,ny*nz,&
                      MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)
  
  !-- save the 2D decomposition to disk
  if (mpirank==0) then
    open(14,file='decomposition',status='replace')
    do n=0,ny-1
      do L=0,nz-1
        write(14,'(i4.4,1x)',advance='no') yz_plane_ranks(n,L)
      enddo
      write(14,'(A)',advance='yes') " "
    enddo
    close(14)
  endif
end subroutine fft_initialize


subroutine fft_free
  !====================================================================
  !     Free memory allocated for FFT
  !====================================================================
  use vars
  use p3dfft
  use fftw3_descriptors
  implicit none

  integer :: j

  !-- Clean 3d workspace
  call p3dfft_clean

  !-- Clean 1d workspaces
  do j = 1,3
     call dfftw_destroy_plan(Desc_Handle_1D_f(j))
     call dfftw_destroy_plan(Desc_Handle_1D_b(j))
  enddo

  deallocate(yz_plane_ranks)

end subroutine fft_free


subroutine coftxyz(f,fk)
  !====================================================================
  !     Calculation of the Fourier-coefficients of a real function
  !     along x(1st index),y(2nd index),and z(3rd index)
  !====================================================================
  use vars
  use p3dfft
  use mpi
  implicit none

  real(kind=pr),intent(in) ::  f(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(out) ::  fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: norm

  ! Compute forward FFT
  call p3dfft_ftran_r2c(f,fk,'fff')

  ! Normalize
  norm = 1.d0 / dble(nx*ny*nz)
  fk(:,:,:) = fk(:,:,:) * norm
end subroutine coftxyz


subroutine cofitxyz(fk,f)
  !====================================================================
  !     Calculation of a real function from its Fourier coefficients
  !     along x(1st index),y(2nd index),and z(3rd index)
  !====================================================================
  use vars
  use p3dfft
  use mpi
  implicit none

  complex(kind=pr),intent(in) ::  fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr),intent(out) ::  f(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  ! Compute backward FFT
  call p3dfft_btran_c2r(fk,f,'fff')
end subroutine cofitxyz


!----------------------------------------------------------------
! wavenumber functions: return the kx,ky,kz wavenumbers
! as a function of the array index
!----------------------------------------------------------------

real(kind=pr) function wave_x( ix )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: ix
  wave_x = scalex*dble(ix)
end function

real(kind=pr) function wave_y( iy )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: iy
  wave_y = scaley*dble(modulo(iy+ny/2,ny)-ny/2)
end function

real(kind=pr) function wave_z( iz )
  use vars ! for scale and precision statement
  implicit none
  integer, intent (in) :: iz
  wave_z = scalez*dble(modulo(iz+nz/2,nz)-nz/2)
end function

end module p3dfft_wrapper