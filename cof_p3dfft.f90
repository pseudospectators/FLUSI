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


subroutine fft_initialize
  !====================================================================
  !     Allocate memory and initialize FFT
  !====================================================================
  use mpi_header ! Module incapsulates mpif.
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

  ! ------ Three-dimensional FFT ------
  !-- Set up dimensions
  !-- !!! Note that p3dfft_setup permutes mpidims, if DIMS_C is not defined !!!
  if(ny>mpisize) then
     mpidims(1) = 1
     mpidims(2) = mpisize / mpidims(1)
  else
     mpidims(1) = 0
     mpidims(2) = 0
     call MPI_Dims_create(mpisize,nmpidims,mpidims,mpicode)
     if(mpidims(1) > mpidims(2)) then
        mpidims(1) = mpidims(2)
        mpidims(2) = mpisize / mpidims(1)
     endif
  endif

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
     print *, 'wrong mpidims: change mpisize'
     stop
  endif

  !-- Stop if the sizes are odd or smaller than 4
  if(nx<4 .or. ny<4 .or. nz<4 .or. modulo(nx,2)==1 .or. modulo(ny,2)==1 .or. modulo(nz,2)==1 ) then
     if(mpirank == 0 ) then
        print *, 'nx, ny, nz must be even and not smaller than 4'
     endif
     stop
  endif

  !-- Initialize P3FFT
  call p3dfft_setup(mpidims,nx,ny,nz,.false.)

  !-- Get local sizes
  call get_dims(ra,rb,rs,1)  ! real blocks
  call get_dims(ca,cb,cs,2)  ! complex blocks
  ra(:) = ra(:) - 1
  rb(:) = rb(:) - 1
  ca(:) = ca(:) - 1
  cb(:) = cb(:) - 1

  !-- Allocate domain partitioning tables and gather sizes from all
  !-- processes (only for real arrays)
  allocate(ra_table(1:3,0:mpisize-1),rb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER(ra,3,mpiinteger,ra_table,3,mpiinteger,MPI_COMM_WORLD,&
       mpicode)
  call MPI_ALLGATHER(rb,3,mpiinteger,rb_table,3,mpiinteger,MPI_COMM_WORLD,&
       mpicode)

  ! ------ Multiple one-dimensional FFTs ------
  !-- Create Cartesian topology for one-dimensional transforms
  call MPI_CART_CREATE(MPI_COMM_WORLD,nmpidims,mpidims,(/.false.,.false./),&
       .false.,mpicommcart,mpicode)
  ! call MPI_CART_CREATE(MPI_COMM_WORLD,nmpidims,mpidims,(/0,0/),0,mpicommcart,mpicode)
  call MPI_CART_COORDS(mpicommcart,mpirank,nmpidims,mpicoords,mpicode)

  !-- Restrict communications to slabs:
  !-- Communicate between processes which have the same second(for idir=1) or third(for idir=2) index in the cartesian topology
  do idir = 1,2
     subcart(3-idir) = .true.
     subcart(idir) = .false.
     call MPI_CART_SUB(mpicommcart,subcart,mpicommslab(idir),mpicode)
  enddo

  !-- Create FFTW plans for all possible sizes
  ! extents of arrays before and after transpose x->y
  call trextents(1,(/nx,ny,nz/),mpidims,mpicoords,ka,kb,ks,kat,kbt,kst)

  ! allocate plan for transform with x leading
  L = nx
  n = ks(2)*ks(3)
  allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
  call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(1),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
       FFTW_ESTIMATE)
  call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(1),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
       FFTW_ESTIMATE)
  deallocate(f,ft )

  ! allocate plan for transform with y leading
  L = ny
  n = kst(2)*kst(3)
  allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
  call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(2),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
       FFTW_ESTIMATE)
  call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(2),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
       FFTW_ESTIMATE)
  deallocate(f,ft )

  ! extents of arrays before and after transpose y->z
  call trextents(2,(/ny,nx,nz/),mpidims,mpicoords,ka,kb,ks,kat,kbt,kst)

  ! allocate plan for transform with z leading
  L = nz
  n = kst(2)*kst(3)
  allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
  call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(3),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
       FFTW_ESTIMATE)
  call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(3),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
       FFTW_ESTIMATE)
  deallocate(f,ft )

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
  if(mpirank ==0) then
     write(*,'(A)') "!!! cof_p3dfft.f90: calling p3dfft_clean"
  endif

  call p3dfft_clean

  !-- Clean 1d workspaces
  do j = 1,3
     if(mpirank ==0) then
        write(*,'("!!! cof_p3dfft.f90: cleaning 1d workspace j=",i1)') j
     endif
     call dfftw_destroy_plan(Desc_Handle_1D_f(j))
     call dfftw_destroy_plan(Desc_Handle_1D_b(j))
  enddo

  !-- Deallocate domain partitioning tables
  deallocate(ra_table,rb_table )

end subroutine fft_free


subroutine coftxyz(f,fk)
  !====================================================================
  !     Calculation of the Fourier-coefficients of a real function
  !     along x(1st index),y(2nd index),and z(3rd index)
  !====================================================================
  use vars
  use p3dfft
  use mpi_header
  implicit none

  real(kind=pr),intent(in) ::  f(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(out) ::  fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr),save :: t1

  t1 = MPI_wtime()
  ! Compute forward FFT
  call p3dfft_ftran_r2c(f,fk)

  ! Normalize
  fk(:,:,:) = fk(:,:,:) / dble(nx*ny*nz)

  time_fft  = time_fft  + MPI_wtime() - t1  ! for global % of FFTS
  time_fft2 = time_fft2 + MPI_wtime() - t1  ! for % FFT in cal_nlk only

  ! Filter Nyquist frequency
  !if((ca(2)<=nx/2) .and.(cb(2)>=nx/2) ) then
  !    fk(:,nx/2,:) = 0.0d0
  !endif
  !if((ca(3)<=ny/2) .and.(cb(3)>=ny/2) ) then
  !    fk(:,:,ny/2) = 0.0d0
  !endif
  !if((ca(1)<=nz/2) .and.(cb(1)>=nz/2) ) then
  !    fk(nz/2,:,:) = 0.0d0
  !endif

end subroutine coftxyz


subroutine cofitxyz(fk,f)
  !====================================================================
  !     Calculation of a real function from its Fourier coefficients
  !     along x(1st index),y(2nd index),and z(3rd index)
  !====================================================================
  use vars
  use p3dfft
  use mpi_header
  implicit none

  complex(kind=pr),intent(in) ::  fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr),intent(out) ::  f(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),save :: t1
  t1 = MPI_wtime()

  ! Compute backward FFT
  call p3dfft_btran_c2r(fk,f)


  time_ifft  = time_ifft  + MPI_wtime() - t1
  time_ifft2 = time_ifft2 + MPI_wtime() - t1
end subroutine cofitxyz


subroutine cofts(iplan,f,fk,L,n)
  !====================================================================
  !     Calculation of the Fourier coefficients of n real functions
  !     aligned in columns
  !     FILTERING
  !     FFTW3 library is used
  !====================================================================
  use vars
  use fftw3_descriptors
  implicit none
  include 'fftw3.f'

  integer,intent(in) :: iplan,L,n
  integer :: j
  real(kind=pr),dimension(:,:),allocatable :: ft
  real(kind=pr),dimension(0:L-1,0:n-1),intent(in) ::  f
  real(kind=pr),dimension(0:L-1,0:n-1),intent(out) ::  fk
  real(kind=pr) :: norm

  allocate(ft(0:L+1,0:n-1) )

  call dfftw_execute_dft_r2c(Desc_Handle_1D_f(iplan),f,ft)

  ! Output
  norm = 1.0 / real(L)
  do j=0,n-1
     fk(:,j) = ft(0:L-1,j) * norm
  end do

  deallocate(ft )

  !      last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofts


subroutine cofits(iplan,fk,f,L,n)
  !====================================================================
  !     Calculation of n real functions from their Fourier coefficients
  !     aligned in columns
  !     FILTERING
  !     FFTW3 library is used
  !====================================================================
  use vars
  use fftw3_descriptors
  implicit none
  include 'fftw3.f'

  integer,intent(in) :: iplan,L,n
  integer :: j
  real(kind=pr),dimension(:,:),allocatable :: ft
  real(kind=pr),dimension(0:L-1,0:n-1),intent(in) ::  fk
  real(kind=pr),dimension(0:L-1,0:n-1),intent(out) ::  f

  allocate(ft(0:L+1,0:n-1) )

  ! Input
  do j=0,n-1
     ft(0:L-1,j) = fk(:,j)
     ft(L:L+1,j) = 0.0
  end do

  call dfftw_execute_dft_c2r(Desc_Handle_1D_b(iplan),ft,f)

  deallocate(ft )

  ! last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofits

! Compute extents of transposed arrays
subroutine trextents(idir,n,mpidims,mpicoords,ka,kb,ks,kat,kbt,kst )
  implicit none

  integer,intent(in) :: idir
  integer,dimension(1:3),intent(in) :: n
  integer,dimension(2),intent(in) :: mpidims,mpicoords
  integer,dimension(1:3),intent(out) :: ka,kb,ks,kat,kbt,kst

  ka(1) = 0
  kb(1) = n(1)-1
  ka(2) = mpicoords(2)*n(2)/mpidims(2)
  kb(2) =(mpicoords(2)+1)*n(2)/mpidims(2)-1
  ka(3) = mpicoords(1)*n(3)/mpidims(1)
  kb(3) =(mpicoords(1)+1)*n(3)/mpidims(1)-1
  ks(:) = kb(:)-ka(:)+1

  kat(1) = 0
  kbt(1) = n(idir+1)-1
  kat(1+idir) = mpicoords(3-idir)*n(1)/mpidims(3-idir)
  kbt(1+idir) =(mpicoords(3-idir)+1)*n(1)/mpidims(3-idir)-1
  kat(4-idir) = ka(4-idir)
  kbt(4-idir) = kb(4-idir)
  kst(:) = kbt(:)-kat(:)+1

end subroutine trextents
