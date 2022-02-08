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
  integer(kind=8),dimension(1:3),save :: Desc_Handle_1D_f,Desc_Handle_1D_b

end module fftw3_descriptors

!---------------------------------------------------------------------
! module that contains wrappers for P3DFFT
!---------------------------------------------------------------------

module p3dfft_wrapper
  use vars
  implicit none

  ! this flag tells us whether we are using p3dfft or just initialized the domain
  ! decomposition, in which case we do not use it
  logical, save, private :: using_p3dfft = .true.

  contains


! Compute the FFT of the real-valued 3D array inx and save the output
! in the complex-valued 3D array outk.
subroutine fft(inx, outk)
    use mpi
    use p3dfft
    use vars ! For precision specficiation and array sizes
    real(kind=pr),intent(in)::inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
    real(kind=pr) :: t1
    real(kind=pr) :: norm
    integer(kind=8) :: npoints

    if (using_p3dfft .eqv. .false.) then
      call abort(33343,'P3DFFT is not initialized, you cannot perform FFTs')
    endif

    t1 = MPI_wtime()
    ! Compute forward FFT
    outk = 0.0d0
    call p3dfft_ftran_r2c( inx, outk, 'fff')

    ! Normalize
    !  npoints = int(nx,kind=int64) * int(ny,kind=int64) * int(nz,kind=int64)
    npoints = int(nx,kind=8) * int(ny,kind=8) * int(nz,kind=8)
    norm = 1.d0 / dble(npoints)
    outk = outk * norm

    ! save timing
    call toc( "FFT", MPI_wtime() - t1)
end subroutine fft


! Compute the inverse FFT of the complex-valued 3D array ink and save the
! output in the real-valued 3D array outx.
subroutine ifft(ink, outx)
    use mpi
    use p3dfft
    use vars ! For precision specficiation and array sizes
    complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
    real(kind=pr),intent(out)::outx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    real(kind=pr) :: t1
    t1 = MPI_wtime()

    if (using_p3dfft .eqv. .false.) then
      call abort(33349,'P3DFFT is not initialized, you cannot perform FFTs')
    endif

    ! Compute backward FFT
    call p3dfft_btran_c2r(ink, outx, 'fff')

    ! save timing
    call toc( "iFFT", MPI_wtime() - t1)
end subroutine ifft


! Compute the FFT of the real-valued 3D array inx and save the output
! in the complex-valued 3D array outk.
subroutine fft3(inx, outk)
    use mpi
    use p3dfft
    use vars ! For precision specficiation and array sizes
    real(kind=pr),intent(in)::inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
    complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

    real(kind=pr) :: t1
    real(kind=pr) :: norm
    integer(kind=8) :: npoints

    t1 = MPI_wtime()

    if (using_p3dfft .eqv. .false.) then
      call abort(33343,'P3DFFT is not initialized, you cannot perform FFTs')
    endif

    outk = 0.0d0
    call p3dfft_ftran_r2c(inx(:,:,:,1), outk(:,:,:,1), 'fff')
    call p3dfft_ftran_r2c(inx(:,:,:,2), outk(:,:,:,2), 'fff')
    call p3dfft_ftran_r2c(inx(:,:,:,3), outk(:,:,:,3), 'fff')

    ! Normalize
    npoints = int(nx,kind=8) * int(ny,kind=8) * int(nz,kind=8)
    norm = 1.d0 / dble(npoints)
    outk = outk * norm

    ! save timing
    call toc( "FFT", MPI_wtime() - t1)

end subroutine fft3


! Compute the inverse FFT of the complex-valued 3D array ink and save the
! output in the real-valued 3D array outx.
subroutine ifft3(ink, outx)
    use mpi
    use p3dfft
    use vars ! For precision specficiation and array sizes
    complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
    real(kind=pr),intent(out)::outx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)

    real(kind=pr) :: t1
    real(kind=pr) :: norm
    integer(kind=8) :: npoints

    t1 = MPI_wtime()

    if (using_p3dfft .eqv. .false.) then
      call abort(33343,'P3DFFT is not initialized, you cannot perform iFFTs')
    endif

    call p3dfft_btran_c2r(ink(:,:,:,1), outx(:,:,:,1), 'fff')
    call p3dfft_btran_c2r(ink(:,:,:,2), outx(:,:,:,2), 'fff')
    call p3dfft_btran_c2r(ink(:,:,:,3), outx(:,:,:,3), 'fff')

    ! save timing
    call toc( "iFFT", MPI_wtime() - t1)
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
  integer :: mpicode,idir,L,n,ix,iy,iz,nxc
  integer,dimension(1:3) :: ka,kb,ks,kat,kbt,kst
  logical,dimension(2) :: subcart
  real(kind=pr),dimension(:,:),allocatable :: f,ft
  integer, dimension(:,:), allocatable :: yz_plane_local
  real(kind=pr) :: k

  ! setup a few globals we need throughout the code. the most important information
  ! is the domain size
  ! the scaling factors are used to rescale the FFTs (which are defined on a 2*pi
  ! domain). this is important when computing derivatives.
  scalex = 2.d0*pi/xl
  scaley = 2.d0*pi/yl
  scalez = 2.d0*pi/zl
  ! the grid spacing is always equidistant (since the fourier basis is not nicely
  ! diagonalizing on any stretched grid). note this is also done in params.f90, but
  ! in some cases the params.f90 routines are not called (namely postprocessing)
  ! so just in case, we repeat this here
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  ! yes, we can perform ffts (see also decomposition_initialize)
  using_p3dfft = .true.

  !-----------------------------------------------------------------------------
  !------ Three-dimensional FFT                                           ------
  !-----------------------------------------------------------------------------
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
  endif

  !-- Set up dimensions. It is very important that mpidims(2) > mpidims(1)
  ! because P3Dfft crashes otherwise. This means a 1D decomposition is always
  ! along the z direction in real space.
  if (nz>mpisize) then
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
     decomposition="2D"
  endif

  if (root) write(*,'("mpidims= ",i3,1x,i3)') mpidims
  if (root) write(*,'("Using ",A," decomposition!")') trim(adjustl(decomposition))

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
    call abort(33340,"wrong mpidims: change mpisize")
  endif

  !-- Initialize P3DFFT
  if (iDealias==1) then
    !-- Pruned fft for dealiasing. Only stable in the x direction (Why?)
    nxc = ceiling((3.d0/4.d0)*dble(nx)) ! 3/4*nx is more likely an integer than 3/4*nx
    if (nxc < 2) then ! Do not use pruned for small nx
      nxc = nx
    endif
    call p3dfft_setup(mpidims,nx,ny,nz,MPI_COMM_WORLD,nxc,ny,nz,overwrite=.false.)
    if (mpirank==0) then
      write(*,*) "Using pruned FFT"
    endif
  else
    !-- No pruned if no dealiasing
    call p3dfft_setup(mpidims,nx,ny,nz,MPI_COMM_WORLD,overwrite=.false.)
  endif

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

  if ( rb(2)-ra(2)+1<2*ng .or. rb(3)-ra(3)+1<2*ng ) then
    if (mpirank==0) write(*,*) "Too many CPUs: the ghosts span more than one CPU"
    if (mpirank==0) write(*,*) "y", rb(2)-ra(2)+1, "z", rb(3)-ra(3)+1
    call abort(33341, "Too many CPUs: the ghosts span more than one CPU")
  endif

  !-- Allocate domain partitioning tables and gather sizes from all processes
  !-- (only for real arrays)
  allocate ( ra_table(1:3,0:mpisize-1), rb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER (ra, 3, MPI_INTEGER, ra_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  call MPI_ALLGATHER (rb, 3, MPI_INTEGER, rb_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  !-- cmplx arrays
  allocate ( ca_table(1:3,0:mpisize-1), cb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER (ca, 3, MPI_INTEGER, ca_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  call MPI_ALLGATHER (cb, 3, MPI_INTEGER, cb_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)

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


!   !-----------------------------------------------------------------------------
!   ! ------ Multiple one-dimensional FFTs                                  ------
!   !-----------------------------------------------------------------------------
!   !-- Create Cartesian topology for one-dimensional transforms
!   call MPI_CART_CREATE(MPI_COMM_WORLD,nmpidims,mpidims,(/.false.,.false./),&
!        .false.,mpicommcart,mpicode)
!   ! call MPI_CART_CREATE(MPI_COMM_WORLD,nmpidims,mpidims,(/0,0/),0,mpicommcart,mpicode)
!   call MPI_CART_COORDS(mpicommcart,mpirank,nmpidims,mpicoords,mpicode)
!
!   !-- Restrict communications to slabs:
!   !-- Communicate between processes which have the same second(for idir=1) or third(for idir=2) index in the cartesian topology
!   do idir = 1,2
!      subcart(3-idir) = .true.
!      subcart(idir) = .false.
!      call MPI_CART_SUB(mpicommcart,subcart,mpicommslab(idir),mpicode)
!   enddo
!
!   !-- Create FFTW plans for all possible sizes
!   ! extents of arrays before and after transpose x->y
!   call trextents(1,(/nx,ny,nz/),mpidims,mpicoords,ka,kb,ks,kat,kbt,kst)
!
!   ! allocate plan for transform with x leading
!   L = nx
!   n = ks(2)*ks(3)
!   allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
!   call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(1),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
!        FFTW_ESTIMATE)
!   call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(1),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
!        FFTW_ESTIMATE)
!   deallocate(f,ft )
!
!   ! allocate plan for transform with y leading
!   L = ny
!   n = kst(2)*kst(3)
!   allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
!   call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(2),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
!        FFTW_ESTIMATE)
!   call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(2),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
!        FFTW_ESTIMATE)
!   deallocate(f,ft )
!
!   ! extents of arrays before and after transpose y->z
!   call trextents(2,(/ny,nx,nz/),mpidims,mpicoords,ka,kb,ks,kat,kbt,kst)
!
!   ! allocate plan for transform with z leading
!   L = nz
!   n = kst(2)*kst(3)
!   allocate(f(0:L-1,0:n-1),ft(0:L+1,0:n-1) )
!   call dfftw_plan_many_dft_r2c(Desc_Handle_1D_f(3),1,L,n,f,0,1,L,ft,0,1,L/2+1,&
!        FFTW_ESTIMATE)
!   call dfftw_plan_many_dft_c2r(Desc_Handle_1D_b(3),1,L,n,ft,0,1,L/2+1,f,0,1,L,&
!        FFTW_ESTIMATE)
!   deallocate(f,ft )
!   deallocate (yz_plane_local)
if (mpirank==0) write(*,'(A)') "------------------------------p3dfft init DONE---------------------------"
end subroutine fft_initialize


!-------------------------------------------------------------------------------

subroutine setup_cart_groups
  ! Setup 1d communicators
  use mpi
  use vars
  use p3dfft

  implicit none
  integer :: mpicolor,mpikey,mpicode
  integer :: mpicommtmp1,mpicommtmp2
  logical, dimension(2) :: subcart
  logical :: mpiperiods(2),periods(1),reorder
!  integer :: mpiperiods(2),periods(1),reorder
  integer :: one=1,two=2,dims(1) ! Required for MPI_CART_GET, MPI_CART_CREATE in openmpi

  if (mpirank==0) write(*,*) " setting up cart_groups..."

  ! Set parameters
  periods(1)=.true. ! This should be an array - if not, openmpi fails
  reorder=.false.
!  periods(1)=1
!  reorder=0

  ! Get Cartesian topology information
  call MPI_CART_GET(mpicommcart,two,mpidims,mpiperiods,mpicoords,mpicode)

 ! As the name implies, MPI_Comm_split creates new communicators by "splitting"
 ! a communicator into a group of sub-communicators based on the input values
 ! color and key. The first argument, mpicommcart, is the communicator that will be
 ! used as the basis for the new communicators.
 ! The second argument, color, determines to which new communicator each processes
 ! will belong. All processes which pass in the same value for color are assigned
 ! to the same communicator.
 ! The third argument, key, determines the ordering (rank) within each new communicator.
 ! The process which passes in the smallest value for color will be rank 0, the next
 ! smallest will be rank 1, and so on.
 ! The fourth argument, newcomm is how MPI returns the new communicator back to the user.


  ! Communicator for line in y direction
  ! all ranks that have the same z-value
  mpicolor = mpicoords(2)
  ! and they are sorted by their y-value
  mpikey = mpicoords(1)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp1,mpicode)
  dims(1) = mpidims(1)
  call MPI_CART_CREATE(mpicommtmp1,one,dims,periods,reorder,mpicommz,mpicode)


  ! Communicator for line in z direction
  mpicolor = mpicoords(1)
  mpikey = mpicoords(2)
  call MPI_COMM_SPLIT (mpicommcart,mpicolor,mpikey,mpicommtmp2,mpicode)
  dims(1) = mpidims(2)
  call MPI_CART_CREATE(mpicommtmp2,one,dims,periods,reorder,mpicommy,mpicode)


  mpicommslab = (/mpicommy,mpicommz/)

  if (mpirank==0) write(*,*) " done setting up cart_groups!"
end subroutine setup_cart_groups



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
     write(*,*) "*** cleaning FFT"
  endif

  if (using_p3dfft) then
    call p3dfft_clean
    !-- Clean 1d workspaces
    do j = 1,3
       call dfftw_destroy_plan(Desc_Handle_1D_f(j))
       call dfftw_destroy_plan(Desc_Handle_1D_b(j))
    enddo
  endif

  if (allocated(yz_plane_ranks)) deallocate(yz_plane_ranks)
  if (allocated(rb_table)) deallocate(rb_table)
  if (allocated(ra_table)) deallocate(ra_table)
  if (allocated(cb_table)) deallocate(cb_table)
  if (allocated(ca_table)) deallocate(ca_table)
end subroutine fft_free

!-------------------------------------------------------------------------------

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
  norm = 1.0d0 / real(L)
  do j=0,n-1
     fk(:,j) = ft(0:L-1,j) * norm
  end do

  deallocate(ft )

  !      last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofts

!-------------------------------------------------------------------------------

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
     ft(L:L+1,j) = 0.0d0
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
  use mpi
  implicit none

  integer, intent (in) :: idir, mpireal
  integer, dimension (1:3), intent (in) :: n, ka, kb, ks, kat, kbt, kst
  integer, dimension (2), intent (in) :: mpidims, mpicoords, mpicommslab
  real (kind=pr), dimension (ka(1):kb(1),ka(2):kb(2),ka(3):kb(3)), intent (inout) :: datain
  real (kind=pr), dimension (kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)), intent (inout) :: dataout
  integer :: jj, kk
  integer :: mpicode, mpirealsize, mpiblock
  integer, dimension(:), allocatable :: types, counts
  integer, dimension(:), allocatable :: displs
  real (kind=pr), dimension (:,:,:), allocatable :: datatemp

  !--Check direction
  if ( (idir/=1) .and. (idir/=2) ) then
     call abort(3334333, 'tr2dplus: illegal idir')
  endif

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


!-------------------------------------------------------------------------------
!     Initialize domain decomposition, do not use p3dfft
!-------------------------------------------------------------------------------
! This is very useful in postprocessing, since in many cases we do not actually
! intent to do FFTs, for example when computing the energy of a field, or when
! extrcating subsets, dry runs, etc.
! These routines provide the domain decomposition for REAL data (since complex
! data would imply using FFTs), that is, for all procs, the local bounds
! ra(1:3) and rb(1:3), as well as the communicators used for ghost node synchronizaiton
!-------------------------------------------------------------------------------
subroutine decomposition_initialize(force2D_decomp)
  use mpi ! Module incapsulates mpif.
  use vars
  implicit none

  integer,parameter :: nmpidims = 2
  integer :: mpicode,idir,L,n
  integer,dimension(1:3) :: ka,kb,ks,kat,kbt,kst
  logical,dimension(2) :: subcart
  logical, intent(in), optional :: force2D_decomp
  real(kind=pr),dimension(:,:),allocatable :: f,ft

  ! setup a few globals we need throughout the code. the most important information
  ! is the domain size
  ! the scaling factors are used to rescale the FFTs (which are defined on a 2*pi
  ! domain). this is important when computing derivatives.
  scalex = 2.d0*pi/xl
  scaley = 2.d0*pi/yl
  scalez = 2.d0*pi/zl
  ! the grid spacing is always equidistant (since the fourier basis is not nicely
  ! diagonalizing on any stretched grid). note this is also done in params.f90, but
  ! in some cases the params.f90 routines are not called (namely postprocessing)
  ! so just in case, we repeat this here
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  ! default case, no decomposition
  decomposition="none"
  ! to avoid user errors, we set the flag to false and yell if we try to perform
  ! ffts anyways:
  using_p3dfft = .false.

  if (mpirank==0) then
    write(*,'(A)') "-----------------------------------decomposition init (NO P3DFFT)---------------------------"
    write(*,'("Initializing decomposition, n=(",3(i4,1x),")")') nx,ny,nz
    write(*,'("DECOMPOSITION domain size:     ",3(g15.8,1x))') xl,yl,zl
    write(*,'("DECOMPOSITION scale factors:   ",3(g15.8,1x))') scalex, scaley, scalez
    write(*,'("DECOMPOSITION lattice spacing: ",3(g15.8,1x))') dx,dy,dz
  endif

  !-- Set up dimensions. It is very important that mpidims(2) > mpidims(1)
  ! because P3Dfft crashes otherwise. This means a 1D decomposition is always
  ! along the z direction in real space.
  if (present(force2D_decomp)) then
      if ((nz>mpisize.or.nx==1).and.(force2D_decomp.eqv..false.)) then
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
  else
      if (nz>mpisize.or.nx==1) then
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
  endif

  if (root) write(*,'("mpidims= ",i3,1x,i3)') mpidims
  if (root) write(*,'("Using ",A," decomposition!")') trim(adjustl(decomposition))

  !-- Check dimensions
  if(mpidims(1)*mpidims(2)/=mpisize) then
     call abort(3333334,'wrong mpidims: change mpisize')
  endif

  !-- Set subdomain bounds
  !-- Get Cartesian topology info
  !-- Get local sizes
  call p3dfft_stub(mpidims,nx,ny,nz,MPI_COMM_WORLD,mpitaskid,mpitasks,mpicommcart,ra,rb,rs)
  ra(:) = ra(:) - 1
  rb(:) = rb(:) - 1

  !-- extents of real arrays that have ghost points. We add ghosts in all
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
    call abort(3334509, "Too many CPUs: the ghosts span more than one CPU")
  endif

  !-- Allocate domain partitioning tables and gather sizes from all processes
  !-- (only for real arrays)
  ! TODO: These tables are currently not used for communication between subdomains,
  ! but may be still useful for development/debugging purposes.
  allocate ( ra_table(1:3,0:mpisize-1), rb_table(1:3,0:mpisize-1) )
  call MPI_ALLGATHER (ra, 3, MPI_INTEGER, ra_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)
  call MPI_ALLGATHER (rb, 3, MPI_INTEGER, rb_table, 3, MPI_INTEGER, MPI_COMM_WORLD, mpicode)

  if (mpirank==0) write(*,'(A)') "------------------------------decomposition init DONE---------------------------"
end subroutine decomposition_initialize


!----------------------------------------------------------------
! domain decomposition routine derived from P3DFFT
!----------------------------------------------------------------
subroutine p3dfft_stub(dims_in,nx,ny,nz,mpi_comm_in,mpi_taskid,mpi_tasks,mpi_comm_out,istart,iend,isize)
  implicit none

  integer, intent (in) :: nx,ny,nz,mpi_comm_in
  integer, intent (out) :: istart(3),iend(3),isize(3)
  integer, intent (out) :: mpi_taskid,mpi_tasks,mpi_comm_out
  integer, intent (in) :: dims_in(2)

  integer :: i,j,k,ierr,mpicomm,mpi_comm_cart
  integer :: ipid,jpid,iproc,jproc
  integer :: numtasks,taskid
  integer :: cartid(2),dims(2)
  integer, dimension (:), allocatable :: jist,jisz,jien,kjst,kjsz,kjen
  logical :: periodic(2),remain_dims(2)

  if(nx .le. 0 .or. ny .le. 0 .or. nz .le. 0) then
     print *,'Invalid dimensions :',nx,ny,nz
     call abort(333952, 'p3dfft_stub()::Invalid dimensions')
  endif

  mpicomm = mpi_comm_in
  call MPI_COMM_SIZE (mpicomm,numtasks,ierr)
  call MPI_COMM_RANK (mpicomm,taskid,ierr)

  if(dims_in(1) .le. 0 .or. dims_in(2) .le. 0 .or.  dims_in(1)*dims_in(2) .ne. numtasks) then
     print *,'Invalid processor geometry: ',dims,' for ',numtasks, 'tasks'
     call abort(333456, 'p3dfft_stub()::Invalid processor geometry')
  endif

  if(taskid .eq. 0) then
     print *,'Using stride-1 layout'
  endif

  iproc = dims_in(1)
  jproc = dims_in(2)
  dims(1) = dims_in(2)
  dims(2) = dims_in(1)

  periodic(1) = .false.
  periodic(2) = .false.
! creating cartesian processor grid
  call MPI_Cart_create(mpicomm,2,dims,periodic,.false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
  call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
! process with a linear id of 5 may have cartid of (3,1)

  ipid = cartid(2)
  jpid = cartid(1)

  allocate (jist(0:iproc-1))
  allocate (jisz(0:iproc-1))
  allocate (jien(0:iproc-1))
  allocate (kjst(0:jproc-1))
  allocate (kjsz(0:jproc-1))
  allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny,nz) => (iproc,jproc)
!
  call MapDataToProc(ny,iproc,jist,jien,jisz)
  call MapDataToProc(nz,jproc,kjst,kjen,kjsz)

! These are local array indices for each processor
  istart(1) = 1
  iend(1) = nx
  isize(1) = nx
  istart(2) = jist(ipid)
  iend(2) = jien(ipid)
  isize(2) = jisz(ipid)
  istart(3) = kjst(jpid)
  iend(3) = kjen(jpid)
  isize(3) = kjsz(jpid)

  deallocate(jist,jisz,jien,kjst,kjsz,kjen)

  mpi_taskid = taskid
  mpi_tasks = numtasks
  mpi_comm_out = mpi_comm_cart

end subroutine p3dfft_stub


!----------------------------------------------------------------
! calculate subdomain bounds
!----------------------------------------------------------------
subroutine MapDataToProc (data,proc,st,en,sz)
  implicit none
  integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
  integer i,size,nl,nu

  size=data/proc
  nu = data - size * proc
  nl = proc - nu
  st(0) = 1
  sz(0) = size
  en(0) = size
  do i=1,nl-1
     st(i) = st(i-1) + size
     sz(i) = size
     en(i) = en(i-1) + size
  enddo
  size = size + 1
  do i=nl,proc-1
     st(i) = en(i-1) + 1
     sz(i) = size
     en(i) = en(i-1) + size
  enddo
  en(proc-1)= data
  sz(proc-1)= data-st(proc-1)+1
end subroutine

end module p3dfft_wrapper
