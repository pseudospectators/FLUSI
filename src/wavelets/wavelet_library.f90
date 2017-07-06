module flusi_wavelet_lib

  use vars, only : pr, ra, rb, nx, ny, nz, on_proc, mpidims, mpicoords, mpicommslab, abort
  use mpi
  use p3dfft_wrapper

  implicit none

  ! largest number of filter coefficients allowed in program:
  integer, parameter :: MaxFlt = 30
  !
  type orth_wavelet
    real(kind=pr), dimension(:), allocatable :: HD ! low pass decomposition filter
    real(kind=pr), dimension(:), allocatable :: GD ! high pass decomposition filter
    real(kind=pr), dimension(:), allocatable :: HR ! low pass reconstruction filter
    real(kind=pr), dimension(:), allocatable :: GR ! high pass reconstruction filter
    character(len=80) :: name
  end type

contains


  ! 1D fast wavelet transform, periodized, orthogonal
  subroutine FWT1_PO(u, wc, wavelet)
    implicit none
    ! derived datatype containing the filters
    type(orth_wavelet),intent(in) :: wavelet
    ! signal to be transformed
    real(kind=pr), dimension(1:), intent(in) :: u
    ! resulting wavelet coefficients
    real(kind=pr), dimension(1:), intent(inout) :: wc
    ! signal length
    integer :: N, J0, nc, Jmin, j
    ! buffer
    real(kind=pr), allocatable, dimension(:) :: buffer1, buffer2

    ! get length and level of input signal
    call dyadlength(u, N, J0)

    nc = N
    Jmin = 0
    wc = u

    ! FWT goes down in levels starting from the signal sampling
    ! (finest level)
    do j = J0-1, Jmin, -1
        allocate( buffer1(1:nc), buffer2(1:nc) )

        ! low-pass filter
        call filter1P(wc(1:nc), buffer1, wavelet%HD, lbound(wavelet%HD,dim=1), ubound(wavelet%HD,dim=1))
        ! high-pass filter
        call filter1P(wc(1:nc), buffer2, wavelet%GD, lbound(wavelet%GD,dim=1), ubound(wavelet%GD,dim=1))

        ! decimation by 2
        wc(1:nc/2) = buffer1(1:nc-1:2)
        wc(nc/2+1:nc) = buffer2(1:nc-1:2)

        nc = nc / 2
        deallocate(buffer1, buffer2)
    enddo

    call init_empty_file('yo.txt')
    open(14,file='yo.txt',status='unknown',position='append')
    do j = 1, N
      write(14,'(es15.8)') wc(j)
    enddo
    close(14)

  end subroutine FWT1_PO


  !-----------------------------------------------------------------------------
  ! 3D fast wavelet transform, periodized, orthogonal
  !-----------------------------------------------------------------------------
  subroutine FWT3_PO(u, wc, wavelet)
    implicit none
    ! derived datatype containing the filters
    type(orth_wavelet),intent(in) :: wavelet
    ! signal to be transformed
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: u
    ! resulting wavelet coefficients
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
    ! signal length
    integer :: N, J0, nc, Jmin, j, ix,iy,iz, idir, mpicode
    ! buffer
    real(kind=pr), allocatable, dimension(:) :: buffer1, buffer2
    real(kind=pr), allocatable, dimension(:,:,:) :: work
    integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst

    N = nx ! this works currently only for nx = ny = nz
    J0 = ceiling( log(dble(N)) / log(2.d0) )

    ! the first iteration acts on the complete data
    nc = N
    Jmin = 1
    ! initialize wavelet coefficients
    wc(:,:,:) = u(:,:,:)

    ! -------------------------------------------------------
    ! FWT goes down in levels starting from the finest level
    ! -------------------------------------------------------
    do j = J0-1, Jmin, -1
        ! an inplace filtering is not possible, so we allocate buffers for a line
        ! (of current size)
        allocate( buffer1(1:nc), buffer2(1:nc) )

        ! ---------------------
        ! --- x direction -----
        ! ---------------------
        ! Filters are applied in the first coordinate, as it is always contiguous in memory
        ! and not split among MPI processes
        call WaveDecomposition_dim1( wc, nc, wavelet, buffer1, buffer2 )


        ! ---------------------
        ! --- y direction -----
        ! ---------------------
        if (mpidims(2) == 1) then
          ! CASE A): index 2 is contiguous and not split among procs.
          call WaveDecomposition_dim2( wc, nc, wavelet, buffer1, buffer2 )

        else
          ! transposition: exchange x-y data
          ! idir = 1 - permute 1 and 2 indices
          ! idir = 2 - permute 1 and 3 indices
          idir = 1
          call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
          allocate( work(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
          work = 0.d0
          call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                       MPI_DOUBLE_PRECISION, wc, work )


          ! at this point we have work(iy,ix,iz)
          ! again, the first dimension is filtered, but this now is the y-data
          call WaveDecomposition_dim1( work, nc, wavelet, buffer1, buffer2 )

          ! transposition: exchange x-y data
          ! idir = 1 - permute 1 and 2 indices
          ! idir = 2 - permute 1 and 3 indices
          idir = 1
          call trextents(idir, (/ny,nx,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
          call subtr ( idir, (/ny,nx,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                       MPI_DOUBLE_PRECISION, work, wc )
          ! at this point we have again wc(ix,iy,iz), which is what we started with
        endif


        idir = 2
        call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
        if(allocated(work)) deallocate(work)
        allocate( work(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
        call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                     MPI_DOUBLE_PRECISION, wc, work )

        ! --- z direction -----
        ! at this point we have work(iz,iy,ix)
        ! again, the first dimension is filtered, but this now is the z-data
        call WaveDecomposition_dim1( work, nc, wavelet, buffer1, buffer2 )

        ! final permutation back to the original dataset
        idir = 2
        call trextents(idir, (/nz,ny,nx/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
        call subtr ( idir, (/nz,ny,nx/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                     MPI_DOUBLE_PRECISION, work, wc )
        if(allocated(work)) deallocate(work)
        ! at this point we have wc(ix,iy,iz)

        nc = nc / 2
        write(*,*) nc
        deallocate(buffer1, buffer2)
    enddo


    call save_field_hdf5(0.d0,'wc_00.h5',wc)

  end subroutine FWT3_PO


  !-----------------------------------------------------------------------------
  ! Along the first dimension of array u, apply low- and high-pass filters, decimate
  ! their result and sort them into the first nc elements of the array.
  ! The code acts only on the first nc datapoints (which is nc=nx in the beginning)
  ! in the second iteration, nc=nx/2 and so forth.
  !
  ! Input:  uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
  ! Output: hhhhhhhhhhhhhhhhgggggggggggggggg
  !         | low pass     ||  high pass   |
  !-----------------------------------------------------------------------------
  subroutine WaveDecomposition_dim1( wc, nc, wavelet, buffer1, buffer2 )
    implicit none
    integer, intent(in) :: nc
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
    real(kind=pr), dimension(1:nc), intent(inout) :: buffer1, buffer2
    type(orth_wavelet), intent(in) :: wavelet
    integer :: ix,iy,iz

    ! Filters are applied in the first coordinate, as it is always contiguous in memory
    ! and not split among MPI processes
    do iy = 0, nc-1 ! note loop bounds 0,nc-1 (and not 0,ny-1)
      do iz = 0, nc-1
        ! is this on my local bounds? (MPI-sense, some procs may idle)
        if (on_proc((/0,iy,iz/))) then
          ! low-pass filter (scaling function)
          call filter1P(wc(0:nc-1,iy,iz), buffer1, wavelet%HD, lbound(wavelet%HD,dim=1), ubound(wavelet%HD,dim=1))
          ! high-pass filter (these guys are the details)
          call filter1P(wc(0:nc-1,iy,iz), buffer2, wavelet%GD, lbound(wavelet%GD,dim=1), ubound(wavelet%GD,dim=1))

          ! decimation by 2, sort into array
          wc(0:nc/2-1,iy,iz) = buffer1(1:nc-1:2)
          wc(nc/2:nc-1,iy,iz) = buffer2(1:nc-1:2)
        endif
      enddo
    enddo
  end subroutine WaveDecomposition_dim1


  ! Along the second  dimension of array u, apply low- and high-pass filters, decimate
  ! their result and sort them into the first nc elements of the array.
  ! NOTE: This procedure MUST NOT be called if the second index is split among MPI procs
  subroutine WaveDecomposition_dim2( wc, nc, wavelet, buffer1, buffer2 )
    implicit none
    integer, intent(in) :: nc
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
    real(kind=pr), dimension(1:nc), intent(inout) :: buffer1, buffer2
    type(orth_wavelet), intent(in) :: wavelet
    integer :: ix,iy,iz

    if (mpidims(2)>1) then
      call abort(99919991, 'WaveDecomposition_dim2')
    endif

    ! NOTE: we assume dimension 2 to be contiguous!!
    do ix = 0, nc-1
      do iz = 0, nc-1
        ! is this on my local bounds? (MPI-sense, some procs may idle)
        if (on_proc((/ix,0,iz/))) then
          ! low-pass filter
          call filter1P(wc(ix,0:nc-1,iz), buffer1, wavelet%HD, lbound(wavelet%HD,dim=1), ubound(wavelet%HD,dim=1))
          ! high-pass filter
          call filter1P(wc(ix,0:nc-1,iz), buffer2, wavelet%GD, lbound(wavelet%GD,dim=1), ubound(wavelet%GD,dim=1))

          ! decimation by 2
          wc(ix,0:nc/2-1,iz) = buffer1(1:nc-1:2)
          wc(ix,nc/2:nc-1,iz) = buffer2(1:nc-1:2)
        endif
      enddo
    enddo
  end subroutine WaveDecomposition_dim2

  ! sets up the coiflet coefficients for the given order "order"
  subroutine setup_coiflet_coefs( order, wavelet )
    implicit none
    ! coiflet parameter
    integer, intent(in) :: order
    ! derived datatype containing the filters
    type(orth_wavelet),intent(inout) :: wavelet

    if (allocated(wavelet%HD)) then
      write(*,*) 'wavelet seems to be already initialized...error.'
      stop
    endif

    ! source: http://wavelets.pybytes.com/wavelet/coif2/
    ! and fields-code // FRA_CVE code from japan.
    wavelet%name = "coiflet12"
    allocate( wavelet%HD(-4:7) )
    allocate( wavelet%GD(-6:5) )
    allocate( wavelet%HR(-4:7) )
    allocate( wavelet%GR(-6:5) )

    ! low pass decomposition filter
    wavelet%HD = (/1.638733646318000d-02, -4.146493678197000d-02,-6.737255472230000d-02,&
    3.861100668230900d-01,8.127236354496100d-01,4.170051844237800d-01,-7.648859907826000d-02,&
    -5.943441864647000d-02,2.368017194688000d-02,5.611434819370000d-03,-1.823208870910000d-03,&
    -7.205494453700000d-04/)

    ! Low-pass reconstruction filter (same as decomposition above)
    wavelet%HR = (/1.638733646318000d-02, -4.146493678197000d-02,-6.737255472230000d-02,&
    3.861100668230900d-01,8.127236354496100d-01,4.170051844237800d-01,-7.648859907826000d-02,&
    -5.943441864647000d-02,2.368017194688000d-02,5.611434819370000d-03,-1.823208870910000d-03,&
    -7.205494453700000d-04/)

    ! high-pass filter. In wavelab, this is MirrorFilt(qmf) but note in wavlélab
    ! the coefficients are reversed. (why?)
    ! this gives detail-coefficients
    wavelet%GD = (/ 7.205494453700000d-04,-1.823208870910000d-03,-5.611434819370000d-03,&
    2.368017194688000d-02,5.943441864647000d-02,-7.648859907826000d-02,-4.170051844237800d-01,&
    8.127236354496100d-01,-3.861100668230900d-01,-6.737255472230000d-02,4.146493678197000d-02,&
    1.638733646318000d-02/)

    ! high-pass reconstruction filter
    wavelet%GR = (/7.205494453700000d-04,-1.823208870910000d-03,-5.611434819370000d-03,&
    2.368017194688000d-02,5.943441864647000d-02,-7.648859907826000d-02,-4.170051844237800d-01,&
    8.127236354496100d-01,-3.861100668230900d-01,-6.737255472230000d-02,4.146493678197000d-02,&
    1.638733646318000d-02/)

    write(*,*) "HD",wavelet%HD
    write(*,*) "GD",wavelet%GD
  end subroutine



  ! one-dimensional, periodized filter
  subroutine filter1P(u, u_filtered, filter, i1, i2)
    implicit none
    ! the actual filter
    real(kind=pr), dimension(i1:i2), intent(in) :: filter
    ! filter bound indices
    integer, intent(in) :: i1, i2
    ! the signal to be filtered
    real(kind=pr), dimension(1:), intent(in) :: u
    ! the resulting filtered signal
    real(kind=pr), dimension(1:), intent(inout) :: u_filtered

    integer :: N, i, j

    N = size(u)
    u_filtered = 0.d0

    ! apply filter f to periodic signal u, i.e. convolute with filter
    do i = 1, N
        do j = i1, i2
            u_filtered(i) = u_filtered(i) + u( perindex(i+j,N) ) * filter(j)
        enddo
    enddo

  end subroutine filter1P



  ! periodic index
  ! uses one-based indexing
  integer function perindex(i,n)
    implicit none
    integer, intent(in) :: i, n
    perindex = i
    do while (perindex<1 .or. perindex>n)
        if (perindex<1) then
            perindex = perindex+n
        elseif (perindex>n) then
            perindex = perindex-n
        endif
    end do

  end function



  subroutine dyadlength(u, N, J)
    implicit none
    ! the signal
    real(kind=pr), intent(in), dimension(:) :: u
    !
    integer, intent(out) :: N
    !
    integer, intent(out) :: J
    ! dyadlength -- Find length and dyadic length of array
    !  Usage
    !    [n,J] = dyadlength(x)
    !  Inputs
    !    x    array of length n = 2^J (hopefully)
    !  Outputs
    !    n    length(x)
    !    J    least power of two greater than n
    !
    !  Side Effects
    !    A warning is issued if n is not a power of 2.
    !
    !  See Also
    !    quadlength, dyad, dyad2ix
    !
    N = size(u,1)
    J = ceiling( log(dble(N)) / log(2.d0) )
  end subroutine
end module
