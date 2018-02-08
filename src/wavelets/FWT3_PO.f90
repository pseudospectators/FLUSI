!-----------------------------------------------------------------------------
! 3D fast wavelet transform, periodized, orthogonal
!-----------------------------------------------------------------------------
subroutine FWT3_PO(u, wc, wavelet, Jmin, nc)
  implicit none
  ! derived datatype containing the filters
  type(orth_wavelet),intent(in) :: wavelet
  ! signal to be transformed
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: u
  ! resulting wavelet coefficients
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
  ! coarse level and number of points on coarse level
  integer, intent(inout) :: Jmin, nc
  ! signal length
  integer :: N, J0, j, ix,iy,iz, idir, mpicode
  ! buffer
  real(kind=pr), allocatable, dimension(:) :: buffer1, buffer2
  real(kind=pr), allocatable, dimension(:,:,:) :: work
  integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst

  N = nx ! this works currently only for nx = ny = nz
  if (nx /= ny .or. ny /= nz .or. nx /= nz) then
    call abort(55526,"Wavelet error: condition nx=ny=nz is violated.")
  endif
  ! finest level
  J0 = ceiling( log(dble(N)) / log(2.d0) )

  ! the first iteration acts on the complete data
  nc = N
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
      call WaveDecomposition_dim1( wc, ra, rb, nc, wavelet, buffer1, buffer2 )


      ! ---------------------
      ! --- y direction -----
      ! ---------------------
      if (mpidims(2) == 1) then
        ! CASE A): index 2 is contiguous and not split among procs.
        call WaveDecomposition_dim2( wc, ra, rb, nc, wavelet, buffer1, buffer2 )

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
        call WaveDecomposition_dim1( work, kat, kbt, nc, wavelet, buffer1, buffer2 )

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
      call WaveDecomposition_dim1( work, kat, kbt, nc, wavelet, buffer1, buffer2 )

      ! final permutation back to the original dataset
      idir = 2
      call trextents(idir, (/nz,ny,nx/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
      call subtr ( idir, (/nz,ny,nx/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                   MPI_DOUBLE_PRECISION, work, wc )
      if(allocated(work)) deallocate(work)
      ! at this point we have wc(ix,iy,iz)

      nc = nc / 2

      ! if the number of points becomes odd, we've reached the end of possible
      ! divisions by 2, so we stop here
      if ( modulo(nc,2) /= 0 ) then

        ! return minimum level (note this connotation of level is a bit strange
        ! if you do not have a power-2 grid)
        Jmin = j
        ! the inverse transform needs to know how many points it starts with
        ! which is the last division that gives an even number
        nc = nc * 2
        ! we still have to deallocate the buffers
        deallocate(buffer1, buffer2)
        ! break loop
        exit
      endif

      deallocate(buffer1, buffer2)
  enddo

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
  subroutine WaveDecomposition_dim1( wc, qa, qb, nc, wavelet, buffer1, buffer2 )
    implicit none
    integer, intent(in) :: nc
    ! note this routine is one of the few where the decomposition has to be general
    integer, dimension(1:3), intent(in) :: qa, qb
    real(kind=pr), dimension(qa(1):qb(1),qa(2):qb(2),qa(3):qb(3)), intent(inout) :: wc
    real(kind=pr), dimension(1:nc), intent(inout) :: buffer1, buffer2
    type(orth_wavelet), intent(in) :: wavelet
    integer :: ix,iy,iz

    ! Filters are applied in the first coordinate, as it is always contiguous in memory
    ! and not split among MPI processes
    do iy = 0, nc-1 ! note loop bounds 0,nc-1 (and not 0,ny-1)
      do iz = 0, nc-1
        ! is this on my local bounds? (MPI-sense, some procs may idle)
        if (iy>=qa(2) .and. iy<=qb(2) .and. iz>=qa(3) .and. iz<=qb(3)) then
          ! low-pass filter (scaling function)
          call filter1P(wc(0:nc-1,iy,iz), buffer1, wavelet%HD, lbound(wavelet%HD,dim=1), ubound(wavelet%HD,dim=1))
          ! high-pass filter (these guys are the details)
          call filter1P(wc(0:nc-1,iy,iz), buffer2, wavelet%GD, lbound(wavelet%GD,dim=1), ubound(wavelet%GD,dim=1))

          ! decimation by 2, sort into array
          wc(0:nc/2-1,iy,iz) = buffer1(1:nc-1:2)  ! bot: 0:nc/2-1 (zero-based)
          wc(nc/2:nc-1,iy,iz) = buffer2(1:nc-1:2) ! top: nc/2:nc-1 (zero-based)
        endif
      enddo
    enddo
  end subroutine WaveDecomposition_dim1


  ! Along the second  dimension of array u, apply low- and high-pass filters, decimate
  ! their result and sort them into the first nc elements of the array.
  ! NOTE: This procedure MUST NOT be called if the second index is split among MPI procs
  subroutine WaveDecomposition_dim2( wc, qa, qb, nc, wavelet, buffer1, buffer2 )
    implicit none
    integer, intent(in) :: nc
    ! note this routine is one of the few where the decomposition has to be general
    integer, dimension(1:3), intent(in) :: qa, qb
    real(kind=pr), dimension(qa(1):qb(1),qa(2):qb(2),qa(3):qb(3)), intent(inout) :: wc
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
        if (ix>=qa(1) .and. ix<=qb(1) .and. iz>=qa(3) .and. iz<=qb(3)) then
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
