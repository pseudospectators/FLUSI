!-----------------------------------------------------------------------------
! 3D fast wavelet transform, periodized, orthogonal
!-----------------------------------------------------------------------------
subroutine IWT3_PO(wc, u, wavelet, Jmin, ncoarse)
  implicit none
  ! derived datatype containing the filters
  type(orth_wavelet),intent(in) :: wavelet
  ! wavelet coefficients
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
  ! output of reconstruction
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: u
  ! coarse level and number of points on coarse level
  integer, intent(inout) :: Jmin, ncoarse
  ! signal length
  integer :: N, J0, j, ix, iy, iz, idir, mpicode, nc
  ! buffer
  real(kind=pr), allocatable, dimension(:) :: buffer1, buffer2, buffer3
  real(kind=pr), allocatable, dimension(:,:,:) :: work
  integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst

  N = nx ! this works currently only for nx = ny = nz
  J0 = ceiling( log(dble(N)) / log(2.d0) )

  ! the first iteration acts on the coarest level (given in argument, but do not
  ! modify its value! copy here for case when several IFWT are performed and only
  ! one FWT)
  ! nc = 2**(Jmin+1)
  nc = ncoarse
  ! initialize wavelet coefficients
  u(:,:,:) = wc(:,:,:)


  if (nx /= ny .or. ny /= nz .or. nx /= nz) then
    call abort(55527,"Wavelet error: condition nx=ny=nz is violated.")
  endif

  ! -------------------------------------------------------
  ! IWT starts from coarse level and goes up level by level
  ! -------------------------------------------------------
  do j = Jmin, J0-1
      ! an inplace filtering is not possible, so we allocate buffers for a line
      ! (of current size)
      allocate( buffer1(1:nc), buffer2(1:nc), buffer3(1:nc) )

      ! ---------------------
      ! --- x direction -----
      ! ---------------------
      ! Filters are applied in the first coordinate, as it is always contiguous in memory
      ! and not split among MPI processes
      call WaveReconstruction_dim1( u, ra, rb, nc, wavelet, buffer1, buffer2, buffer3 )

      ! ---------------------
      ! --- y direction -----
      ! ---------------------
      if (mpidims(2) == 1) then
        ! CASE A): index 2 is contiguous and not split among procs.
        call WaveReconstruction_dim2( u, ra, rb, nc, wavelet, buffer1, buffer2, buffer3 )

      else
        ! transposition: exchange x-y data
        ! idir = 1 - permute 1 and 2 indices
        ! idir = 2 - permute 1 and 3 indices
        idir = 1
        call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
        allocate( work(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
        work = 0.d0
        call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                     MPI_DOUBLE_PRECISION, u, work )

        ! at this point we have work(iy,ix,iz)
        ! again, the first dimension is filtered, but this now is the y-data
        call WaveReconstruction_dim1( work, kat, kbt, nc, wavelet, buffer1, buffer2, buffer3 )

        ! transposition: exchange x-y data
        ! idir = 1 - permute 1 and 2 indices
        ! idir = 2 - permute 1 and 3 indices
        idir = 1
        call trextents(idir, (/ny,nx,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
        call subtr ( idir, (/ny,nx,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                     MPI_DOUBLE_PRECISION, work, u )
        ! at this point we have again u(ix,iy,iz), which is what we started with
      endif


      idir = 2
      call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
      if(allocated(work)) deallocate(work)
      allocate( work(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
      call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                   MPI_DOUBLE_PRECISION, u, work )

      ! ---------------------
      ! --- z direction -----
      ! ---------------------
      ! at this point we have work(iz,iy,ix)
      ! again, the first dimension is filtered, but this now is the z-data
      call WaveReconstruction_dim1( work, kat, kbt, nc, wavelet, buffer1, buffer2, buffer3 )

      ! final permutation back to the original dataset
      idir = 2
      call trextents(idir, (/nz,ny,nx/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
      call subtr ( idir, (/nz,ny,nx/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                   MPI_DOUBLE_PRECISION, work, u )
      if(allocated(work)) deallocate(work)
      ! at this point we have u(ix,iy,iz)

      nc = nc * 2
      deallocate(buffer1, buffer2, buffer3)
  enddo

end subroutine IWT3_PO


!-----------------------------------------------------------------------------
! Reconstruction from low- and high pass filtered coefficients.
! Data is first upsampled, then filtered with the reconstruction filters.
! Note reconstruction filters are reverse of decomposition filters.
!
!         | low pass     ||  high pass   |
! Input:  hhhhhhhhhhhhhhhhgggggggggggggggg
! Output: uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
!-----------------------------------------------------------------------------
subroutine WaveReconstruction_dim1( wc, qa, qb, nc, wavelet, buffer1, buffer2, buffer3 )
  implicit none
  integer, intent(in) :: nc
  ! note this routine is one of the few where the decomposition has to be general
  integer, dimension(1:3), intent(in) :: qa, qb
  real(kind=pr), dimension(qa(1):qb(1),qa(2):qb(2),qa(3):qb(3)), intent(inout) :: wc
  real(kind=pr), dimension(1:nc), intent(inout) :: buffer1, buffer2, buffer3
  type(orth_wavelet), intent(in) :: wavelet
  integer :: ix,iy,iz

  ! Filters are applied in the first coordinate, as it is always contiguous in memory
  ! and not split among MPI processes
  do iy = 0, nc-1 ! note loop bounds 0,nc-1 (and not 0,ny-1)
    do iz = 0, nc-1
      ! is this on my local bounds? (MPI-sense, some procs may idle)
      if (iy>=qa(2) .and. iy<=qb(2) .and. iz>=qa(3) .and. iz<=qb(3)) then
        ! bot: 0:nc/2-1 (zero-based)
        ! top: nc/2:nc-1 (zero-based)

        ! fill upsampling buffer for low-pass filter: every second point
        buffer3 = 0.d0
        buffer3(1:nc-1:2) = wc(0:nc/2-1,iy,iz) ! bot, UpDyadLo
        ! apply low-pass filter to upsampled signal
        call filter1P( buffer3, buffer1, wavelet%HR, lbound(wavelet%HR,dim=1), ubound(wavelet%HR,dim=1))

        ! fill upsampling buffer for high-pass filter: every second point
        buffer3 = 0.d0
        buffer3(1:nc-1:2) = wc(nc/2:nc-1,iy,iz) ! top, UpDyadHi
        ! apply high-pass filter to upsampled signal
        call filter1P( buffer3, buffer2, wavelet%GR, lbound(wavelet%GR,dim=1), ubound(wavelet%GR,dim=1))

        wc(0:nc-1,iy,iz) = buffer1 + buffer2
      endif
    enddo
  enddo
end subroutine WaveReconstruction_dim1

!-----------------------------------------------------------------------------
! Reconstruction from low- and high pass filtered coefficients.
! Data is first upsampled, then filtered with the reconstruction filters.
! Note reconstruction filters are reverse of decomposition filters.
!
! Routine works along second index of wc, requiring it to be contiguous
!
!         | low pass     ||  high pass   |
! Input:  hhhhhhhhhhhhhhhhgggggggggggggggg
! Output: uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
!-----------------------------------------------------------------------------
subroutine WaveReconstruction_dim2( wc, qa, qb, nc, wavelet, buffer1, buffer2, buffer3 )
  implicit none
  integer, intent(in) :: nc
  ! note this routine is one of the few where the decomposition has to be general
  integer, dimension(1:3), intent(in) :: qa, qb
  real(kind=pr), dimension(qa(1):qb(1),qa(2):qb(2),qa(3):qb(3)), intent(inout) :: wc
  real(kind=pr), dimension(1:nc), intent(inout) :: buffer1, buffer2, buffer3
  type(orth_wavelet), intent(in) :: wavelet
  integer :: ix,iy,iz

  ! Filters are applied in the first coordinate, as it is always contiguous in memory
  ! and not split among MPI processes
  do ix = 0, nc-1 ! note loop bounds 0,nc-1 (and not 0,ny-1)
    do iz = 0, nc-1
      ! is this on my local bounds? (MPI-sense, some procs may idle)
      if (iy>=qa(2) .and. iy<=qb(2) .and. iz>=qa(3) .and. iz<=qb(3)) then
        ! bot: 0:nc/2-1 (zero-based)
        ! top: nc/2:nc-1 (zero-based)

        ! fill upsampling buffer for low-pass filter: every second point
        buffer3 = 0.d0
        buffer3(1:nc-1:2) = wc(ix,0:nc/2-1,iz) ! bot, UpDyadLo
        ! apply low-pass filter to upsampled signal
        call filter1P( buffer3, buffer1, wavelet%HR, lbound(wavelet%HR,dim=1), ubound(wavelet%HR,dim=1))

        ! fill upsampling buffer for high-pass filter: every second point
        buffer3 = 0.d0
        buffer3(1:nc-1:2) = wc(ix,nc/2:nc-1,iz) ! top, UpDyadHi
        ! apply high-pass filter to upsampled signal
        call filter1P( buffer3, buffer2, wavelet%GR, lbound(wavelet%GR,dim=1), ubound(wavelet%GR,dim=1))

        wc(ix,0:nc-1,iz) = buffer1 + buffer2
      endif
    enddo
  enddo
end subroutine WaveReconstruction_dim2
