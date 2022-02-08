! -------------------------------------------------------------------------
! Coherent Scalar Extraction of a 3D scalar field. The field is given
! by an NxNxN array, and the wavelet filter is given by the wavelet struct.
!
! Note N must be dyadic (16,32,64,128,256 etc) for the Fast Wavelet
! Transform
! -------------------------------------------------------------------------
! Algorithm:
!   In general, the coherent part is what remains after denoising. In order
!   to denoise, one has to know the variance of the noise (=incoherent
!   part), but that is a priori unknown. Therefore, an iterative method is
!   employed. We start assuming the entire signal is just noise, denoise,
!   and take the variance of what remains as better estimate.
! -------------------------------------------------------------------------
! Original version for 1D signals by Antonio Merulla, Michele Caldoro
! Modified for 3D CVE, 11/2015, T.Engels
! Ported to flusi (from matlab) 07/2017, T. Engels
! -------------------------------------------------------------------------
subroutine coherent_scalar_extraction( wavelet, maxiter, u_tot, u_coh, u_inc )
  implicit none
  type(orth_wavelet), intent(in) :: wavelet
  integer, intent(in) :: maxiter
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: u_tot, u_coh
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout), optional :: u_inc
  integer :: i, Jmin, nc
  real(kind=pr) :: Z, threshold_this, threshold_old, rate
  real(kind=pr) :: Z_tot, Z_coh, Z_inc
  real(kind=pr), dimension(:,:,:), allocatable :: wc, wc_tmp

  allocate( wc(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( wc_tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )

  ! forward wavelet transform, is executed only once
  if (mpirank==0) write(*,*) 'performing forward FWT of total field...'
  ! we set the minimum level to 1, but if the grid becomes odd-sized, we stop
  ! the wavelet transform, so Jmin is dictated also by FWT3_PO
  Jmin = 0
  call FWT3_PO(u_tot, wc, wavelet, Jmin, nc)
  ! call IWT3_PO(wc, u_coh, wavelet, Jmin, nc)

  ! the initial guess for the incoherent part is the entire signal itself
  ! (surely an overestimate)
  u_inc = u_tot

  ! loop variables initialization
  i=1
  threshold_this = 0.d0
  threshold_old = 9.0d9

  do while ( (abs(threshold_this-threshold_old)>1.0e-6*threshold_old) .and. (i<=maxiter) )
      ! current threshold becomes preceeding threshold
      threshold_old = threshold_this

      ! estimate new threshold (donoho threshold) using the current guess for
      ! the incoherent signal
      ! Eqn (6) in Kadoch et al "On the role of vortical structures for turbulent mixing using direct
      ! numerical simulation and wavelet-based coherent vorticity extraction"
      ! J. Turb.(12)20: 1-17 (2011)
      Z = sum( u_inc**2 ) / 2.d0
      Z = mpisum( Z )
      Z = Z / (dble(nx)*dble(ny)*dble(nz))
      threshold_this = sqrt( 4.d0 * Z * log( (dble(nx)*dble(ny)*dble(nz)) ) )

      ! apply the threshold, i.e. delete wavelet coefficients smaller than thresh_corr
      ! then inverse transform yields the coherent signal. substract this
      ! from the total signal to get the incoherent one
      call hard_thresholding( wc, wc_tmp, threshold_this)
      ! compute compression rate
      call compression_rate( wc_tmp, rate )
      ! coherent field is what remains after removing insignificant coeffs
      call IWT3_PO(wc_tmp, u_coh, wavelet, Jmin, nc)
      u_inc = u_tot - u_coh
      i = i+1

      if (mpirank==0) then
        write(*,'("Iteration ",i2," coherent part now contains ",f12.3,"% zeros, threshold=",g12.4)') i, 100.d0*rate, threshold_this
      endif
  end do

  if (mpirank==0) then
    write(*,*) 'Final threshold: ',threshold_old
    write(*,'("Did ",i3," iterations")') i
  endif

  Z_tot = mpisum(sum( u_tot**2 )) / 2.d0 / (dble(nx)*dble(ny)*dble(nz))
  Z_coh = mpisum(sum( u_coh**2 )) / 2.d0 / (dble(nx)*dble(ny)*dble(nz))
  Z_inc = mpisum(sum( u_inc**2 )) / 2.d0 / (dble(nx)*dble(ny)*dble(nz))

  ! compute (integral) enstrophy
  if (mpirank==0) then
    write(*,'("Z_tot = ", es15.8)') Z_tot
    write(*,'("Z_coh = ", es15.8)') Z_coh
    write(*,'("Z_inc = ", es15.8)') Z_inc
  endif

  deallocate(wc, wc_tmp)
end subroutine coherent_scalar_extraction



! -------------------------------------------------------------------------
! Coherent Vortex Extraction of a 3D vorticity field. The field is given
! by an NxNxNx3 array, and the wavelet filter is given by the wavelet struct.
!
! Note N must be dyadic (16,32,64,128,256 etc) for the Fast Wavelet
! Transform
! -------------------------------------------------------------------------
! Algorithm:
!   In general, the coherent part is what remains after denoising. In order
!   to denoise, one has to know the variance of the noise (=incoherent
!   part), but that is a priori unknown. Therefore, an iterative method is
!   employed. We start assuming the entire signal is just noise, denoise,
!   and take the variance of what remains as better estimate.
! -------------------------------------------------------------------------
! Original version for 1D signals by Antonio Merulla, Michele Caldoro
! Modified for 3D CVE, 11/2015, T.Engels
! Ported to flusi (from matlab) 07/2017, T. Engels
! -------------------------------------------------------------------------
subroutine coherent_vortex_extraction( wavelet, maxiter, vor_tot, vor_coh, vor_inc )
  implicit none
  type(orth_wavelet), intent(in) :: wavelet
  integer, intent(in) :: maxiter
  real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3), intent(inout) :: vor_tot, vor_coh, vor_inc
  integer :: i, k, Jmin, nc
  real(kind=pr) :: Z, threshold_this, threshold_old, rate(1:3)
  real(kind=pr) :: Z_tot, Z_coh, Z_inc, t0
  real(kind=pr), dimension(:,:,:,:), allocatable :: wc
  real(kind=pr), dimension(:,:,:), allocatable :: wc_tmp

  allocate( wc(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )
  allocate( wc_tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )

  ! forward wavelet transform, is executed only once
  if (mpirank==0) write(*,*) 'performing forward FWT of all components of input field...'
  ! we set the minimum level to 1, but if the grid becomes odd-sized, we stop
  ! the wavelet transform, so Jmin is dictated also by FWT3_PO
  Jmin = 1

  t0 = MPI_wtime()
  do k = 1,3
    call FWT3_PO(vor_tot(:,:,:,k), wc(:,:,:,k), wavelet, Jmin, nc)
  enddo
  if (mpirank==0) write(*,*) 'Done. Elapsed time: ', MPI_wtime()-t0

  ! the initial guess for the incoherent part is the entire signal itself
  ! (surely an overestimate)
  vor_inc = vor_tot

  ! loop variables initialization
  i=1
  threshold_this = 0.d0
  threshold_old = 9.0d9

  do while ( (abs(threshold_this-threshold_old)>1.0e-6*threshold_old) .and. (i<=maxiter) )
      ! current threshold becomes preceeding threshold
      threshold_old = threshold_this

      ! estimate new threshold (donoho threshold) using the current guess for
      ! the incoherent signal
      ! Eqn (6) in Kadoch et al "On the role of vortical structures for turbulent mixing using direct
      ! numerical simulation and wavelet-based coherent vorticity extraction"
      ! J. Turb.(12)20: 1-17 (2011)
      wc_tmp = vor_inc(:,:,:,1)**2 + vor_inc(:,:,:,2)**2 + vor_inc(:,:,:,3)**2
      Z = sum( wc_tmp ) / 2.d0
      Z = mpisum( Z )
      Z = Z / (dble(nx)*dble(ny)*dble(nz))
      threshold_this = sqrt( (4.d0/3.d0) * Z * log( (dble(nx)*dble(ny)*dble(nz)) ) )

      do k = 1, 3
        ! apply the threshold, i.e. delete wavelet coefficients smaller than thresh_corr
        ! then inverse transform yields the coherent signal. substract this
        ! from the total signal to get the incoherent one
        call hard_thresholding( wc(:,:,:,k), wc_tmp, threshold_this )
        ! compute compression rate
        call compression_rate( wc_tmp, rate(k) )
        ! coherent field is what remains after removing insignificant coeffs
        call IWT3_PO( wc_tmp, vor_coh(:,:,:,k), wavelet, Jmin, nc )
        ! compute incoherent part
        vor_inc(:,:,:,k) = vor_tot(:,:,:,k) - vor_coh(:,:,:,k)
      enddo

      i = i+1

      if (mpirank==0) then
        write(*,'("Iteration ",i2," coherent part contains ",3(f12.3,"% ")," zeros, threshold=",g12.4)') &
        i, 100.d0*rate, threshold_this
      endif
  end do

  if (mpirank==0) then
    write(*,*) 'Final threshold: ', threshold_this
    write(*,'("Did ",i3," iterations")') i
  endif

  ! compute (integral) enstrophy
  wc_tmp = vor_tot(:,:,:,1)**2 + vor_tot(:,:,:,2)**2 + vor_tot(:,:,:,3)**2
  Z_tot = sum( wc_tmp ) / 2.d0
  Z_tot = mpisum( Z_tot )
  Z_tot = Z_tot / (dble(nx)*dble(ny)*dble(nz))

  wc_tmp = vor_coh(:,:,:,1)**2 + vor_coh(:,:,:,2)**2 + vor_coh(:,:,:,3)**2
  Z_coh = sum( wc_tmp ) / 2.d0
  Z_coh = mpisum( Z_coh )
  Z_coh = Z_coh / (dble(nx)*dble(ny)*dble(nz))

  wc_tmp = vor_inc(:,:,:,1)**2 + vor_inc(:,:,:,2)**2 + vor_inc(:,:,:,3)**2
  Z_inc = sum( wc_tmp ) / 2.d0
  Z_inc = mpisum( Z_inc )
  Z_inc = Z_inc / (dble(nx)*dble(ny)*dble(nz))

  if (mpirank==0) then
    write(*,'("Z_tot = ", es15.8)') Z_tot
    write(*,'("Z_coh = ", es15.8)') Z_coh
    write(*,'("Z_inc = ", es15.8)') Z_inc
  endif

  deallocate(wc, wc_tmp)
end subroutine coherent_vortex_extraction
