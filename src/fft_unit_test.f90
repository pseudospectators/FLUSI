!-------------------------------------------------------------------------------
! FFT unit test
!-------------------------------------------------------------------------------
! Computes some derivatives that we can also compute exactly to ensure
! proper functioning of the FFT wrappers
! Input:
!       u: real valued work array
!       uk: complex valued work array
! Output:
!       none
! Side effects:
!       kills run if test fails (if the error is bigger than 1e-13)
!       says "hooray" if everything is fine
!-------------------------------------------------------------------------------
subroutine FFT_unit_test ( u, uk )
  use mpi
  use vars
  use p3dfft_wrapper
  ! input: real work array
  real(kind=pr),intent(inout):: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  ! input: complex work array
  complex(kind=pr), intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: iz,ix,iy
  logical :: fail
  real (kind=pr) :: kx, err, ky, kz, x, y, z
  if (mpirank==0) write(*,*) "----------------------"
  if (mpirank==0) write(*,*) "starting FFT unit test"
  fail = .false.

  !-----------------------------------------------------------------------------
  ! identity test
  !-----------------------------------------------------------------------------
  if (mpirank==0) write(*,*) "Testing if ifft(fft(u))=u ..."
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x = dble(ix)*dx
        y = dble(iy)*dy
        z = dble(iz)*dz
        ! Only one Fourier mode. Modified for the use of pruned FFT.
        u(ix,iy,iz) = dsin( y*2.d0*pi/yl ) * dcos( z*2.d0*pi/zl ) * dsin( x*2.d0*pi/xl )
      enddo
    enddo
  enddo
  call fft( inx=u, outk=uk )
  call ifft( ink=uk, outx=u )

  err = 0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x = dble(ix)*dx
        y = dble(iy)*dy
        z = dble(iz)*dz
        err = err + (u(ix,iy,iz) - dsin( y*2.d0*pi/yl ) * dcos( z*2.d0*pi/zl ) * dsin( x*2.d0*pi/xl ))**2
      enddo
    enddo
  enddo
  if (mpirank==0) then
    write(*,*) "FFT identity test done. error=", err
    if ( err > 1.0d-13 ) then
      write (*,*) "Very bad: FFT identity test failed."
      fail = .true.
    endif
  endif
  !-----------------------------------------------------------------------------
  ! derivative in X direction
  !-----------------------------------------------------------------------------
  if (nx>1) then
    ! fill real work array with a sine wave
    do ix=ra(1),rb(1)
      u(ix,:,:) = dsin( dble(ix)*dx*2.d0*pi/xl )
    enddo

    ! to fourier space
    call fft(u, uk)

    ! compute gradient
    do ix=ca(3), cb(3)
      kx = wave_x( ix )
      uk(:,:,ix) = uk(:,:,ix) * kx *dcmplx(0.d0,1.d0)
    enddo

    ! to x space
    call ifft(uk, u)

    ! compute error
    err = 0.d0
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          err = err + (u(ix,iy,iz)-dcos( dble(ix)*dx*2.d0*pi/xl )*2.d0*pi/xl)**2
        enddo
      enddo
    enddo
    err = err*dx*dy*dz

    if (mpirank==0) then
      write(*,*) "FFT unit (x) test done. error=", err
      if ( err > 1.0d-13 ) then
        write (*,*) "Very bad: FFT unit test failed."
        fail = .true.
      endif
    endif
  endif
  !-----------------------------------------------------------------------------
  ! derivative in Y direction
  !-----------------------------------------------------------------------------
  if (ny>1) then
    ! fill real work array with a sine wave
    do iy=ra(2),rb(2)
      u(:,iy,:) = dsin( dble(iy)*dy*2.d0*pi/yl )
    enddo

    ! to fourier space
    call fft(u, uk)

    ! compute gradient
    do iy=ca(2), cb(2)    ! ky : 0..ny/2-1 ,then, -ny/2..-1
      ky = wave_y(iy)
      uk(:,iy,:) = uk(:,iy,:)*ky*dcmplx(0.d0,1.d0)
    enddo

    ! to x space
    call ifft(uk, u)

    ! compute error
    err = 0.d0
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          err = err + (u(ix,iy,iz)-dcos( dble(iy)*dy*2.d0*pi/yl )*2.d0*pi/yl)**2
        enddo
      enddo
    enddo
    err = err*dx*dy*dz

    if (mpirank==0) then
      write(*,*) "FFT unit (y) test done. error=", err
      if ( err > 1.0d-13 ) then
        write (*,*) "Very bad: FFT unit test failed."
        fail = .true.
      endif
    endif
  endif
  !-----------------------------------------------------------------------------
  ! derivative in Z direction
  !-----------------------------------------------------------------------------
  if (nz>1) then
    ! fill real work array with a sine wave
    do iz=ra(3),rb(3)
      u(:,:,iz) = dsin( dble(iz)*dz*2.d0*pi/zl )
    enddo

    ! to fourier space
    call fft(u, uk)

    ! compute gradient
    do iz=ca(1),cb(1)  ! kz : 0..nz/2-1 ,then, -nz/2..-1
      kz = wave_z(iz)
      uk(iz,:,:) = uk(iz,:,:)*kz*dcmplx(0.d0,1.d0)
    enddo

    ! to x space
    call ifft(uk, u)

    ! compute error
    err = 0.d0
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          err = err + (u(ix,iy,iz)-dcos( dble(iz)*dz*2.d0*pi/zl )*2.d0*pi/zl)**2
        enddo
      enddo
    enddo
    err = err*dx*dy*dz

    if (mpirank==0) then
      write(*,*) "FFT unit (z) test done. error=", err
      if ( err > 1.0d-13 ) then
        write (*,*) "Very bad: FFT unit test failed."
        fail = .true.
      endif
    endif
  endif
  if (fail) call abort(20202,"At least one FFT test failed!")
  if (mpirank==0) write(*,*) "----------------------"
end subroutine FFT_unit_test
