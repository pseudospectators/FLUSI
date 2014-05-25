!-------------------------------------------------------------------------------
! right hand side for passive scalars. all input arrays are scalars (and not 
! vectors), in future revisions this will be called for *each* passive scalar
! while for starters we have only one.
!-------------------------------------------------------------------------------
subroutine cal_nlk_scalar( time, it, u, uk, nlk, workc1, work )
  use mpi
  use p3dfft_wrapper
  use fsi_vars
  implicit none

  real(kind=pr),intent(in) :: time
  integer, intent(in) :: it
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::workc1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)) 
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  real(kind=pr) :: kx,ky,kz,chi,w,D,maxi,t1
  real(kind=pr) :: k(1:3) ! wavenumber vector
  complex(kind=pr) :: imag
  integer :: ix,iy,iz,id,mpicode
  t1 = MPI_wtime()
  
  !-----------------------------------------------------------------------------
  ! Sanity test: if the scalar values are out of range, we skip the entire 
  ! right hand side, but we do not abort (since the fluid is fine, it's just
  ! the scalar that fails.
  !-----------------------------------------------------------------------------
  if ((modulo(it,10)==0).and.(compute_scalar)) then
     call ifft( ink=uk, outx=work(:,:,:,1) )
     !-- global maximum
     call MPI_ALLREDUCE(maxval(work(:,:,:,1)),maxi,1,MPI_DOUBLE_PRECISION,&
          MPI_MAX,MPI_COMM_WORLD,mpicode)
     !-- skip passive scalar from now on, if out of range    
     if ( maxi > 1.0e2 ) then
        if (mpirank==0) then
          write(*,'("WARNING! PASSIVE SCALAR FAILED! time=",es12.4)') time
          write(*,'("THIS MEANS WE WILL NO LONGER COMPUTE IT FROM NOW ON!!!")')
        endif
        !-- kill run if scalar is considered crutial for results
        if (stop_on_fail=="yes") stop
        !-- otherwise, warn and continue without scalar
        compute_scalar = .false. ! callees will skip
     endif    
  endif
  
  !-- initialization
  nlk = dcmplx(0.d0,0.d0)
  imag = dcmplx(0.d0,1.d0)
  
  !------------------------------------------------------------------------------
  ! compute right hand side of passive scalar, component by component
  !------------------------------------------------------------------------------  
  do id = 1,3
    ! step 1: compute id-component of gradient theta in f-space
    do iz=ca(1),cb(1)
      kz=wave_z(iz)
      do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
          kx=wave_x(ix)
          k = (/kx, ky ,kz /)
          workc1(iz,iy,ix) = uk(iz,iy,ix) * imag * k(id)
        enddo
      enddo
    enddo
    
    ! step 2: 
    call ifft ( ink=workc1, outx=work(:,:,:,1))
    work(:,:,:,2) = work(:,:,:,1)
    
    ! step 3
    do ix=ra(1),rb(1)
      do iy=ra(2),rb(2)
        do iz=ra(3),rb(3) 
          chi = mask(ix,iy,iz)*eps
          w = -( (1.d0-chi)*u(ix,iy,iz,id) + chi*us(ix,iy,iz,id) )
          work(ix,iy,iz,1) = work(ix,iy,iz,1)*w
        enddo
      enddo
    enddo
    ! step 4
    call fft ( inx=work(:,:,:,1), outk=workc1 )
    ! step 5
    nlk = nlk + workc1
    
    !-------------------------------------------------------------
    

    ! step 6
    do ix=ra(1),rb(1)
      do iy=ra(2),rb(2)
        do iz=ra(3),rb(3) 
          chi = mask(ix,iy,iz)*eps
          !-- inside obstacle, reduced fluxes
          D = kappa*(1.d0-chi) + eps_scalar*chi
          work(ix,iy,iz,2) = work(ix,iy,iz,2) * D
        enddo
      enddo
    enddo
    
    ! step 7
    call fft( inx=work(:,:,:,2), outk=workc1 )
    
    ! step 8: compute id-component of gradient theta in f-space
    do iz=ca(1),cb(1)
      kz=wave_z(iz)
      do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
          kx=wave_x(ix)
          k = (/kx, ky ,kz /)
          workc1(iz,iy,ix) = workc1(iz,iy,ix) * imag * k(id)
        enddo
      enddo
    enddo
    ! step 9
    nlk = nlk + workc1   
  enddo
  time_scalar = time_scalar + MPI_wtime() - t1
end subroutine cal_nlk_scalar