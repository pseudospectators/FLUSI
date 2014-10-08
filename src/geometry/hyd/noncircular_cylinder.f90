!-------------------------------------------------------------------------------
! Non-circular cylinders.
!-------------------------------------------------------------------------------
! The radius is described as a function of theta by Fourier coefficients, i.e.
! R(theta) = a0 + ai*sin() + bi*cos()
! It is intended to run this subroutine only once (or twice, but not in every 
! time step, which is why we read the parameters from the shape.in file here
! (as opposed to do so once on initialization)
!-------------------------------------------------------------------------------
subroutine noncircular_cylinder(mask,mask_color,us)
  use vars
  implicit none
  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  integer :: iz, iy, nfft, mpicode,k
  real(kind=pr) :: y,z,tmp, R, N_smooth, a0, theta, R0, safety
  real(kind=pr),dimension(:),allocatable :: ai,bi

  !-----------------------------------------------------------------------------
  ! Phase 1: read in shape information from file
  !-----------------------------------------------------------------------------
  ! learn how many Fourier coefficients to expect
  if (mpirank==0) then
    ! open input file
    call check_file_exists('shape.in')
    open(37, file='shape.in', form='formatted', status='old')
    read(37,*) nfft    
    write(*,'("Non-circular cylinder with nfft=",i2)') nfft
  endif
  ! BCAST nfft to all procs
  call MPI_BCAST( nfft,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
  ! allocate fourier coefficent arrays
  allocate(ai(1:nfft),bi(1:nfft))
  
  if (mpirank==0) then
    read(37,*) a0
    read(37,*) ai
    read(37,*) bi
    write(*,*) ai
    write(*,*) bi
    close(37)
  endif
  
  call MPI_BCAST( a0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
  call MPI_BCAST( ai,nfft,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
  call MPI_BCAST( bi,nfft,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode )
  
  
  !-----------------------------------------------------------------------------
  ! Phase 2: construct mask function
  !-----------------------------------------------------------------------------
  N_smooth = 1.5d0
  safety = 2.d0*N_smooth*max(dz,dy)

  do iz=ra(3), rb(3)
    do iy=ra(2),rb(2)    
      y = dble(iy)*dy
      z = dble(iz)*dz
      R = dsqrt( (y-y0)**2 + (z-z0)**2 )
      
      ! bounding box check
      if (R <= 1.5d0) then
        theta = pi + atan2(y-y0,z-z0)
        ! construct R0(theta)
        R0 = a0/2.d0
        do k = 1,nfft
          R0 = R0 + ai(k)*dcos(dble(k)*theta) + bi(k)*dsin(dble(k)*theta)
        enddo
        
        if ( R <= R0+safety ) then
          call SmoothStep (tmp, R, R0, N_smooth*max(dz,dy))
          mask (:,iy,iz) = tmp

          ! assign color "1" where >0 indicates something "useful"
          if (tmp > 1.0e-12) mask_color(:,iy,iz) = 1
        endif
      endif
      enddo
  enddo
  
  deallocate(ai,bi)
end subroutine noncircular_cylinder