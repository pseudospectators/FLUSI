!-------------------------------------------------------------------------------
! Prescribed beam deformation (given fseries of theta)
! read the Fourier coefficients for theta from disk (only once of course)
! and for any time t evaluate theta, theta_dot from the Fourier series
! then, integrate the position vector and be happy.
! this is for PASSIVE FSI
! To DO: read nfft from file as well
!        read filenames from params file
!-------------------------------------------------------------------------------
subroutine prescribed_beam ( time, dt, beam_solid )! note this is actuall only ONE beam!!
  implicit none
  real(kind=pr), intent (in) :: dt, time
  type(solid), intent(inout) :: beam_solid
  real(kind=pr),dimension(:,:),allocatable,save :: data_ai, data_bi
  integer :: i, mpicode, ns_file
  integer, save :: nfft

  !-----------------------------------------------------------------------------
  ! on first coll, allocate memory and read data
  !-----------------------------------------------------------------------------
  if (.not. allocated(data_ai)) then
    ! initialization code

    if (mpirank==0) then
      ! read from disk.
      open(37, file="beam_theta_prescribed_ai.in", form='formatted', status='old')
      open(38, file="beam_theta_prescribed_bi.in", form='formatted', status='old')

      read(37,*) ns_file, nfft
      write(*,'("beam_theta_prescribed_ai.in: ns=",i3," nfft=",i3)') ns_file, nfft
      read(38,*) ns_file, nfft
      write(*,'("beam_theta_prescribed_bi.in: ns=",i3," nfft=",i3)') ns_file, nfft

      if (ns_file /= ns) then
        call abort(99443,"beam resolution and prescribed.in resolution different")
      endif
    endif

    call MPI_BCAST(nfft,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode)

    allocate(data_ai(0:ns-1,nfft+1)) !+1 for zero mode
    allocate(data_bi(0:ns-1,nfft))

    if (mpirank==0) then
      do i = 0, ns-1 ! beam points are in rows
        read(37,*) data_ai(i,:)
        read(38,*) data_bi(i,:)
      enddo

      close(37)
      close(38)
    endif
    call MPI_BCAST(data_ai,ns*(nfft+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
    call MPI_BCAST(data_bi,ns*nfft,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
  endif

  !-----------------------------------------------------------------------------
  ! reading done, evaluate Fourier series for all points on the beam
  !-----------------------------------------------------------------------------
  do i=0,ns-1
    call fseries_eval(time+dt, beam_solid%theta(i) , beam_solid%theta_dot(i) ,&
                      data_ai(i,1), data_ai(i,2:nfft+1), data_bi(i,1:nfft))
  enddo

  ! note in the "true" solvers, the beam in the routine is advanced from t to
  ! t_n+1, so this is what we do here, too.
  call integrate_position ( time+dt, beam_solid )

end subroutine prescribed_beam
