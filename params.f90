subroutine params (tdrag, ifield, infile, tsave, ismth, &
     tstart, nt)
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none

  integer, intent (out) :: ismth, nt
  integer, dimension (3), intent (out) :: ifield
  real (kind=pr), intent (out)         :: tdrag, tsave, tstart
  integer                              :: mpicode
  integer, dimension (11)              :: comm_int
  real (kind=pr), dimension (11)       :: comm_real
  character (len=80)                   :: infile

  !-- Read parameters used for computations 
  !-- both with and without obstacles
  ! Node 0 reads stdin
  if ( mpirank == 0 ) then
     read *, nx
     read *, ny
     read *, nz
     read *, nt
     read *, tmax
     read *, cfl
     read *, tsave
     read *, ifield (1)
     read *, ifield (2)
     read *, ifield (3)
     read *, xl
     read *, yl
     read *, zl
     read *, nu
     read *, idealis
     read *, ihypvisc
     read *, iobst
     read *, inicond
     read *, tdrag

     comm_int(1:11) = (/ nx, ny, nz, nt, ifield(1:3), idealis, ihypvisc, iobst, inicond /)
     comm_real(1:8) = (/ tmax, cfl, tsave, xl, yl, zl, nu, tdrag /)
  endif

  ! Broadcast to other processes
  call MPI_BCAST( comm_int, 11, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( comm_real, 8, mpireal, 0, MPI_COMM_WORLD, mpicode )

  if ( mpirank /= 0 ) then
     nx = comm_int(1)
     ny = comm_int(2)
     nz = comm_int(3)
     nt = comm_int(4)
     tmax = comm_real(1)
     cfl = comm_real(2)
     tsave = comm_real(3)
     ifield (1) = comm_int(5)
     ifield (2) = comm_int(6)
     ifield (3) = comm_int(7)
     xl = comm_real(4)
     yl = comm_real(5)
     zl = comm_real(6)
     nu = comm_real(7)
     idealis = comm_int(8)
     ihypvisc = comm_int(9)
     iobst = comm_int(10)
     inicond = comm_int(11)
     tdrag = comm_real(8)
  endif

  !-- Stop if the sizes are odd or smaller than 4
  if ( nx<4 .or. ny<4 .or. nz<4 .or. modulo(nx,2)==1 .or. modulo(ny,2)==1 .or. modulo(nz,2)==1 ) then
    if ( mpirank == 0 ) then
      print *, 'nx, ny, nz must be even and not smaller than 4'
    endif
    stop
  endif

  !-- Set up various parameters
  pi     = 4.0 * atan (1.0)
  scalex = 2.0*pi / xl
  scaley = 2.0*pi / yl
  scalez = 2.0*pi / zl
  
  dx = xl / dble (nx)
  dy = yl / dble (ny)
  dz = zl / dble (nz)

  !-- Initialize fft
  call fft_initialize 

  !-- Set up mask for dealiasing
  !!!!!! TO DO !!!!!!!
!  allocate ( dealiase (0:nx-1, 0:ny-1, 0:nz-1) )
!  dealiase = 1.0d0
!  if (idealis == 1) call dealiase_mask

  !-- Initialize variables in case of no obstacles
  x0    = 0.0
  y0    = 0.0
  z0    = 0.0
  eps   = 0.0
  size  = 0.0
  Ux    = 0.0
  Uy    = 0.0
  Uz    = 0.0
  Ax    = 0.0
  Ay    = 0.0
  Az    = 0.0
  imask = 0
  ismth = 0
  imove = 0

  !-- Set up parameters for computations
  !-- with obstacles 
  if (iobst > 0) then

     ! Node 0 reads stdin
     if ( mpirank == 0 ) then
        read *, imask
        read *, ismth
        read *, imove
        read *, eps
        read *, x0
        read *, y0
        read *, z0
        read *, size
        read *, Ux
        read *, Uy
        read *, Uz
        read *, Ax
        read *, Ay
        read *, Az

        comm_int(1:3) = (/ imask, ismth, imove /)
        comm_real(1:11) = (/ eps, x0, y0, z0, size, Ux, Uy, Uz, Ax, Ay, Az /)
     end if     

     ! Broadcast to other processes
     call MPI_BCAST( comm_int, 3, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
     call MPI_BCAST( comm_real, 11, mpireal, 0, MPI_COMM_WORLD, mpicode )

     if ( mpirank /= 0 ) then
        imask = comm_int(1)
        ismth = comm_int(2)
        imove = comm_int(3)
        eps = comm_real(1)
        x0 = comm_real(2)
        y0 = comm_real(3)
        z0 = comm_real(4)
        size = comm_real(5)
        Ux = comm_real(6)
        Uy = comm_real(7)
        Uz = comm_real(8)
        Ax = comm_real(9)
        Ay = comm_real(10)
        Az = comm_real(11)
     endif

  endif

  !-- Set up initial conditions
  tstart = 0.0 ! starting time
  
  iSaveVelocity  = ifield(1)
  iSaveVorticity = ifield(2)
  iSavePress     = ifield(3)
  iSaveMask 	 = 1

  if (inicond == 3) then
     if ( mpirank == 0 ) then
       read *, tstart
       read *, infile
     endif

     ! Broadcast to other processes
     call MPI_BCAST( tstart, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
     call MPI_BCAST( infile, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpicode )
  end if

  !-- Print out header
  if ( mpirank == 0 ) then
     write (*,*) '--------------------------------------------'
     write (*, '(" nx = ",i4, 1x, "ny = ", i4, 1x, "nz = ", i4)') nx, ny, nz
     write (*, '(" xl = ", f6.2, 1x, "yl = ", f6.2, 1x, "zl = ", f6.2)') xl, yl, zl
     write (*, '(" nt = ", i6, 1x)') nt
     write (*, '(" tmax = ", f6.2)') tmax
     write (*, '(" viscosity = ", es11.4)') nu
     write (*,*) '--------------------------------------------'
     write (*, '(" calculate drag every", es11.4, " time steps")') tdrag
     write (*, '(" save fields every   ", es11.4)') tsave
     if (ifield (1) == 1) write (*,*) 'save pressure'
     if (ifield (2) == 1) write (*,*) 'save velocity'
     if (ifield (3) == 1) write (*,*) 'save vorticity'
     write (*,*) '--------------------------------------------'
     if (iobst > 0) then 
        write (*,*) 'with obstacle'
        write (*, '(" x0 = ", f6.2, 1x, "y0 = ", f6.2, 1x, "z0 = ", f6.2, 1x, &
            & "size = ", f6.2)') x0, y0, z0, size
        write (*, '(" epsilon_eff   = ", es11.4)') eps * xl / real (nx)
        write (*, '(" epsilon       = ", es11.4)') eps
        write (*, '(" Ux = ", f6.2, " Uy = ", f6.2, " Uz = ", f6.2)') Ux, Uy, Uz
        write (*, '(" Ax = ", f6.2, " Ay = ", f6.2, " Az = ", f6.2)') Ax, Ay, Az
        write (*, '(" Re = ", es11.4)') sqrt(Ux**2+Uy**2+Uz**2) * size / nu
     else
        write (*,*) 'without obstacle'
     end if
     write (*,*) ' '
     write (*,*) '--------------------------------------------'
     write (*,*) ' '
  endif

end subroutine params









