module helpers
  use vars
  implicit none

  ! module global variables


contains

  !-----------------------------------------------------------------------------
  ! This function returns, to a given filename, the corresponding dataset name
  ! in the hdf5 file, following flusi conventions (folder/ux_0000.h5 -> "ux")
  !-----------------------------------------------------------------------------
  character(len=strlen)  function get_dsetname(fname)
    implicit none
    character(len=*), intent(in) :: fname
    ! extract dsetname (from "/" until "_", excluding both)
    get_dsetname  = fname  ( index(fname,'/',.true.)+1:index( fname, '_',.true. )-1 )
    return
  end function get_dsetname



  !-------------------------------------------------------------------------------
  ! evaluate a fourier series given by the coefficents a0,ai,bi
  ! at the time "time", return the function value "u" and its
  ! time derivative "u_dt". Uses assumed-shaped arrays, requires an interface.
  !-------------------------------------------------------------------------------
  subroutine fseries_eval(time,u,u_dt,a0,ai,bi)
    use vars
    implicit none

    real(kind=pr), intent(in) :: a0, time
    real(kind=pr), intent(in), dimension(:) :: ai,bi
    real(kind=pr), intent(out) :: u, u_dt
    real(kind=pr) :: c,s,f
    integer :: nfft, i

    nfft=size(ai)

    ! frequency factor
    f = 2.d0*pi

    u = 0.5d0*a0
    u_dt = 0.d0

    do i=1,nfft
        s = dsin(f*dble(i)*time)
        c = dcos(f*dble(i)*time)
        ! function value
        u    = u + ai(i)*c + bi(i)*s
        ! derivative (in time)
        u_dt = u_dt + f*dble(i)*(-ai(i)*s + bi(i)*c)
    enddo
  end subroutine fseries_eval


  !-------------------------------------------------------------------------------
  ! evaluate hermite series, given by coefficients ai (function values)
  ! and bi (derivative values) at the locations x. Note that x is assumed periodic;
  ! do not include x=1.0.
  ! a valid example is x=(0:N-1)/N
  !-------------------------------------------------------------------------------
  subroutine hermite_eval(time,u,u_dt,ai,bi)
    use vars
    implicit none

    real(kind=pr), intent(in) :: time
    real(kind=pr), intent(in), dimension(:) :: ai,bi
    real(kind=pr), intent(out) :: u, u_dt
    real(kind=pr) :: dt,h00,h10,h01,h11,t
    integer :: n, j1,j2

    n=size(ai)

    dt = 1.d0 / dble(n)
    j1 = floor(time/dt) + 1
    j2 = j1 + 1
    ! periodization
    if (j2 > n) j2=1
    ! normalized time (between two data points)
    t = (time-dble(j1-1)*dt) /dt

    ! values of hermite interpolant
    h00 = (1.d0+2.d0*t)*((1.d0-t)**2)
    h10 = t*((1.d0-t)**2)
    h01 = (t**2)*(3.d0-2.d0*t)
    h11 = (t**2)*(t-1.d0)

    ! function value
    u = h00*ai(j1) + h10*dt*bi(j1) &
      + h01*ai(j2) + h11*dt*bi(j2)

    ! derivative values of basis functions
    h00 = 6.d0*t**2 - 6.d0*t
    h10 = 3.d0*t**2 - 4.d0*t + 1.d0
    h01 =-6.d0*t**2 + 6.d0*t
    h11 = 3.d0*t**2 - 2.d0*t

    ! function derivative value
    u_dt = (h00*ai(j1) + h10*dt*bi(j1) &
         + h01*ai(j2) + h11*dt*bi(j2) ) / dt
  end subroutine hermite_eval


function mpisum( a )
  use mpi
  use vars
  implicit none
  real(kind=pr) :: a_loc, mpisum
  real(kind=pr),intent(in) :: a
  integer :: mpicode
  a_loc=a
  call MPI_ALLREDUCE (a_loc,mpisum,1, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
end function

end module helpers
