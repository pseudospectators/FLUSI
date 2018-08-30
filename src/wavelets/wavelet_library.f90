module flusi_wavelet_lib

  use vars, only : pr, ra, rb, nx, ny, nz, on_proc, mpidims, mpicoords, mpicommslab, abort, mpirank, root
  use mpi
  use p3dfft_wrapper
  use module_helpers, only : mpisum

  implicit none

  ! derived datatype for the orthogonal wavelet. note it is defined by four filter banks
  ! (high/low pass or decomposition/reconstruction), and we give it a name
  type orth_wavelet
    real(kind=pr), dimension(:), allocatable :: HD ! low pass decomposition filter
    real(kind=pr), dimension(:), allocatable :: GD ! high pass decomposition filter
    real(kind=pr), dimension(:), allocatable :: HR ! low pass reconstruction filter
    real(kind=pr), dimension(:), allocatable :: GR ! high pass reconstruction filter
    character(len=80) :: name
  end type

contains

  ! include files, note they have to be marked as dependencies in the makefile
  include "FWT3_PO.f90"
  include "IWT3_PO.f90"
  include "coherent_vortex_extraction.f90"


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

    ! ****** decomposition filter ******

    ! low pass decomposition filter (in WaveLab 850, this is QMF)
    wavelet%HD = (/1.638733646318000d-02, -4.146493678197000d-02,-6.737255472230000d-02,&
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

    ! ****** reconstruction filter ******

    allocate( wavelet%HR(-7:4) )
    allocate( wavelet%GR(-5:6) )

    ! Low-pass reconstruction filter (same as decomposition above)
    wavelet%HR = (/1.638733646318000d-02, -4.146493678197000d-02,-6.737255472230000d-02,&
    3.861100668230900d-01,8.127236354496100d-01,4.170051844237800d-01,-7.648859907826000d-02,&
    -5.943441864647000d-02,2.368017194688000d-02,5.611434819370000d-03,-1.823208870910000d-03,&
    -7.205494453700000d-04/)

    ! high-pass reconstruction filter
    wavelet%GR = (/7.205494453700000d-04,-1.823208870910000d-03,-5.611434819370000d-03,&
    2.368017194688000d-02,5.943441864647000d-02,-7.648859907826000d-02,-4.170051844237800d-01,&
    8.127236354496100d-01,-3.861100668230900d-01,-6.737255472230000d-02,4.146493678197000d-02,&
    1.638733646318000d-02/)


    call reverse( wavelet%HR )
    call reverse( wavelet%GR )
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

  subroutine dyad(J, i1, i2)
    implicit none
    integer, intent(in) :: J
    integer, intent(out) :: i1, i2

    ! dyad -- Index entire j-th dyad of 1-d wavelet xform
    !  Usage
    !    ix = dyad(j);
    !  Inputs
    !    j     integer
    !  Outputs
    !    ix    list of all indices of wavelet coeffts at j-th level
    ! Adpoted from wavelab 850
    ! i = (2^(j)+1):(2^(j+1))
    ! zero-based indexing
    i1 = 2**J
    i2 = 2**(J+1)-1
  end subroutine dyad


  subroutine reverse( x )
    implicit none
    real(kind=pr), dimension(1:), intent(inout) :: x
    real(kind=pr), dimension(:), allocatable :: tmp
    integer :: i, N

    N = size(x)
    allocate( tmp(1:N) )
    tmp = x

    do i=1,N
      x(i) = tmp(N-(i-1))
    enddo

    deallocate(tmp)
  end subroutine

  ! out of place hard thresholding
  subroutine hard_thresholding( wc_in, wc_out, threshold)
    implicit none
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc_in
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc_out
    real(kind=pr), intent(in) :: threshold

    wc_out = wc_in

    where ( abs(wc_out) < threshold)
      wc_out = 0.d0
    end where
  end subroutine

  ! count zeros in field
  subroutine compression_rate( wc, rate )
    implicit none
    real(kind=pr), dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)), intent(inout) :: wc
    real(kind=pr), intent(out) :: rate
    integer :: i, N, ix, iy, iz

    i = 0
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          if (abs(wc(ix,iy,iz))<1e-10) then
            i = i +1
          endif
        enddo
      enddo
    enddo
    rate = dble(i)
    rate = mpisum(rate)
    rate = rate / ( dble(nx)*dble(ny)*dble(nz) )
  end subroutine

end module
