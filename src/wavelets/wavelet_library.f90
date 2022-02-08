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

    ! --------------------------------------------------------------------------

    ! Haar wavelet
    ! wavelet%name = "haar"
    ! ! ****** decomposition filter ******
    ! allocate( wavelet%HD(-1:1) ) ! TILDE
    ! wavelet%HD = (/0.0d0, 1.0d0/dsqrt(2.0d0), 1.0d0/dsqrt(2.0d0)/) ! TILDE
    ! allocate( wavelet%GD(-1:1) ) ! TILDE
    ! wavelet%GD = (/ 0.0d0, 1.0d0, -1.0d0/) / sqrt(2.0d0) ! TILDE
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-1:1) )
    ! wavelet%HR = wavelet%HD
    ! allocate( wavelet%GR(-1:1) )
    ! wavelet%GR = wavelet%GD

    ! --------------------------------------------------------------------------

    ! ! 2nd order interpolating wavelet, no lifted (as used in wabbit)
    ! ! with 2nd order multiresolution
    ! wavelet%name = "CDF2,0"
    ! ! ****** decomposition filter ****** ! TILDE
    ! allocate( wavelet%HD(-1:1) ) ! TILDE
    ! wavelet%HD = (/0.0d0, 1.0d0, 0.0d0/) ! TILDE
    ! allocate( wavelet%GD(0:2) ) ! TILDE
    ! wavelet%GD = (/ -0.5d0, 1.0d0, -0.5d0/) ! TILDE
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-1:1) )
    ! wavelet%HR = (/0.5d0, 1.0d0, 0.5d0/)
    ! allocate( wavelet%GR(-1:1) )
    ! wavelet%GR = (/0.0d0, 0.0d0, 1.0d0/)

    ! --------------------------------------------------------------------------

    ! ! fourth order interpolating wavelet, as used in wabbit
    ! ! deslauriers-dubuc 4th order.
    ! wavelet%name = "CDF4,0"
    ! ! ****** decomposition filter ******
    ! allocate( wavelet%HD(-1:1) )  ! H TILDE
    ! wavelet%HD = (/0.0d0, 1.0d0, 0.0d0/) ! TILDE
    ! allocate( wavelet%GD(-2:4) ) ! G TILDE
    ! wavelet%GD = (/+1.0d0/16.0d0, 0.0d0, -9.0d0/16.0d0, 1.0d0, -9.0d0/16.0d0, 0.0d0, +1.0d0/16.0d0 /) ! TILDE
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-3:3) )
    ! wavelet%HR = (/-1.0d0/16.0d0, 0.0d0, +9.0d0/16.0d0, 1.0d0, +9.0d0/16.0d0, 0.0d0, -1.0d0/16.0d0/) ! H
    ! allocate( wavelet%GR(-1:1) )
    ! wavelet%GR = (/0.0d0, 0.0d0, 1.0d0/) ! G
    ! --------------------------------------------------------------------------

    ! fourth order interpolating wavelet, as used in wabbit
    ! deslauriers-dubuc 4th order.
    ! wavelet%name = "CDF4,2"
    ! ! ****** decomposition filter ******
    ! allocate( wavelet%HD(-4:4) )  ! H TILDE
    ! wavelet%HD = (/2.0d0**-6, 0.0d0, -2.0d0**-3, 2.0d0**-2, 23.0d0*2.0d0**-5, 2.0d0**-2, -2.0d0**-3, 0.0d0, 2.0d0**-6/) ! TILDE
    ! allocate( wavelet%GD(-2:4) ) ! G TILDE
    ! wavelet%GD = (/+1.0d0/16.0d0, 0.0d0, -9.0d0/16.0d0, 1.0d0, -9.0d0/16.0d0, 0.0d0, +1.0d0/16.0d0 /) ! TILDE
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-3:3) )
    ! wavelet%HR = (/-1.0d0/16.0d0, 0.0d0, +9.0d0/16.0d0, 1.0d0, +9.0d0/16.0d0, 0.0d0, -1.0d0/16.0d0/) ! H
    ! allocate( wavelet%GR(-5:3) )
    ! wavelet%GR = (/2.0d0**-6, -0.0d0, -2.0d0**-3, -2.0d0**-2, 23.0d0*2.0d0**-5, -2.0d0**-2, -2.0d0**-3, -0.0d0, 2.0d0**-6/)  !G

    wavelet%name = "CDF4,4"
    ! ****** decomposition filter ******
    allocate( wavelet%HD(-6:6) )  ! H TILDE
    wavelet%HD = (/ -2.0d0**-9, &
                     0.0d0, &
                     9.0d0*2.0d0**-8, &
                     -2.0d0**-5, &
                     -63.0d0*2.0d0**-9, &
                     9.0d0*2.0d0**-5, &
                     87.0d0*2.0d0**-7 , &
                     9.0d0*2.0d0**-5, &
                     -63.0d0*2.0d0**-9, -2.0d0**-5, 9.0d0*2.0d0**-8, 0.0d0, -2.0d0**-9/) ! TILDE
    allocate( wavelet%GD(-2:4) ) ! G TILDE
    wavelet%GD = (/+1.0d0/16.0d0, 0.0d0, -9.0d0/16.0d0, 1.0d0, -9.0d0/16.0d0, 0.0d0, +1.0d0/16.0d0 /) ! TILDE
    ! ****** reconstruction filter ******
    allocate( wavelet%HR(-3:3) )
    wavelet%HR = (/-1.0d0/16.0d0, 0.0d0, +9.0d0/16.0d0, 1.0d0, +9.0d0/16.0d0, 0.0d0, -1.0d0/16.0d0/) ! H
    allocate( wavelet%GR(-7:5) )
    wavelet%GR = (/ -2.0d0**-9, -0.0d0, 9.0d0*2.0d0**-8, +2.0d0**-5, -63.0d0*2.0d0**-9, -9.0d0*2.0d0**-5, 87.0d0*2.0d0**-7 , &
    -9.0d0*2.0d0**-5, -63.0d0*2.0d0**-9, +2.0d0**-5, 9.0d0*2.0d0**-8, -0.0d0, -2.0d0**-9/)  !G


    ! ! --------------------------------------------------------------------------
    ! !
    ! ! lifted 2nd order interpolating wavelet, no lifted (as used in wabbit)
    ! wavelet%name = "CDF2,0+2"
    ! ! ****** decomposition filter ****** ! TILDE
    ! allocate( wavelet%HD(-2:2) )  ! H TILDE
    ! wavelet%HD =  (-1.0d0)*(/+1.0d0/8.0d0, -1.0d0/4.0d0, -3.0d0/4.0d0, -1.0d0/4.0d0, +1.0d0/8.0d0/) ! H TILDE
    ! allocate( wavelet%GD(0:2) ) ! TILDE
    ! wavelet%GD = (/ -0.5d0, 1.0d0, -0.5d0/) ! TILDE
    !
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-1:1) )
    ! wavelet%HR = (/0.5d0, 1.0d0, 0.5d0/)
    ! allocate( wavelet%GR(-3:1) )
    ! wavelet%GR = (/-1.0d0/8.0d0, -1.0d0/4.0d0, +3.0d0/4.0d0, -1.0d0/4.0d0, -1.0d0/8.0d0/) ! G


    ! ! --------------------------------------------------------------------------
    ! !
    ! ! coiflet
    ! wavelet%name = "coiflet12"
    ! ! ****** decomposition filter ****** ! TILDE
    ! ! low-pass filter H is in Wavelab850:
    ! !       qmf = MakeONFilter('Coiflet',2)
    ! ! its sum is =sqrt(2)
    ! ! ATTENTION: the coefficients in wikipedia https://en.wikipedia.org/wiki/Coiflet ,
    ! ! here called C12, are normalized to =2, so multiply them with sqrt(2)/2
    ! allocate( wavelet%HD(-4:7) )
    ! wavelet%HD = (/1.638733646318000d-02, -4.146493678197000d-02,-6.737255472230000d-02,&
    ! 3.861100668230900d-01,8.127236354496100d-01,4.170051844237800d-01,-7.648859907826000d-02,&
    ! -5.943441864647000d-02,2.368017194688000d-02,5.611434819370000d-03,-1.823208870910000d-03,&
    ! -7.205494453700000d-04/)
    !
    ! call mirror_filter( wavelet%HD, lbound(wavelet%HD,dim=1), ubound(wavelet%HD,dim=1), wavelet%GD)
    !
    ! ! allocate( wavelet%GD(-6:5) ) ! TILDE
    ! ! wavelet%GD = (/ 7.205494453700000d-04,-1.823208870910000d-03,-5.611434819370000d-03,&
    ! ! 2.368017194688000d-02,5.943441864647000d-02,-7.648859907826000d-02,-4.170051844237800d-01,&
    ! ! 8.127236354496100d-01,-3.861100668230900d-01,-6.737255472230000d-02,4.146493678197000d-02,&
    ! ! 1.638733646318000d-02/) ! TILDE
    !
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-7:4) )
    ! wavelet%HR = wavelet%HD
    !
    ! allocate( wavelet%GR(-5:6) )
    ! wavelet%GR = wavelet%GD

    ! ! --------------------------------------------------------------------------
    ! !
    ! ! coiflet using coeffs from http://wavelets.pybytes.com/wavelet/coif2/
    ! ! wavelet%name = "coiflet12"
    ! ! ****** decomposition filter ****** ! TILDE
    ! allocate( wavelet%HD(0:11) )
    ! wavelet%HD = (/-0.0007205494453645122d0,&
    ! -0.0018232088707029932d0,&
    ! 0.0056114348193944995d0,&
    ! 0.023680171946334084d0,&
    ! -0.0594344186464569d0,&
    ! -0.0764885990783064d0,&
    ! 0.41700518442169254d0,&
    ! 0.8127236354455423d0,&
    ! 0.3861100668211622d0,&
    ! -0.06737255472196302d0,&
    ! -0.04146493678175915d0,&
    ! 0.016387336463522112d0&
    ! /)
    !
    ! allocate( wavelet%GD(0:11) ) ! TILDE
    ! wavelet%GD = (/-0.016387336463522112d0,&
    ! -0.04146493678175915d0,&
    ! 0.06737255472196302d0,&
    ! 0.3861100668211622d0,&
    ! -0.8127236354455423d0,&
    ! 0.41700518442169254d0,&
    ! 0.0764885990783064d0,&
    ! -0.0594344186464569d0,&
    ! -0.023680171946334084d0,&
    ! 0.0056114348193944995d0,&
    ! 0.0018232088707029932d0,&
    ! -0.0007205494453645122d0&
    ! /) ! TILDE
    !
    ! ! ****** reconstruction filter ******
    ! allocate( wavelet%HR(-11:0) )
    ! wavelet%HR = wavelet%HD
    !
    ! allocate( wavelet%GR(-11:0) )
    ! wavelet%GR = wavelet%GD

    ! --------------------------------------------------------------------------

    call reverse( wavelet%HR )
    call reverse( wavelet%GR )

    if (root) then
        write(*,*) "-------------------wavelet----------------------"
        write(*,*) "HD (H_tilde) bounds:", LBOUND(wavelet%HD,dim=1), UBOUND(wavelet%HD,dim=1)
        write(*,*) "sum(HD)=", sum(wavelet%HD)
        write(*,'(100(es12.4,1x))') wavelet%HD
        write(*,*) "------------------------------------------------"
        write(*,*) "GD (G_tilde) bounds:", LBOUND(wavelet%GD,dim=1), UBOUND(wavelet%GD,dim=1)
        write(*,*) "sum(GD)=", sum(wavelet%GD)
        write(*,'(100(es12.4,1x))') wavelet%GD
        write(*,*) "------------------------------------------------"
        write(*,*) "HR (H) bounds:", LBOUND(wavelet%HR,dim=1), UBOUND(wavelet%HR,dim=1)
        write(*,*) "sum(HR)=", sum(wavelet%HR)
        write(*,'(100(es12.4,1x))') wavelet%HR
        write(*,*) "------------------------------------------------"
        write(*,*) "GR (G) bounds:", LBOUND(wavelet%GR,dim=1), UBOUND(wavelet%GR,dim=1)
        write(*,*) "sum(GR)=", sum(wavelet%GR)
        write(*,'(100(es12.4,1x))') wavelet%GR
        write(*,*) "------------------------------------------------"
    endif

  end subroutine

  ! given a g filter, compute the h filter.
  subroutine mirror_filter(g, a, b, h)
      implicit none
      integer, intent(in) :: a,b
      real(kind=pr), intent(in) :: g(a:b)
      real(kind=pr), intent(out), allocatable :: h(:)

      integer :: i

      allocate( h(-b+1:-a+1) )

      ! wavelab850matlab
	  ! y = -( (-1).^(1:length(x)) ).*x;
      ! so wavelab really flips every second sign, regardless of what the lower
      ! bound index is (which is thus supposed ZERO)

      h = 0.0d0

      ! it is still unclear to me which coeffs do need sign inversion.
      ! for coiflet, it seems to be the even ones. wikipedia says "every second one"
      do i = -b+1, -a+1
          if ((1-i)>=a .and. (1-i<=b)) then
              h(i) = (-1.0d0)**(1.0d0-dble(i)) * g(1-i)
          endif
      enddo

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
