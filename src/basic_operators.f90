!-------------------------------------------------------------------------------
! this module contains basic elementary operators, and as it is a module
! function overloading can be used.
! List of functions:
!       * curl
!       * divergence
!       * max/min operations on real fields
!       * check fields for NaN
!-------------------------------------------------------------------------------
module basic_operators
  !-- interface for curl operators
  interface curl
    module procedure curl, curl_inplace, curl3, curl3_inplace
  end interface

  !-- interface for maxima of fields
  interface field_max_magnitude
    module procedure field_max_magnitude, field_max_magnitude3
  end interface

  !-- check fields for NaN
  interface checknan
    module procedure checknan_cmplx, checknan_real
  end interface

  interface dealias
    module procedure dealias, dealias1, dealias3
  end interface

!!!!!!!!!!!
 contains
!!!!!!!!!!!


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl(out1,out2,out3,in1,in2,in3)
  use p3dfft_wrapper
  use vars
  implicit none

  ! input field in Fourier space
  complex(kind=pr),intent(in)::in1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in)::in2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in)::in3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  ! output field in Fourier space
  complex(kind=pr),intent(out)::out1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(out)::out2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(out)::out3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  ! Compute curl of given field in Fourier space:
  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)

           out1(iz,iy,ix)=imag*(ky*in3(iz,iy,ix) -kz*in2(iz,iy,ix))
           out2(iz,iy,ix)=imag*(kz*in1(iz,iy,ix) -kx*in3(iz,iy,ix))
           out3(iz,iy,ix)=imag*(kx*in2(iz,iy,ix) -ky*in1(iz,iy,ix))
        enddo
     enddo
  enddo
end subroutine curl


! Given three components of a fields in Fourier space, compute the
! curl in physical space.  Arrays are 3-dimensional.
subroutine curl_inplace(fx,fy,fz)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  ! Field in Fourier space
  complex(kind=pr),intent(inout)::fx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::fy(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::fz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  complex(kind=pr) :: t1,t2,t3 ! temporary loop variables
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  ! Compute curl of given field in Fourier space:
  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)

           t1=fx(iz,iy,ix)
           t2=fy(iz,iy,ix)
           t3=fz(iz,iy,ix)

           fx(iz,iy,ix)=imag*(ky*t3 - kz*t2)
           fy(iz,iy,ix)=imag*(kz*t1 - kx*t3)
           fz(iz,iy,ix)=imag*(kx*t2 - ky*t1)
        enddo
     enddo
  enddo
end subroutine curl_inplace


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl3_inplace(fk)
  use vars
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(inout)::fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  call curl_inplace( fk(:,:,:,1), fk(:,:,:,2), fk(:,:,:,3) )
end subroutine curl3_inplace


! Given three components of an input fields in Fourier space,
! the curl in physical space. The precision is reduced
! to second order in space, since this kind of filtering may help in postprocessing
subroutine curl3_2nd_inplace(fk)
  use vars
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(inout)::fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  call curl_2nd( fk(:,:,:,1), fk(:,:,:,2), fk(:,:,:,3) )
end subroutine curl3_2nd_inplace


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional. The precision is reduced
! to second order in space, since this kind of filtering may help in postprocessing
subroutine curl_2nd (fx,fy,fz)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  ! Field in Fourier space
  complex(kind=pr),intent(inout)::fx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::fy(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::fz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  complex(kind=pr) :: t1,t2,t3 ! temporary loop variables
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  ! Compute curl of given field in Fourier space:
  do ix=ca(3),cb(3)
     kx=dsin(dx*wave_x(ix))/dx ! (reduced to 2nd order)
     do iy=ca(2),cb(2)
        ky=dsin(dy*wave_y(iy))/dy ! (reduced to 2nd order)
        do iz=ca(1),cb(1)
           kz=dsin(dz*wave_z(iz))/dz ! (reduced to 2nd order)

           t1=fx(iz,iy,ix)
           t2=fy(iz,iy,ix)
           t3=fz(iz,iy,ix)

           fx(iz,iy,ix)=imag*(ky*t3 - kz*t2)
           fy(iz,iy,ix)=imag*(kz*t1 - kx*t3)
           fz(iz,iy,ix)=imag*(kx*t2 - ky*t1)
        enddo
     enddo
  enddo

end subroutine curl_2nd


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl3(ink,outk)
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

   ! Compute curl of given field in Fourier space:
  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)

           outk(iz,iy,ix,1)=imag*(ky*ink(iz,iy,ix,3) -kz*ink(iz,iy,ix,2))
           outk(iz,iy,ix,2)=imag*(kz*ink(iz,iy,ix,1) -kx*ink(iz,iy,ix,3))
           outk(iz,iy,ix,3)=imag*(kx*ink(iz,iy,ix,2) -ky*ink(iz,iy,ix,1))
        enddo
     enddo
  enddo
end subroutine curl3


! compute the curl of vector field u with PERIODIC Finite-differences of
! 2nd or 4th order. no Fourier transform is used, the curl here is very fast
! but may be ess precise than fourier.
! Ghost nodes are used for intra-processor communication
subroutine curl_FD( u, rotu, order )
  use p3dfft_wrapper
  use ghosts
  implicit none

  ! input/output field in x-space
  ! NOTE ghost nodes are used!
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout)::rotu(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  character(len=*), intent(in) :: order

  integer :: ix,iy,iz
  real(kind=pr) :: dxinv,dyinv,dzinv,a1,a2,a4,a5
  real(kind=pr) :: uxdy,uxdz,uydx,uydz,uzdx,uzdy

  select case(order)
  case('centered_2nd')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute curl(u) using second order centered period FD
      !-------------------------------------------------------------------------
      dxinv = 1.d0/(2.d0*dx)
      dyinv = 1.d0/(2.d0*dy)
      dzinv = 1.d0/(2.d0*dz)

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dyinv
              uxdz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dzinv
              uydx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dxinv
              uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
              uzdx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dxinv
              uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv

              rotu(ix,iy,iz,1) = uzdy - uydz
              rotu(ix,iy,iz,2) = uxdz - uzdx
              rotu(ix,iy,iz,3) = uydx - uxdy
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
            uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv

            rotu(ix,iy,iz,1) = uzdy - uydz
            rotu(ix,iy,iz,2) = 0.d0
            rotu(ix,iy,iz,3) = 0.d0
          enddo
        enddo
      endif

  case('centered_4th')
      call synchronize_ghosts( u, 3 )
      !-------------------------------------------------------------------------
      ! compute curl(u) using second order centered period FD
      !-------------------------------------------------------------------------
      a1 = 1.d0/12.d0
      a2 =-2.d0/3.d0
      a4 = 2.d0/3.d0
      a5 =-1.d0/12.d0

      dxinv = 1.d0/dx
      dyinv = 1.d0/dy
      dzinv = 1.d0/dz

      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdy = (a1*u(ix,iy-2,iz,1)+a2*u(ix,iy-1,iz,1)+a4*u(ix,iy+1,iz,1)+a5*u(ix,iy+2,iz,1))*dyinv
              uxdz = (a1*u(ix,iy,iz-2,1)+a2*u(ix,iy,iz-1,1)+a4*u(ix,iy,iz+1,1)+a5*u(ix,iy,iz+2,1))*dzinv
              uydx = (a1*u(ix-2,iy,iz,2)+a2*u(ix-1,iy,iz,2)+a4*u(ix+1,iy,iz,2)+a5*u(ix+2,iy,iz,2))*dxinv
              uydz = (a1*u(ix,iy,iz-2,2)+a2*u(ix,iy,iz-1,2)+a4*u(ix,iy,iz+1,2)+a5*u(ix,iy,iz+2,2))*dzinv
              uzdx = (a1*u(ix-2,iy,iz,3)+a2*u(ix-1,iy,iz,3)+a4*u(ix+1,iy,iz,3)+a5*u(ix+2,iy,iz,3))*dxinv
              uzdy = (a1*u(ix,iy-2,iz,3)+a2*u(ix,iy-1,iz,3)+a4*u(ix,iy+1,iz,3)+a5*u(ix,iy+2,iz,3))*dyinv

              rotu(ix,iy,iz,1) = uzdy - uydz
              rotu(ix,iy,iz,2) = uxdz - uzdx
              rotu(ix,iy,iz,3) = uydx - uxdy
            enddo
          enddo
        enddo
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydz = (a1*u(ix,iy,iz-2,2)+a2*u(ix,iy,iz-1,2)+a4*u(ix,iy,iz+1,2)+a5*u(ix,iy,iz+2,2))*dzinv
            uzdy = (a1*u(ix,iy-2,iz,3)+a2*u(ix,iy-1,iz,3)+a4*u(ix,iy+1,iz,3)+a5*u(ix,iy+2,iz,3))*dyinv

            rotu(ix,iy,iz,1) = uzdy - uydz
            rotu(ix,iy,iz,2) = 0.d0
            rotu(ix,iy,iz,3) = 0.d0
          enddo
        enddo
      endif

  case default
      call abort(5556,'invalid METHOD in curl_FD:'//order)
  end select
end subroutine curl_FD



! Compute the curl of vector field u in physical space, using a u
! discretization and Finite-differences. NOTE NO ghost nodes are used. The data
! is transposed and the derivatives are always computed on the contiguous index
! This is the same as we have in the wavelet transform.
subroutine curl_FD_nonper( u, rotu, order )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space NOTE: no ghost nodes
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::rotu(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  character(len=*), intent(in) :: order

  real(kind=pr), allocatable, dimension(:,:,:) :: uzdy, uydz, uxdz, uzdx, uydx, uxdy

  allocate(uzdy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uydz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uxdz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uzdx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uydx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uxdy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))


  call diffy_nonper(u(:,:,:,3),uzdy)
  call diffz_nonper(u(:,:,:,2),uydz)
  call diffz_nonper(u(:,:,:,1),uxdz)
  call diffx_nonper(u(:,:,:,3),uzdx)
  call diffx_nonper(u(:,:,:,2),uydx)
  call diffy_nonper(u(:,:,:,1),uxdy)

  rotu(:,:,:,1) = uzdy - uydz
  rotu(:,:,:,2) = uxdz - uzdx
  rotu(:,:,:,3) = uydx - uxdy

  deallocate(uzdy, uydz, uxdz, uzdx, uydx, uxdy)
end subroutine curl_FD_nonper


! compute the curl of vector field u with PERIODIC Finite-differences of
! 2nd or 4th order. no Fourier transform is used, the curl here is very fast
! but may be ess precise than fourier.
! Ghost nodes are used for intra-processor communication
subroutine Q_FD( u, Q, order )
  use p3dfft_wrapper
  use ghosts
  implicit none

  ! input/output field in x-space
  ! NOTE ghost nodes are used!
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout)::Q(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  character(len=*), intent(in) :: order

  integer :: ix,iy,iz
  real(kind=pr) :: dxinv,dyinv,dzinv,a1,a2,a4,a5
  real(kind=pr) :: uxdx,uxdy,uxdz,uydx,uydy,uydz,uzdx,uzdy,uzdz
  real(kind=pr) :: A(1:3,1:3)

  select case(order)
  case('centered_2nd')
    call synchronize_ghosts( u, 3 )
    !-------------------------------------------------------------------------
    ! compute Q(u) using second order centered period FD
    !-------------------------------------------------------------------------
    dxinv = 1.d0/(2.d0*dx)
    dyinv = 1.d0/(2.d0*dy)
    dzinv = 1.d0/(2.d0*dz)

    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          uxdx = (u(ix+1,iy,iz,1) - u(ix-1,iy,iz,1))*dxinv
          uxdy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dyinv
          uxdz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dzinv

          uydx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dxinv
          uydy = (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dyinv
          uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv

          uzdx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dxinv
          uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv
          uzdz = (u(ix,iy,iz+1,3) - u(ix,iy,iz-1,3))*dzinv

          A(1,:) =  (/ uxdx, uxdy, uxdz/)
          A(2,:) =  (/ uydx, uydy, uydz/)
          A(3,:) =  (/ uzdx, uzdy, uzdz/)

          Q(ix,iy,iz) = -0.5d0*( sum( A*transpose(A) ) )
        enddo
      enddo
    enddo

  case('centered_4th')
    call synchronize_ghosts( u, 3 )
    !-------------------------------------------------------------------------
    ! compute Q(u) using second order centered period FD
    !-------------------------------------------------------------------------
    a1 = 1.d0/12.d0
    a2 =-2.d0/3.d0
    a4 = 2.d0/3.d0
    a5 =-1.d0/12.d0

    dxinv = 1.d0/dx
    dyinv = 1.d0/dy
    dzinv = 1.d0/dz

    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          uxdx = (a1*u(ix-2,iy,iz,1)+a2*u(ix-1,iy,iz,1)+a4*u(ix+1,iy,iz,1)+a5*u(ix+2,iy,iz,1))*dxinv
          uxdy = (a1*u(ix,iy-2,iz,1)+a2*u(ix,iy-1,iz,1)+a4*u(ix,iy+1,iz,1)+a5*u(ix,iy+2,iz,1))*dyinv
          uxdz = (a1*u(ix,iy,iz-2,1)+a2*u(ix,iy,iz-1,1)+a4*u(ix,iy,iz+1,1)+a5*u(ix,iy,iz+2,1))*dzinv

          uydx = (a1*u(ix-2,iy,iz,2)+a2*u(ix-1,iy,iz,2)+a4*u(ix+1,iy,iz,2)+a5*u(ix+2,iy,iz,2))*dxinv
          uydy = (a1*u(ix,iy-2,iz,2)+a2*u(ix,iy-1,iz,2)+a4*u(ix,iy+1,iz,2)+a5*u(ix,iy+2,iz,2))*dyinv
          uydz = (a1*u(ix,iy,iz-2,2)+a2*u(ix,iy,iz-1,2)+a4*u(ix,iy,iz+1,2)+a5*u(ix,iy,iz+2,2))*dzinv

          uzdx = (a1*u(ix-2,iy,iz,3)+a2*u(ix-1,iy,iz,3)+a4*u(ix+1,iy,iz,3)+a5*u(ix+2,iy,iz,3))*dxinv
          uzdy = (a1*u(ix,iy-2,iz,3)+a2*u(ix,iy-1,iz,3)+a4*u(ix,iy+1,iz,3)+a5*u(ix,iy+2,iz,3))*dyinv
          uzdz = (a1*u(ix,iy,iz-2,3)+a2*u(ix,iy,iz-1,3)+a4*u(ix,iy,iz+1,3)+a5*u(ix,iy,iz+2,3))*dzinv

          A(1,:) =  (/ uxdx, uxdy, uxdz/)
          A(2,:) =  (/ uydx, uydy, uydz/)
          A(3,:) =  (/ uzdx, uzdy, uzdz/)

          Q(ix,iy,iz) = -0.5d0*( sum( A*transpose(A) ) )
        enddo
      enddo
    enddo

  case default
      call abort(5556,'invalid METHOD in curl_FD:'//order)
  end select
end subroutine Q_FD



! Compute the Q-criterion of vector field u in physical space, using a u
! discretization and Finite-differences. NOTE NO ghost nodes are used. The data
! is transposed and the derivatives are always computed on the contiguous index
! This is the same as we have in the wavelet transform.
subroutine Q_FD_nonper( u, Q, order )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space NOTE: no ghost nodes
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::Q(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*), intent(in) :: order

  real(kind=pr), allocatable, dimension(:,:,:) :: uxdx, uxdy, uxdz
  real(kind=pr), allocatable, dimension(:,:,:) :: uydx, uydy, uydz
  real(kind=pr), allocatable, dimension(:,:,:) :: uzdx, uzdy, uzdz
  real(kind=pr) :: A(1:3,1:3)
  integer :: ix, iy, iz

  allocate(uxdx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uxdy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uxdz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uydx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uydy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uydz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uzdx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uzdy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uzdz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call diffx_nonper(u(:,:,:,1),uxdx)
  call diffy_nonper(u(:,:,:,1),uxdy)
  call diffz_nonper(u(:,:,:,1),uxdz)

  call diffx_nonper(u(:,:,:,2),uydx)
  call diffy_nonper(u(:,:,:,2),uydy)
  call diffz_nonper(u(:,:,:,2),uydz)

  call diffx_nonper(u(:,:,:,3),uzdx)
  call diffy_nonper(u(:,:,:,3),uzdy)
  call diffz_nonper(u(:,:,:,3),uzdz)

  ! source: https://github.com/ganglere/matlab/blob/master/VortexID.m

  do iz = ra(3),rb(3)
    do iy = ra(2),rb(2)
      do ix = ra(1),rb(1)

        A(1,:) =  (/ uxdx(ix,iy,iz), uxdy(ix,iy,iz), uxdz(ix,iy,iz)/)
        A(2,:) =  (/ uydx(ix,iy,iz), uydy(ix,iy,iz), uydz(ix,iy,iz)/)
        A(3,:) =  (/ uzdx(ix,iy,iz), uzdy(ix,iy,iz), uzdz(ix,iy,iz)/)

        Q(ix,iy,iz) = -0.5d0*( sum( A*transpose(A) ) )
      enddo
    enddo
  enddo

  deallocate( uxdx, uxdy, uxdz)
  deallocate( uydx, uydy, uydz)
  deallocate( uzdx, uzdy, uzdz)
end subroutine Q_FD_nonper


! Compute derivative in x-direction of a scalar field using a non-periodic
! discretization and Finite-differences. NOTE NO ghost nodes are used. The data
! is transposed and the derivatives are always computed on the contiguous index
! This is the same as we have in the wavelet transform.
! As the x-index is always contiguous, this is the easiest case and no MPI communication
! is involved.
subroutine diffx_nonper( u, udx )
  use vars
  use p3dfft_wrapper
  implicit none
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::udx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz
  real(kind=pr), dimension(:,:), allocatable :: D1

  allocate(D1(0:nx-1,0:nx-1))
  call deriv_matrix_1d_nonper(D1, dx)

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      udx(:,iy,iz)= matmul(D1,u(:,iy,iz))
    enddo
  enddo

  deallocate(D1)

end subroutine



! Compute derivative in y-direction of a scalar field using a non-periodic
! discretization and Finite-differences. NOTE NO ghost nodes are used. The data
! is transposed and the derivatives are always computed on the contiguous index
! This is the same as we have in the wavelet transform.
! NOTE the y-direction is sometimes contiguous, sometimes not.
subroutine diffy_nonper( u, udx )
  use vars
  use p3dfft_wrapper
  implicit none
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::udx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz,idir
  real(kind=pr), dimension(:,:), allocatable :: D1
  real(kind=pr), allocatable, dimension(:,:,:) :: work
  integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst

  allocate(D1(0:ny-1,0:ny-1))
  call deriv_matrix_1d_nonper(D1, dy)


  if (mpidims(2) == 1) then
    ! CASE A): index 2 is contiguous and not split among procs.
    do iz=ra(3),rb(3)
      do ix=ra(1),rb(1)
        udx(ix,:,iz)= matmul(D1,u(ix,:,iz))
      enddo
    enddo

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
    do iz=kat(3),kbt(3)
      do ix=kat(2),kbt(2)
        work(:,ix,iz)= matmul(D1,work(:,ix,iz))
      enddo
    enddo

    ! transposition: exchange x-y data
    ! idir = 1 - permute 1 and 2 indices
    ! idir = 2 - permute 1 and 3 indices
    idir = 1
    call trextents(idir, (/ny,nx,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
    call subtr ( idir, (/ny,nx,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
                 MPI_DOUBLE_PRECISION, work, udx )
    ! at this point we have again udx(ix,iy,iz), which is what we started with
  endif
  if(allocated(work)) deallocate(work)
  deallocate(D1)
end subroutine





! Compute derivative in z-direction of a scalar field using a non-periodic
! discretization and Finite-differences. NOTE NO ghost nodes are used. The data
! is transposed and the derivatives are always computed on the contiguous index
! This is the same as we have in the wavelet transform.
! NOTE the z-direction is never contiguous, if more than 1CPU is used.
subroutine diffz_nonper( u, udx )
  use vars
  use p3dfft_wrapper
  implicit none

  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::udx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz,idir
  real(kind=pr), dimension(:,:), allocatable :: D1
  real(kind=pr), allocatable, dimension(:,:,:) :: work
  integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst

  allocate(D1(0:nz-1,0:nz-1))
  call deriv_matrix_1d_nonper(D1, dz)

  ! transposition: exchange x-y data
  ! idir = 1 - permute 1 and 2 indices
  ! idir = 2 - permute 1 and 3 indices
  idir = 2
  call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  allocate( work(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
  call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
               MPI_DOUBLE_PRECISION, u, work )

  ! --- z direction -----
  ! at this point we have work(iz,iy,ix)
  do iy=kat(2),kbt(2)
    do ix=kat(3),kbt(3)
      work(:,iy,ix)= matmul(D1,work(:,iy,ix))
    enddo
  enddo

  ! final permutation back to the original dataset
  idir = 2
  call trextents(idir, (/nz,ny,nx/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  call subtr ( idir, (/nz,ny,nx/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, &
               MPI_DOUBLE_PRECISION, work, udx )
  if(allocated(work)) deallocate(work)
  ! at this point we have wc(ix,iy,iz)
  deallocate(D1)
end subroutine


! create a matrix to compute the first derivative, non-periodic
subroutine deriv_matrix_1d_nonper(D1, h)
  use vars, only : pr
  implicit none
  real(kind=pr), intent(inout) :: D1(0:,0:), h
  integer :: nn,N,i

  nn = size(D1,1)
  D1=0.d0
  do i=0,nn-1 !main diag
    D1(i,i) = 0.d0 !main
  enddo
  do i=0,nn-2 !first upper/lower diag
    D1(i,i+1) = 1.d0  !up1
    D1(i+1,i) =-1.d0  !lo1
  enddo
  N = nn-1
  D1(0,0)=  -3.d0
  D1(0,1)=   4.d0
  D1(0,2)=  -1.d0
  D1(N,N-2)= 1.d0
  D1(N,N-1)=-4.d0
  D1(N,N)=   3.d0
  D1 = D1/(2.0d0*h)

end subroutine


!-------------------------------------------------------------------------------
! compute divergence of vector valued field (ink, 4D array)
! and returns it in outk
!-------------------------------------------------------------------------------
subroutine divergence( ink, outk )
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none
  ! input vector field in Fourier space
  complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  ! output scalar field in Fourier space
  complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz

  do ix=ca(3),cb(3)
    !-- wavenumber in x-direction
    kx = wave_x(ix)
    do iy=ca(2),cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do iz=ca(1),cb(1)
        !-- wavenumber in z-direction
        kz = wave_z(iz)
        outk(iz,iy,ix)=kx*ink(iz,iy,ix,1)+ky*ink(iz,iy,ix,2)+kz*ink(iz,iy,ix,3)
        outk(iz,iy,ix)=dcmplx(0.d0,1.d0)*outk(iz,iy,ix)
      enddo
    enddo
  enddo
end subroutine divergence


! computes laplace(ink) for a scalar valued field and returns it in the same array
subroutine laplacien_inplace( ink )
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none
  complex(kind=pr),intent(inout)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2

  do ix=ca(3),cb(3)
    !-- wavenumber in x-direction
    kx = wave_x(ix)
    do iy=ca(2),cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do iz=ca(1),cb(1)
        !-- wavenumber in z-direction
        kz = wave_z(iz)
        k2 = kx*kx + ky*ky + kz*kz
        ink(iz,iy,ix) = -k2*ink(iz,iy,ix)
      enddo
    enddo
  enddo
end subroutine laplacien_inplace


! computes laplace(ink) for a scalar valued field and returns it in the same array
! note wavenumbers are reduced to second order accuracy (for the Q-criterion,
! the result is nicer if filtered)
subroutine laplacien_inplace_filtered( ink )
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none
  complex(kind=pr),intent(inout)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2

  do ix=ca(3),cb(3)
    !-- wavenumber in x-direction  (reduced to 2nd order)
    kx = dsin(dx*wave_x(ix))/dx
    do iy=ca(2),cb(2)
      !-- wavenumber in y-direction  (reduced to 2nd order)
      ky = dsin( dy*wave_y(iy) )/dy
      do iz=ca(1),cb(1)
        !-- wavenumber in z-direction (reduced to 2nd order)
        kz = dsin( dz*wave_z(iz) )/ dz
        k2 = kx*kx + ky*ky + kz*kz
        ink(iz,iy,ix) = -k2*ink(iz,iy,ix)
      enddo
    enddo
  enddo
end subroutine laplacien_inplace_filtered


! returns the globally largest entry of a given (real) field
real(kind=pr) function fieldmax( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: max_local, max_global
  integer :: mpicode

  max_local = maxval(inx)
  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmax = max_global
end function fieldmax


! returns the globally smallest entry of a given (real) field
real(kind=pr) function fieldmin( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: min_local, min_global
  integer :: mpicode

  min_local = minval(inx)
  call MPI_ALLREDUCE (min_local,min_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmin = min_global
end function fieldmin


! returns the mean value of a given field
real(kind=pr) function fieldmean( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: mean_local, mean_global
  integer :: mpicode
  mean_local = sum(inx)
  call MPI_ALLREDUCE (mean_local,mean_global,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmean = mean_global / (dble(nx)*dble(ny)*dble(nz))
end function fieldmean


! returns the globally largest entry of a given vector field
! (L2-norm)
real(kind=pr) function field_max_magnitude3( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr) :: max_local, max_global, value
  integer :: mpicode
  integer ::ix,iy,iz

  max_local = 0.d0
  do iz = ra(3),rb(3)
    do iy = ra(2),rb(2)
      do ix = ra(1),rb(1)
        value = inx(ix,iy,iz,1)*inx(ix,iy,iz,1) + inx(ix,iy,iz,2)*inx(ix,iy,iz,2) &
              + inx(ix,iy,iz,3)*inx(ix,iy,iz,3)
        if (max_local<value) max_local=value
       enddo
    enddo
  enddo

  max_local = dsqrt( max_local )
  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  field_max_magnitude3 = max_global
end function field_max_magnitude3


! returns the globally largest entry of a given vector field
! (L2-norm)
real(kind=pr) function field_max_magnitude( inx1, inx2, inx3 )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: inx2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: inx3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: max_local, max_global
  integer :: mpicode

  max_local = maxval( inx1*inx1 + inx2*inx2  + inx3*inx3 )
  max_local = dsqrt( max_local )

  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value
  field_max_magnitude = max_global
end function field_max_magnitude


!-------------------------------------------------------------------------------
! check a real valued field for NaNs and display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_real( field, msg )
  use vars
  implicit none
  real(kind=pr),intent(in)::field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*),intent(in)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz

  foundnan = 0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        if (is_nan(field(ix,iy,iz))) foundnan = 1
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)

  if (foundnans>0) then
    call abort(201,'NaN in '//msg)
  endif
end subroutine checknan_real


!-------------------------------------------------------------------------------
! check a complex field for NaN's, display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_cmplx( field, msg )
  use vars
  implicit none
  complex(kind=pr),intent(in)::field(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  character(len=*),intent(in)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz

  foundnan = 0
  do ix=ca(3),cb(3)
    do iy=ca(2),cb(2)
      do iz=ca(1),cb(1)
        if (is_nan(aimag(field(iz,iy,ix)))) foundnan = 1
        if (is_nan(real (field(iz,iy,ix)))) foundnan = 1
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)

  if (foundnans>0) then
    call abort(202,'NaN in '//msg)
  endif
end subroutine checknan_cmplx
!-------------------------------------------------------------------------------
subroutine dealias(fk1,fk2,fk3)
  use vars
  use penalization ! mask array etc
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)

  do iz=ca(1),cb(1)
    kzt2=(wave_z(iz)/scalez) / kz_trunc
    kzt2=kzt2*kzt2

    do iy=ca(2),cb(2)
      kyt2=(wave_y(iy)/scaley) / ky_trunc
      kyt2=kyt2*kyt2

      do ix=ca(3),cb(3)
        kxt2=(wave_x(ix)/scalex) / kx_trunc
        kxt2=kxt2*kxt2

        if (kxt2 + kyt2 + kzt2  >= 1.d0) then
          fk1(iz,iy,ix)=0.d0
          fk2(iz,iy,ix)=0.d0
          fk3(iz,iy,ix)=0.d0
        endif

      enddo
    enddo
  enddo

end subroutine dealias


!-------------------------------------------------------------------------------
subroutine dealias1(fk1)
  use vars
  use penalization ! mask array etc
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)

  do iz=ca(1),cb(1)
    kzt2=(wave_z(iz)/scalez) / kz_trunc
    kzt2=kzt2*kzt2

    do iy=ca(2),cb(2)
      kyt2=(wave_y(iy)/scaley) / ky_trunc
      kyt2=kyt2*kyt2

      do ix=ca(3),cb(3)
        kxt2=(wave_x(ix)/scalex) / kx_trunc
        kxt2=kxt2*kxt2

        if (kxt2 + kyt2 + kzt2  >= 1.d0) then
          fk1(iz,iy,ix)=0.d0
        endif

      enddo
    enddo
  enddo

end subroutine dealias1


!-------------------------------------------------------------------------------
subroutine dealias3(fk)
  use vars
  use penalization ! mask array etc
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr) :: kxt2,kyt2,kzt2
  real(kind=pr) :: kx_trunc,ky_trunc,kz_trunc

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)

  do iz=ca(1),cb(1)
    kzt2=(wave_z(iz)/scalez) / kz_trunc
    kzt2=kzt2*kzt2

    do iy=ca(2),cb(2)
      kyt2=(wave_y(iy)/scaley) / ky_trunc
      kyt2=kyt2*kyt2

      do ix=ca(3),cb(3)
        kxt2=(wave_x(ix)/scalex) / kx_trunc
        kxt2=kxt2*kxt2

        if (kxt2 + kyt2 + kzt2  >= 1.d0) then
          fk(iz,iy,ix,:)=0.d0
        endif
      enddo
    enddo
  enddo

end subroutine dealias3

!-------------------------------------------------------------------------------
! From a given vorticity field in Fourier space, compute the velocity field
! using Biot-Savarts Law (which is the inverse of the curl operator)
! INPUT:
!   vork: vorticity in Fourier space
! OUTPUT:
!   uk: velocity field in Fourier space
! NOTES:
!   The inverse of the curl is defined up to an irrotational part (eg constant flow)
!   and we thus set the zeroth Fourier mode of uk to zero
!   In older versions, still (07/2015) present in inicond/fsi/init_fields_fsi.f90
!   we include the FFT in the routine. that is not very nice.
!-------------------------------------------------------------------------------
subroutine Vorticity2Velocity(vork, uk)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr),intent(in ) :: vork(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex(kind=pr),intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix, iy, iz, i
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  complex (kind=pr) :: im, spx,spy,spz
  ! imaginary unit
  im = dcmplx(0.d0,1.d0)

  do ix=ca(3),cb(3)
    kx=wave_x(ix)
    kx2=kx*kx
    do iy=ca(2),cb(2)
      ky=wave_y(iy)
      ky2=ky*ky
      do iz=ca(1),cb(1)
        kz=wave_z(iz)
        kz2=kz*kz

        k_abs_2=kx2+ky2+kz2
        if ( k_abs_2 >= 1.0d-13) then
          ! we first "solve" the poisson eqn
          ! which gives us the streamfunction components
          spx = vork(iz,iy,ix,1) / k_abs_2
          spy = vork(iz,iy,ix,2) / k_abs_2
          spz = vork(iz,iy,ix,3) / k_abs_2
          ! compute curl of streamfunction (=velocity)
          uk(iz,iy,ix,1)=im*(ky*spz - kz*spy)
          uk(iz,iy,ix,2)=im*(kz*spx - kx*spz)
          uk(iz,iy,ix,3)=im*(kx*spy - ky*spx)
        else
          ! set zeo mode to zero (no constant flow)
          uk(iz,iy,ix,1) = dcmplx(0.d0,0.d0)
          uk(iz,iy,ix,2) = dcmplx(0.d0,0.d0)
          uk(iz,iy,ix,3) = dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo
end subroutine Vorticity2Velocity


!-------------------------------------------------------------------------------
! From a given vorticity field in Fourier space, compute the velocity field
! using Biot-Savarts Law (which is the inverse of the curl operator)
! INPUT:
!   uk: vorticity in Fourier space
! OUTPUT:
!   uk: velocity field in Fourier space
! NOTES:
!   The inverse of the curl is defined up to and irrotational part (eg constant flow)
!   and we thus set the zeroth Fourier mode of uk to zero
!-------------------------------------------------------------------------------
subroutine Vorticity2Velocity_inplace(uk)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr),intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix, iy, iz, i
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  complex (kind=pr) :: im, spx,spy,spz
  ! imaginary unit
  im = dcmplx(0.d0,1.d0)

  do ix=ca(3),cb(3)
    kx=wave_x(ix)
    kx2=kx*kx
    do iy=ca(2),cb(2)
      ky=wave_y(iy)
      ky2=ky*ky
      do iz=ca(1),cb(1)
        kz=wave_z(iz)
        kz2=kz*kz

        k_abs_2=kx2+ky2+kz2
        if ( k_abs_2 >= 1.0d-13) then
          ! we first "solve" the poisson eqn
          ! which gives us the streamfunction components
          spx = uk(iz,iy,ix,1) / k_abs_2
          spy = uk(iz,iy,ix,2) / k_abs_2
          spz = uk(iz,iy,ix,3) / k_abs_2
          ! compute curl of streamfunction (=velocity)
          uk(iz,iy,ix,1)=im*(ky*spz - kz*spy)
          uk(iz,iy,ix,2)=im*(kz*spx - kx*spz)
          uk(iz,iy,ix,3)=im*(kx*spy - ky*spx)
        else
          ! set zeo mode to zero (no constant flow)
          uk(iz,iy,ix,1) = dcmplx(0.d0,0.d0)
          uk(iz,iy,ix,2) = dcmplx(0.d0,0.d0)
          uk(iz,iy,ix,3) = dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo
end subroutine Vorticity2Velocity_inplace



!-------------------------------------------------------------------------------
! compute helicity from velocity and vorticity (both in x-space)
! INPUT:
!   u: velocity in x-space
!   vor: vorticity in x-space
! OUTPUT:
!   hel: helicty in x-space
!-------------------------------------------------------------------------------
subroutine helicity(u,vor,hel)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::hel(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix, iy, iz

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! helicity is scalar product of u and vorticity
        hel(ix,iy,iz) = (u(ix,iy,iz,1)*vor(ix,iy,iz,1) &
                      +  u(ix,iy,iz,2)*vor(ix,iy,iz,2) &
                      +  u(ix,iy,iz,3)*vor(ix,iy,iz,3))
      enddo
    enddo
  enddo

end subroutine helicity

!-------------------------------------------------------------------------------
! compute normalized helicity from velocity and vorticity (both in x-space)
! INPUT:
!   u: velocity in x-space
!   vor: vorticity in x-space
! OUTPUT:
!   hel: helicty in x-space
!-------------------------------------------------------------------------------
subroutine helicity_norm(u,vor,hel)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::hel(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  integer :: ix, iy, iz
  real(kind=pr) :: norm

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! normalization is product of norms
        norm = norm2( u(ix,iy,iz,:) ) * norm2( vor(ix,iy,iz,:) )

        if (norm >= 1.0d-5) then
          ! helicity is scalar product of u and vorticity
          hel(ix,iy,iz) = (u(ix,iy,iz,1)*vor(ix,iy,iz,1) &
                        +  u(ix,iy,iz,2)*vor(ix,iy,iz,2) &
                        +  u(ix,iy,iz,3)*vor(ix,iy,iz,3)) / norm
        endif
      enddo
    enddo
  enddo

end subroutine helicity_norm


! gradient of a scalar field in fourier space
subroutine gradient(uk,uk_dx,uk_dy,uk_dz)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),intent(inout):: uk,uk_dx,uk_dy,uk_dz
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)
           uk_dx(iz,iy,ix) = imag*kx*uk(iz,iy,ix)
           uk_dy(iz,iy,ix) = imag*ky*uk(iz,iy,ix)
           uk_dz(iz,iy,ix) = imag*kz*uk(iz,iy,ix)
       enddo
    enddo
  enddo
end subroutine gradient


! compute current dissipation rate, i.e. the los of energy due to friction
subroutine dissipation_rate( uk, epsilon )
  use vars
  use p3dfft_wrapper
  implicit none
  ! input: velocity vector in Fourier space
  complex(kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq),intent(inout):: uk
  ! output: dissipation rate
  real(kind=pr), intent(out)  :: epsilon
  ! local variables
  real(kind=pr) :: epsilon_loc, kx, ky, kz, k2, E
  integer :: ix, iy, iz, mpicode

  ! compute current dissipation rate. There are several ways to do this; one
  ! is computing the enstrophy in x-space, using the vorticity, but we can also
  ! stay with the velocity in k-space and integrate k^2 * E(k)
  ! However, this is exact only if E(k) is not averaged over the wavenumber shell
  ! as it is done when computing the spectrum.

  epsilon_loc = 0.d0

  do ix = ca(3),cb(3)
    kx = wave_x(ix)
    do iy = ca(2),cb(2)
      ky = wave_y(iy)
      do iz = ca(1),cb(1)
        kz = wave_z(iz)
        k2 = (kx*kx)+(ky*ky)+(kz*kz)

        if ( ix==0 .or. ix==nx/2 ) then
          E=dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2)/2. &
           +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2)/2. &
           +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)/2.
        else
          E=dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2) &
           +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2) &
           +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)
        endif

        epsilon_loc = epsilon_loc + k2 * E
      enddo
    enddo
  enddo

  epsilon_loc = 2.d0 * nu * epsilon_loc

  call MPI_ALLREDUCE(epsilon_loc,epsilon,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

end subroutine

end module basic_operators
