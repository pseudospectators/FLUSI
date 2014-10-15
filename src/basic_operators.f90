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

  use ghosts
  use vars

  !-- interface for curl operators
  interface curl
    module procedure curl, curl_inplace, curl3, curl_x
  end interface

  !-- interface for maxima of fields
  interface fieldmaxabs
    module procedure fieldmaxabs, fieldmaxabs3
  end interface

  !-- check fields for NaN 
  interface checknan
    module procedure checknan_cmplx, checknan_real
  end interface
  
  ! divergence (in x- or k- space)
  interface divergence
    module procedure divergence_k, divergence_x
  end interface
 
 
!!!!!!!!!!! 
 contains
!!!!!!!!!!! 
 
subroutine dealias(fk1,fk2,fk3) 
  use vars
  use mpi
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: kx2,ky2,kz2,kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)  

  do iz=ca(1),cb(1)
     kz2 =wave_z(iz)**2
     kzt2=(wave_z(iz)/scalez) / kz_trunc
     kzt2=kzt2*kzt2
     
     do iy=ca(2),cb(2)
        ky2=wave_y(iy)**2
        kyt2=(wave_y(iy)/scaley) / ky_trunc
        kyt2=kyt2*kyt2

        do ix=ca(3),cb(3)
           kx2=wave_x(ix)**2
           kxt2=(wave_x(ix)/scalex) / kx_trunc
           kxt2=kxt2*kxt2

           if ((kxt2 + kyt2 + kzt2  .ge. 1.d0)) then
              fk1(iz,iy,ix)=0.d0
              fk2(iz,iy,ix)=0.d0
              fk3(iz,iy,ix)=0.d0
           endif

        enddo
     enddo
  enddo

end subroutine dealias

subroutine dealias1(fk1)
  use vars
  use mpi
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: kx2,ky2,kz2,kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)  

  do iz=ca(1),cb(1)
     kz2 =wave_z(iz)**2
     kzt2=(wave_z(iz)/scalez) / kz_trunc
     kzt2=kzt2*kzt2
     
     do iy=ca(2),cb(2)
        ky2=wave_y(iy)**2
        kyt2=(wave_y(iy)/scaley) / ky_trunc
        kyt2=kyt2*kyt2

        do ix=ca(3),cb(3)
           kx2=wave_x(ix)**2
           kxt2=(wave_x(ix)/scalex) / kx_trunc
           kxt2=kxt2*kxt2

           if ((kxt2 + kyt2 + kzt2  .ge. 1.d0)) then
              fk1(iz,iy,ix)=0.d0
           endif

        enddo
     enddo
  enddo

end subroutine dealias1
 
 
 
! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl(out1,out2,out3,in1,in2,in3)
  use p3dfft_wrapper
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
  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
           kx=wave_x(ix)

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
  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
           kx=wave_x(ix)
           
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
  use p3dfft_wrapper
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(inout)::fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)
  
  ! Compute curl of given field in Fourier space:
  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
           kx=wave_x(ix)

           fk(iz,iy,ix,1)=imag*(ky*fk(iz,iy,ix,3) -kz*fk(iz,iy,ix,2))
           fk(iz,iy,ix,2)=imag*(kz*fk(iz,iy,ix,1) -kx*fk(iz,iy,ix,3))
           fk(iz,iy,ix,3)=imag*(kx*fk(iz,iy,ix,2) -ky*fk(iz,iy,ix,1))
        enddo
     enddo
  enddo
end subroutine curl3_inplace 


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional. The precision is reduced
! to second order in space, since this kind of filtering may help in postprocessing
subroutine curl_2nd (fx,fy,fz)
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
  do iz=ca(1),cb(1)
     kz=dsin(dz*wave_z(iz))/ dz ! (reduced to 2nd order)
     do iy=ca(2),cb(2)
        ky=dsin(dy*wave_y(iy))/dy ! (reduced to 2nd order)
        do ix=ca(3),cb(3)
           kx=dsin(dx*wave_x(ix))/dx ! (reduced to 2nd order)
           
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
  use p3dfft_wrapper
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)
  
  ! Compute curl of given field in Fourier space:
  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
           kx=wave_x(ix)

           outk(iz,iy,ix,1)=imag*(ky*ink(iz,iy,ix,3) -kz*ink(iz,iy,ix,2))
           outk(iz,iy,ix,2)=imag*(kz*ink(iz,iy,ix,1) -kx*ink(iz,iy,ix,3))
           outk(iz,iy,ix,3)=imag*(kx*ink(iz,iy,ix,2) -ky*ink(iz,iy,ix,1))
        enddo
     enddo
  enddo
end subroutine curl3 


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl_x( u, rotu )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(out)::rotu(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,dxinv,dyinv,dzinv,a1,a2,a4,a5
  real(kind=pr) :: uxdy,uxdz,uydx,uydz,uzdx,uzdy
  ! input vector field in Fourier space
  complex(kind=pr),allocatable,dimension(:,:,:,:)::ink
  ! output scalar field in Fourier space
  complex(kind=pr),allocatable,dimension(:,:,:,:)::outk
  complex(kind=pr) :: imag   ! imaginary unit
  
  
  
  select case(method)
  case('spectral')  
      !-------------------------------------------------------------------------
      ! compute divergence(u) using spectral derivatives. Note that both input
      ! and output of this function are in x-space, which requires 3 fft and
      ! 1 ifft. 
      !-------------------------------------------------------------------------
      allocate(ink (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
      allocate(outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
      ! transform input to Fourier space
      call fft3( inx=u, outk=ink )
      
      imag = dcmplx(0.d0,1.d0)
      ! Compute curl of given field in Fourier space:
      do iz=ca(1),cb(1)
        kz=wave_z(iz)
        do iy=ca(2),cb(2)
            ky=wave_y(iy)
            do ix=ca(3),cb(3)
              kx=wave_x(ix)
              outk(iz,iy,ix,1)=imag*(ky*ink(iz,iy,ix,3) -kz*ink(iz,iy,ix,2))
              outk(iz,iy,ix,2)=imag*(kz*ink(iz,iy,ix,1) -kx*ink(iz,iy,ix,3))
              outk(iz,iy,ix,3)=imag*(kx*ink(iz,iy,ix,2) -ky*ink(iz,iy,ix,1))
            enddo
        enddo
      enddo
      
      ! transform output back to x-space
      call ifft3( ink=outk, outx=rotu(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )
      deallocate(ink,outk)
      
  case('centered_2nd')
      call synchronize_ghosts( u )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
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
      call synchronize_ghosts( u )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
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
      call suicide2('invalid METHOD in curl_x:'//method)
  end select
end subroutine curl_x


!-------------------------------------------------------------------------------
! compute divergence of vector valued field (ink, 4D array)
! and returns it in outk
!-------------------------------------------------------------------------------
subroutine divergence_k( ink, outk )
  use p3dfft_wrapper
  implicit none
  ! input vector field in Fourier space
  complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  ! output scalar field in Fourier space
  complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz

  do iz=ca(1),cb(1)          
    !-- wavenumber in z-direction
    kz = wave_z(iz)       
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)      
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)
        outk(iz,iy,ix)=kx*ink(iz,iy,ix,1)+ky*ink(iz,iy,ix,2)+kz*ink(iz,iy,ix,3)
        outk(iz,iy,ix)=dcmplx(0.d0,1.d0)*outk(iz,iy,ix)
      enddo
    enddo
  enddo
end subroutine divergence_k

!-------------------------------------------------------------------------------
! compute the divergence of the input field u and return it in divu
! depending on the value of "method", different discretization
! is used.
!-------------------------------------------------------------------------------
subroutine divergence_x( u, divu )
  use p3dfft_wrapper
  implicit none
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(inout)::divu(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,dxinv,dyinv,dzinv,uxdx,uydy,uzdz,a1,a2,a4,a5
  ! input vector field in Fourier space
  complex(kind=pr),allocatable,dimension(:,:,:,:)::ink
  ! output scalar field in Fourier space
  complex(kind=pr),allocatable,dimension(:,:,:)::outk
  
  
  
  select case(method)
  case('spectral')  
      !-------------------------------------------------------------------------
      ! compute divergence(u) using spectral derivatives. Note that both input
      ! and output of this function are in x-space, which requires 3 fft and
      ! 1 ifft. 
      !-------------------------------------------------------------------------
      allocate(ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
      allocate(outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
      ! transform input to Fourier space
      call fft3( inx=u, outk=ink )
      do iz=ca(1),cb(1)          
        !-- wavenumber in z-direction
        kz = wave_z(iz)       
        do iy=ca(2), cb(2)
          !-- wavenumber in y-direction
          ky = wave_y(iy)      
          do ix=ca(3), cb(3)
            !-- wavenumber in x-direction
            kx = wave_x(ix)
            outk(iz,iy,ix)=kx*ink(iz,iy,ix,1)+ky*ink(iz,iy,ix,2)+kz*ink(iz,iy,ix,3)
            outk(iz,iy,ix)=dcmplx(0.d0,1.d0)*outk(iz,iy,ix)
          enddo
        enddo
      enddo
      ! transform output back to x-space
      call ifft( ink=outk, outx=divu(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
      deallocate(ink,outk)
      
  case('centered_2nd')
      call synchronize_ghosts_FD( u )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
      !-------------------------------------------------------------------------
      dxinv = 1.d0/(2.d0*dx)
      dyinv = 1.d0/(2.d0*dy)
      dzinv = 1.d0/(2.d0*dz)
      
      if (nx>1) then
        ! three-dimensional simulation
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              uxdx = dxinv*(u(ix+1,iy,iz,1)-u(ix-1,iy,iz,1))
              uydy = dyinv*(u(ix,iy+1,iz,2)-u(ix,iy-1,iz,2))
              uzdz = dzinv*(u(ix,iy,iz+1,3)-u(ix,iy,iz-1,3))
              divu(ix,iy,iz) = uxdx + uydy + uzdz
            enddo
          enddo
        enddo 
      elseif (nx==1) then
        ! two-dimensional simulation
        ix=0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydy = dyinv*(u(ix,iy+1,iz,2)-u(ix,iy-1,iz,2))
            uzdz = dzinv*(u(ix,iy,iz+1,3)-u(ix,iy,iz-1,3))
            divu(ix,iy,iz) = uydy + uzdz
          enddo
        enddo 
      endif
      
  case('centered_4th')
      call synchronize_ghosts_FD( u )
      !-------------------------------------------------------------------------
      ! compute divergence(u) using second order centered period FD
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
              uxdx = (a1*u(ix-2,iy,iz,1)+a2*u(ix-1,iy,iz,1)+a4*u(ix+1,iy,iz,1)+a5*u(ix+2,iy,iz,1))*dxinv
              uydy = (a1*u(ix,iy-2,iz,2)+a2*u(ix,iy-1,iz,2)+a4*u(ix,iy+1,iz,2)+a5*u(ix,iy+2,iz,2))*dyinv
              uzdz = (a1*u(ix,iy,iz-2,3)+a2*u(ix,iy,iz-1,3)+a4*u(ix,iy,iz+1,3)+a5*u(ix,iy,iz+2,3))*dzinv
              divu(ix,iy,iz) = uxdx + uydy + uzdz
            enddo
          enddo
        enddo 
      elseif (nx==1) then
        ! two-dimensional simulation
        ix = 0
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            uydy = (a1*u(ix,iy-2,iz,2)+a2*u(ix,iy-1,iz,2)+a4*u(ix,iy+1,iz,2)+a5*u(ix,iy+2,iz,2))*dyinv
            uzdz = (a1*u(ix,iy,iz-2,3)+a2*u(ix,iy,iz-1,3)+a4*u(ix,iy,iz+1,3)+a5*u(ix,iy,iz+2,3))*dzinv
            divu(ix,iy,iz) = uydy + uzdz
          enddo
        enddo
      endif
      
  case default
    call suicide2('invalid METHOD in divergence:'//method)
  end select
end subroutine divergence_x



! computes laplace(ink) for a scalar valued field and returns it in the same array
subroutine laplacien_inplace( ink )
  use p3dfft_wrapper
  implicit none
  complex(kind=pr),intent(inout)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2

  do iz=ca(1),cb(1)          
    !-- wavenumber in z-direction
    kz = wave_z(iz)
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)
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
  use p3dfft_wrapper
  implicit none
  complex(kind=pr),intent(inout)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2

  do iz=ca(1),cb(1)          
    !-- wavenumber in z-direction (reduced to 2nd order)
    kz = dsin( dz*wave_z(iz) )/ dz
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction  (reduced to 2nd order)
      ky = dsin( dy*wave_y(iy) )/dy
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction  (reduced to 2nd order)
        kx = dsin(dx*wave_x(ix))/dx
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


! returns the globally largest entry of a given vector field
! (L2-norm)
real(kind=pr) function fieldmaxabs3( inx )
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr) :: max_local, max_global, value
  integer :: mpicode
  integer ::ix,iy,iz
  max_local = 0.d0
  do ix = ra(1), rb(1)
    do iy = ra(2), rb(2)
       do iz = ra(3), rb(3)
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
  fieldmaxabs3 = max_global
end function fieldmaxabs3


! returns the globally largest entry of a given vector field
! (L2-norm)
real(kind=pr) function fieldmaxabs( inx1, inx2, inx3 )
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
  fieldmaxabs = max_global
end function fieldmaxabs
  

!-------------------------------------------------------------------------------
! check a real valued field for NaNs and display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_real( field, msg )
  implicit none
  real(kind=pr),intent(in)::field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  character(len=*),intent(in)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz
  foundnan = 0
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)  
        if (is_nan(field(ix,iy,iz))) foundnan = 1
      enddo
    enddo
  enddo  
  
  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)  

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans      
end subroutine checknan_real


!-------------------------------------------------------------------------------
! check a complex field for NaN's, display warning if found
!-------------------------------------------------------------------------------
subroutine checknan_cmplx( field, msg )
  implicit none
  complex(kind=pr),intent(in)::field(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  character(len=*),intent(in)::msg
  integer :: foundnan,foundnans,mpicode,ix,iy,iz
  foundnan = 0
  do iz=ca(1),cb(1)
    do iy=ca(2),cb(2)
      do ix=ca(3),cb(3)
        if (is_nan(aimag(field(iz,iy,ix)))) foundnan = 1
        if (is_nan(real (field(iz,iy,ix)))) foundnan = 1
      enddo
    enddo
  enddo  
  
  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)  

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
end subroutine checknan_cmplx  

!-------------------------------------------------------------------------------
! computes the volume integral of the scalar quantity u 
!-------------------------------------------------------------------------------
real(kind=pr) function  volume_integral( u )
  use p3dfft_wrapper
  implicit none

  ! input/output field in x-space
  real(kind=pr),intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  integer::ix,iy,iz,mpicode
  real(kind=pr)::int_local,dxyz
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        int_local = int_local + u(ix,iy,iz)
      enddo
    enddo
  enddo    
  
  int_local = int_local*dxyz

  call MPI_ALLREDUCE (int_local,volume_integral,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
end function volume_integral


end module basic_operators




