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
    module procedure curl, curl_inplace, curl3
  end interface

  !-- interface for maxima of fields
  interface fieldmaxabs
    module procedure fieldmaxabs, fieldmaxabs3
  end interface

  !-- check fields for NaN
  interface checknan
    module procedure checknan_cmplx, checknan_real
  end interface


!!!!!!!!!!!
 contains
!!!!!!!!!!!


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl(out1,out2,out3,in1,in2,in3)
  use mpi
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
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  ! input/output field in Fourier space
  complex(kind=pr),intent(inout)::fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz
  complex(kind=pr) :: imag   ! imaginary unit
  complex(kind=pr) :: t1,t2,t3 ! temporary loop variables

  imag = dcmplx(0.d0,1.d0)

  ! Compute curl of given field in Fourier space:
  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
           kx=wave_x(ix)

           t1=fk(iz,iy,ix,1)
           t2=fk(iz,iy,ix,2)
           t3=fk(iz,iy,ix,3)

           fk(iz,iy,ix,1)=imag*(ky*t3 - kz*t2)
           fk(iz,iy,ix,2)=imag*(kz*t1 - kx*t3)
           fk(iz,iy,ix,3)=imag*(kx*t2 - ky*t1)
        enddo
     enddo
  enddo
end subroutine curl3_inplace


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
  use mpi
  use p3dfft_wrapper
  use vars
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


! returns the globally largest entry of a given vector field
! (L2-norm)
real(kind=pr) function fieldmaxabs3( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr) :: max_local, max_global, value
  integer :: mpicode
  integer ::ix,iy,iz
  max_local = 0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
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
  fieldmaxabs = max_global
end function fieldmaxabs


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

  call MPI_ALLREDUCE (foundnan,foundnans,1,MPI_INTEGER,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
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

end module basic_operators
