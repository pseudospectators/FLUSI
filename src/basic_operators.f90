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
   do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)

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
  integer :: mpicode, npoints

  ! number of points on local CPU
  npoints = (rb(1)-ra(1)+1)*(rb(2)-ra(2)+1)*(rb(3)-ra(3)+1)

  mean_local = sum(inx) / dble(npoints)
  call MPI_ALLREDUCE (mean_local,mean_global,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  ! return the value
  fieldmean = mean_global
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

  if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
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
  complex(kind=pr),intent(inout) :: fk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
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
          fk(iz,iy,ix,1:3)=0.d0
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
subroutine Vorticity2Velocity(vork,uk)
  use mpi
  use fsi_vars
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
  use fsi_vars
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
  use fsi_vars
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
  use fsi_vars
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
        norm = dsqrt(u(ix,iy,iz,1)**2 + u(ix,iy,iz,2)**2 + u(ix,iy,iz,3)**2) &
             * dsqrt(vor(ix,iy,iz,1)**2 + vor(ix,iy,iz,2)**2 + vor(ix,iy,iz,3)**2)
        ! helicity is scalar product of u and vorticity
        hel(ix,iy,iz) = (u(ix,iy,iz,1)*vor(ix,iy,iz,1) &
                      +  u(ix,iy,iz,2)*vor(ix,iy,iz,2) &
                      +  u(ix,iy,iz,3)*vor(ix,iy,iz,3)) / norm
      enddo
    enddo
  enddo

end subroutine helicity_norm

end module basic_operators
