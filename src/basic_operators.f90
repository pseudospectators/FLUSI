!-------------------------------------------------------------------------------
! this module contains basic elementary operators, and as it is a module
! function overloading can be used.
! List of functions:
!       * curl
!-------------------------------------------------------------------------------
module basic_operators
 !-- interface for curl operators
 interface curl
  module procedure curl, curl_inplace, curl3
 end interface
 
 interface fieldmaxabs
  module procedure fieldmaxabs, fieldmaxabs3
 end interface
 
 
 contains
 
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


subroutine divergence( ink, outk )
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none
  complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
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
  
  
! returns the globally largest value of a given (real) field
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


! returns the globally smallest value of a given (real) field
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


! returns the globally largest value of a given vector field
! (L2-norm)
real(kind=pr) function fieldmaxabs3( inx )
  use mpi
  use vars
  implicit none
  real(kind=pr),intent(in):: inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr) :: max_local, max_global
  integer :: mpicode

  max_local = maxval( inx(:,:,:,1)*inx(:,:,:,1) + inx(:,:,:,2)*inx(:,:,:,2) &
            + inx(:,:,:,3)*inx(:,:,:,3) )
  max_local = dsqrt( max_local )
  call MPI_ALLREDUCE (max_local,max_global,1,&
       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! return the value    
  fieldmaxabs3 = max_global
end function fieldmaxabs3


! returns the globally largest value of a given vector field
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
  
end module basic_operators