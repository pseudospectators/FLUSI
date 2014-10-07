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
