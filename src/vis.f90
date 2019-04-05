!  Calculate visfusive term for time advancement,
!  exp(-nu*k^2*dt)
!  It is real valued, its global size is 0:nz-1, 0:nx/2, 0:ny-1 This
!  is computed only if the time step changes. Is does rarely, because
!  we round it to one digit.  also, dealiasing is done here (multiply
!  aliased avenumbers by zero).
!  Note the distinct fields 1..nf have different viscosities, which are stored
!  in the lin(1:nf) array
subroutine cal_vis(dt,vis)
  use vars
  use mpi
  use p3dfft_wrapper
  implicit none

  integer :: ix,iy,iz
  real(kind=pr),intent(inout) :: vis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(in) :: dt
  real(kind=pr) :: kx2,ky2,kz2,t1,kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc
  integer :: i

  t1=MPI_wtime()

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)

#ifdef __SX__
  do i=1,nf
#endif
  do ix=ca(3),cb(3)
     kx2=wave_x(ix)**2
     kxt2=(wave_x(ix)/scalex) / kx_trunc
     kxt2=kxt2*kxt2

     do iy=ca(2),cb(2)
        ky2=wave_y(iy)**2
        kyt2=(wave_y(iy)/scaley) / ky_trunc
        kyt2=kyt2*kyt2

        do iz=ca(1),cb(1)
           kz2 =wave_z(iz)**2
           kzt2=(wave_z(iz)/scalez) / kz_trunc
           kzt2=kzt2*kzt2

           ! Dealiasing is done here for reasons of efficiency
#ifndef __SX__
           do i=1,nf
#endif
              if ((kxt2 + kyt2 + kzt2  .ge. 1.d0) .and. (iDealias==1)) then
                 vis(iz,iy,ix,i)=0.d0
              else
                 vis(iz,iy,ix,i)=dexp( -dt*lin(i)*(kx2 + ky2 + kz2) )
              endif
#ifndef __SX__
           enddo
#endif

        enddo
     enddo
  enddo
#ifdef __SX__
  enddo
#endif

  call toc("cal_vis", MPI_wtime() - t1)
end subroutine cal_vis
