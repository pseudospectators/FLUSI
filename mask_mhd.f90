! MHD wrapper for different mask functions 
subroutine create_mask_mhd()
  use mpi_header
  use mhd_vars
  implicit none

  mask=0.d0
  if(iPenalization == 1) then
     select case (iMask)
     case("TaylorCouette")
        call  tc_mask_mhd()
     case("smc")
        call smc_mask_mhd()
     case("smclinear")
        call smc_mask_mhd()
     case("smcflat")
        call smc_mask_mhd()
     case("smcnum")
        call smcnum_mask_mhd()
     case default
        if(mpirank == 0) then
           write (*,*) &
                "iMask not properly set for mhd in create_mask_mhd; stopping."
           stop
        endif
     end select
  endif
end subroutine create_mask_mhd

subroutine dealias(fk1,fk2,fk3) 
  use vars
  use mpi_header
  implicit none

  integer :: ix,iy,iz
  complex(kind=pr),intent(inout) :: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) :: fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr) :: kx2,ky2,kz2,t1,kxt2,kyt2,kzt2,kx_trunc,ky_trunc,kz_trunc
  integer :: i

  kx_trunc=(2.d0/3.d0)*dble(nx/2-1)
  ky_trunc=(2.d0/3.d0)*dble(ny/2-1)
  kz_trunc=(2.d0/3.d0)*dble(nz/2-1)  

  do iz=ca(1),cb(1)
     kz2 =(scalez*dble(modulo(iz+nz/2,nz)-nz/2))**2
     kzt2=dble(modulo(iz+nz/2,nz)-nz/2)/kz_trunc
     kzt2=kzt2*kzt2

     do ix=ca(2),cb(2)
        ! kx - x-wavenumber: 0..nx/2
        kx2=(scalex*dble(ix))**2
        kxt2=dble(ix)/kx_trunc
        kxt2=kxt2*kxt2

        do iy=ca(3),cb(3)
           ! ky - y-wavenumber: 0..ny/2-1 ,then, -ny/2..-1
           ky2=(scaley*dble(modulo(iy+ny/2,ny)-ny/2))**2
           kyt2=dble(modulo(iy+ny/2,ny)-ny/2)/ky_trunc
           kyt2=kyt2*kyt2

           if ((kxt2 + kyt2 + kzt2  .ge. 1.d0) .and. (iDealias==1)) then
              fk1(iz,ix,iy)=0.d0
              fk2(iz,ix,iy)=0.d0
              fk3(iz,ix,iy)=0.d0
           endif

        enddo
     enddo
  enddo

end subroutine dealias


! MHD wrapper for setting (possibly velocity-dependent) imposed field.
subroutine update_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none

  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  if(iPenalization == 1) then
     select case (iMask)
     case("TaylorCouette")
        call tc_us_mhd(ub)
     case("smc")
        call smc_us_mhd(ub)
     case("smclinear")
        call smclinear_us_mhd(ub)
     case("smcflat")
        call smcflat_us_mhd(ub)
     case("smcnum")
        !call smc_us_mhd(ub)
        call smcnum_us_mhd(ub)
     case default
        if(mpirank == 0) then
           write (*,*) &
                "iMask not properly set for mhd in update_us_mhd; stopping."
           stop
        endif
     end select
  endif
end subroutine update_us_mhd


! Set the solid velocity for MHD Taylor-Couette flow.
subroutine tc_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz
  
  us=0.d0

  ! Set the velocity field to be the steady-state solution:
  call taylor_couette_u_us(us(:,:,:,1),us(:,:,:,2),us(:,:,:,3))

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x + y*y)
        
        if(r <= R1) then
           do iz=ra(3),rb(3)
              ! Magnetic field:
              us(ix,iy,iz,4)=0.d0
              us(ix,iy,iz,5)=0.d0
              ! FIXME: non-penetration for b?
!!$              us(ix,iy,iz,4)=-omega1*y
!!$              us(ix,iy,iz,5)=omega1*x
           enddo
        endif
        if(r >= R2) then
           do iz=ra(3),rb(3)
              ! NB: We assume that the outer wall is not moving.
              ! Magnetic field: 
              ! FIXME: non-penetration for b?
              us(ix,iy,iz,4)=0.d0
              us(ix,iy,iz,5)=0.d0
           enddo
        endif
     enddo
  enddo

  ! Always penalize the z-component to the axial field.
  us(:,:,:,6)=B0 
end subroutine tc_us_mhd


! Set the mask function for MHD Taylor-Couette flow.
subroutine tc_mask_mhd()
  use mpi_header
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r,x,y
  integer :: ix,iy,iz

  mask=0.d0

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x + y*y)

        if(r <= R1 .or. r >= R2) then
           do iz=ra(3),rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine tc_mask_mhd


! Set the mask function for Sean-Montgomery-Chen flow.
subroutine smc_mask_mhd()
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz

  mask=0.d0

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +ay*y*y)

        if(r >= R1) then
           do iz=ra(3),rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine smc_mask_mhd


! Set the mask function for Sean-Montgomery-Chen flow.
subroutine smcnum_mask_mhd()
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz

  mask=0.d0

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +ay*y*y)

        if(r >= R1) then
           do iz=ra(3),rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine smcnum_mask_mhd



! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smc_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r,x,y
  integer :: ix,iy,iz
  real (kind=pr) :: a,b,c,d,k1,k2,h
  logical, save :: firstcall = .true. 
  
  if (firstcall) then
     firstcall = .false.
     
     ! Velocity is no-slip:
     us(:,:,:,1)=0.d0
     us(:,:,:,2)=0.d0
     us(:,:,:,3)=0.d0

     us(:,:,:,4)=0.d0
     us(:,:,:,5)=0.d0
     us(:,:,:,6)=B0

     k1=Bc*r2/r1
     k2=Bc/r1
     A=(2.d0*k1 -k2*(r2-r3))/(r3*r3*r3 -3.d0*r2*r3*r3 +3.d0*r2*r2*r3 -r2*r2*r2)
     B=(k2 -3.d0*A*(r2*r2 -r3*r3))/(2.d0*r2 -2.d0*r3)
     C=-3.d0*A*r3*r3 -2.d0*B*r3
     D=2.d0*A*r3*r3*r3 +B*r3*r3

     do ix=ra(1),rb(1)
        x=xl*(dble(ix)/dble(nx) -0.5d0)
        do iy=ra(2),rb(2)
           y=yl*(dble(iy)/dble(ny) -0.5d0)

           r=dsqrt(x*x +y*y)

           ! Linear profile:
           if(r >= r1 .and. r < r2) then
              do iz=ra(3),rb(3)
                 us(ix,iy,iz,4)=-Bc*y/r1
                 us(ix,iy,iz,5)=Bc*x/r1
              enddo
           endif

           ! Hermite profile:
           if(r >= r2 .and. r <= r3) then
              h=(A*r*r*r +B*r*r +C*r +D)
              do iz=ra(3),rb(3)
                 us(ix,iy,iz,4)=-h*y/r
                 us(ix,iy,iz,5)=h*x/r
              enddo
           endif

        enddo
     enddo

  endif ! only performed on first call.
end subroutine smc_us_mhd


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smclinear_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r,x,y
  integer :: ix,iy,iz
  real (kind=pr) :: a,b,c,d,k1,k2
  
  ! Velocity is no-slip:
  us(:,:,:,1)=0.d0
  us(:,:,:,2)=0.d0
  us(:,:,:,3)=0.d0

  us(:,:,:,4)=0.d0
  us(:,:,:,5)=0.d0
  us(:,:,:,6)=B0

  k1=Bc*r2/r1
  k2=Bc/r1
  A=(2.d0*k1 -k2*(r2-r3))/(r3*r3*r3 -3.d0*r2*r3*r3 +3.d0*r2*r2*r3 -r2*r2*r2)
  B=(k2 -3.d0*A*(r2*r2 -r3*r3))/(2.d0*r2 -2.d0*r3)
  C=-3.d0*A*r3*r3 -2.d0*B*r3
  D=2.d0*A*r3*r3*r3 +B*r3*r3

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)
        
        ! Linear profile:
        if(r >= r1 .and. r < r2) then
           do iz=ra(3),rb(3)
              us(ix,iy,iz,4)=-Bc*y/r1
              us(ix,iy,iz,5)=Bc*x/r1
           enddo
        endif
        
     enddo
  enddo
end subroutine smclinear_us_mhd


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smcflat_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r,x,y
  integer :: ix,iy,iz
  real (kind=pr) :: a,b,c,d,k1,k2
  
  ! Velocity is no-slip:
  us(:,:,:,1)=0.d0
  us(:,:,:,2)=0.d0
  us(:,:,:,3)=0.d0

  us(:,:,:,4)=0.d0
  us(:,:,:,5)=0.d0
  us(:,:,:,6)=B0

  k1=Bc*r2/r1
  k2=Bc/r1
  A=(2.d0*k1 -k2*(r2-r3))/(r3*r3*r3 -3.d0*r2*r3*r3 +3.d0*r2*r2*r3 -r2*r2*r2)
  B=(k2 -3.d0*A*(r2*r2 -r3*r3))/(2.d0*r2 -2.d0*r3)
  C=-3.d0*A*r3*r3 -2.d0*B*r3
  D=2.d0*A*r3*r3*r3 +B*r3*r3

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)
        
        ! Linear profile:
        if(r >= r1 .and. r < r2) then
           do iz=ra(3),rb(3)
              us(ix,iy,iz,4)=-Bc
              us(ix,iy,iz,5)=Bc
           enddo
        endif
        
     enddo
  enddo
end subroutine smcflat_us_mhd

subroutine bcpoint(on,x,y)
  use mpi_header
  use mhd_vars
  implicit none

  logical, intent(out) :: on
  real(kind=pr), intent(in) :: x,y
  real(kind=pr) :: f
  
  ! f is a sort of search width to find points representing \partial\Omega_f
  ! The larger the value of f, the more points.
  ! Too many points and the boundary has non-zero volume.
  ! NB: f values were determined by hand.
  ! FIXME: find a better way to determine f.
  f=0.d0
  if(nx == 16) f=9.d0
  if(nx == 32) f=15.d0
  if(nx == 64) f=15.d0 ! normally 29?
  if(nx == 128) f=30.d0 ! normally 60, but, needs to be less?
  if(nx == 256) f=115.d0
  if(f == 0.d0) then
     if (mpirank == 0) then
        write(*,*) "Fudge-factor not determined for given resolution."
     endif
     call exit
  endif

  on=.false.
  if(abs(x*x + ay*y*y -r1*r1) < f*dx*dy) on=.true.
end subroutine bcpoint

subroutine bcval(bcx,bcy,x,y)
  use mpi_header
  use mhd_vars
  implicit none
  real(kind=pr), intent(out) :: bcx,bcy
  real(kind=pr), intent(in) :: x,y
  real(kind=pr) n

  n=dsqrt(ay*ay*y*y+x*x)

  bcx=Bc*ay*y/n
  bcy=-Bc*x/n
end subroutine bcval

subroutine setpen(p1,p2,p3)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(out)::p1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out)::p2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out)::p3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  logical onboundary
  integer :: ix,iy,iz
  real(kind=pr) :: x,y,bcx,bcy

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        
        call bcpoint(onboundary,x,y)
        if(onboundary) then
           call bcval(bcx,bcy,x,y)
           !write(*,*) ix,iy
           do iz=ra(3),rb(3)
              p1(ix,iy,iz)=bcx
              p2(ix,iy,iz)=bcy
              p3(ix,iy,iz)=0.d0
           enddo
        else
           do iz=ra(3),rb(3)
              p1(ix,iy,iz)=0.d0
              p2(ix,iy,iz)=0.d0
              p3(ix,iy,iz)=0.d0
           enddo
        endif
     enddo
  enddo
end subroutine setpen

subroutine checkbc(diff,us1,us2,us3)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::us1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::us2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::us3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  logical onboundary
  integer :: ix,iy,iz
  real(kind=pr) :: x,y,bcx,bcy,perror,pnorm,ux,uy
  real(kind=pr),intent(out) :: diff

  diff=0.d0
!  write(*,*) -0.5d0*xl,-0.5d0*yl
!  write(*,*) 0.5d0*xl,0.5d0*yl
  
  do ix=ra(1),rb(1)
     x=xl*(dble(ix)/dble(nx) -0.5d0)


     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        
        call bcpoint(onboundary,x,y)
        if(onboundary) then
           call bcval(bcx,bcy,x,y)

           ! Un-comment to output points to check bcpoint validity.
           !write(*,*) x,y,bcx,bcy

           do iz=ra(3),rb(3)

              ux=us1(ix,iy,iz)
              uy=us2(ix,iy,iz)
              
              pnorm=bcx*bcx + bcy*bcy +1d-16
              perror=(ux-bcx)*(ux-bcx) +(uy-bcy)*(uy-bcy)
              
              perror=perror/pnorm
              if(perror > diff) diff=perror
           enddo
        endif
     enddo
  enddo
  !call abort ! Un-comment to stop after outputting points
end subroutine checkbc

! Compute the source for the pseudotime-stepper in physical space,
! returned in sx, sy, sz
subroutine pseudosource(ux,uy,uz,ukx,uky,ukz,sx,sy,sz)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ux(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::uy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::uz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  complex(kind=pr),intent(inout)::ukx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::uky(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout)::ukz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  real(kind=pr),intent(out)::sx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out)::sy(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out)::sz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  ! Local loop variables
  logical onboundary
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  real(kind=pr) :: x,y,bcx,bcy

  ! Compute gradient

  call fft(ukx,ux)
  call fft(uky,uy)
  call fft(ukz,uz)

  do iz=ca(1),cb(1)
     kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)

           k2=kx*kx +ky*ky +kz*kz
           
           ukx(iz,ix,iy)=-k2*ukx(iz,ix,iy)
           uky(iz,ix,iy)=-k2*uky(iz,ix,iy)
           ukz(iz,ix,iy)=-k2*ukz(iz,ix,iy)
        enddo
     enddo
  enddo
 
  ! call dealias(ukx,uky,ukz)

  call ifft(sx,ukx)
  call ifft(sy,uky)
  call ifft(sz,ukz)

  ! Compute penalisation
  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
  
        call bcpoint(onboundary,x,y)
        if(onboundary) then
           call bcval(bcx,bcy,x,y)
           !write(*,*) ix,iy
           do iz=ra(3),rb(3)
              sx(ix,iy,iz)=-(ux(ix,iy,iz) - bcx)/pseudoeps
              sy(ix,iy,iz)=-(uy(ix,iy,iz) - bcy)/pseudoeps
           enddo
        endif
     enddo
  enddo
  
end subroutine pseudosource

! Compute the Euclideian distance between points (ax,ay,az) and
! (bx,by,bz)
! FIXME: put this in some other "header" file for general use?
subroutine dist(ax,ay,az,bx,by,bz,d)
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: ax,ay,az,bx,by,bz
  real(kind=pr) :: d1,d2,d3
  real(kind=pr), intent(out) :: d

  d1=ax-bx
  d2=ay-by
  d3=az-bz

  d=dsqrt(d1*d1+d2*d2+d3*d3)
end subroutine dist

! Compute l-infinity distance between the real-valued arrays
! (ax,ay,az) and (bx,by,bz)
! FIXME: put this in some other "header" file for general use?
subroutine maxdist(ax,ay,az,bx,by,bz,d)
  use vars
  implicit none
  
  real(kind=pr),intent(in)::ax(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::ay(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::az(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  real(kind=pr),intent(in)::bx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::by(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in)::bz(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: dd
  real(kind=pr),intent(out) :: d
  integer :: ix,iy,iz

  d=0.d0
  do ix=ra(1),rb(1)
     do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
           call  dist(ax(ix,iy,iz),ay(ix,iy,iz),az(ix,iy,iz),&
                bx(ix,iy,iz),by(ix,iy,iz),bz(ix,iy,iz),dd)
           if(dd > d) d=dd
        enddo
     enddo
  enddo
end subroutine maxdist


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smcnum_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr) :: r,x,y,z,kx,ky,kz,k2,mydt,diff,diff0,globdiff,overthing,val
  real(kind=pr) :: vx,vy,vz
  integer :: ix,iy,iz
  integer :: i,myi,mpicode

  logical, save :: firstcall = .true. 
  ! the penalisation field
  real(kind=pr),dimension(:,:,:),allocatable :: usx, usy, usz
  ! for the 2-stage time-stepper
  real(kind=pr),dimension(:,:,:),allocatable :: tusx, tusy, tusz
  ! for computing gradient of penalisation field
  complex(kind=pr),dimension(:,:,:),allocatable :: uskx, usky, uskz
  ! PC time-stepping source buffers
  real(kind=pr),dimension(:,:,:),allocatable :: s1x, s1y, s1z
  real(kind=pr),dimension(:,:,:),allocatable :: s2x, s2y, s2z
  real(kind=pr) s1ms2
  
  ! Local loop variables, which, in modern languages, are declared
  ! locally in the loop
  real(kind=pr) :: bcx,bcy
  complex(kind=pr) :: ux,uy,uz
  complex(kind=pr) :: pkx,pky,pkz
  real (kind=pr) :: peps
  real (kind=pr) :: olddiff
  real (kind=pr) :: newdt
  real (kind=pr) :: lerror, error
  logical keeponkeepingon

  if (firstcall) then
     firstcall = .false.

     ! pseudo time-stepping parameters
     peps=pseudoeps
     mydt=pseudodt
     
     if (mpirank == 0) then
        write(*,*) "Computing penalization field via pseudo time-stepping...."
        write(*,*) "pseudoeps=",pseudoeps
        write(*,*) "pseudodt=",pseudodt
     endif

     myi=0 ! iteration variable

     call allocreal(usx)
     call allocreal(usy)
     call allocreal(usz)

     call allocreal(tusx)
     call allocreal(tusy)
     call allocreal(tusz)

     call allocreal(s1x)
     call allocreal(s1y)
     call allocreal(s1z)
     
     call allocreal(s2x)
     call allocreal(s2y)
     call allocreal(s2z)

     call alloccomplex(uskx)
     call alloccomplex(usky)
     call alloccomplex(uskz)
     
     keeponkeepingon=.true.

     ! initialize penalization field to zero
     usx=0.d0
     usy=0.d0
     usz=0.d0

     olddiff=0.d0
     diff0=0.d0

     if (mpirank == 0) write(*,*) "pseudodt diff0 error:"

     do while(keeponkeepingon) ! Solve for us

        myi=myi+1
        
        ! compute source for first stage:
        call pseudosource(usx,usy,usz,uskx,usky,uskz,s1x,s1y,s1z)
        ! perform the first stage
        tusx = usx +0.5d0*pseudodt*s1x
        tusy = usy +0.5d0*pseudodt*s1y
        tusz = usz +0.5d0*pseudodt*s1z
                
        ! compute source for second stage:
        call pseudosource(tusx,tusy,tusz,uskx,usky,uskz,s2x,s2y,s2z)
        ! perform the second stage
        usx = usx +pseudodt*(s2x)
        usy = usy +pseudodt*(s2y)
        usz = usz +pseudodt*(s2z)

        ! Project onto the solenoidal manifold
        call fft(uskx,usx)
        call fft(usky,usy)
        call fft(uskz,usz)
        call div_field_nul(uskx,usky,uskz)
        call ifft(usx,uskx)
        call ifft(usy,usky)
        call ifft(usz,uskz)
        
        ! output a sample bc point and what it should reach:
        ! ix=8
        ! iy=39
        ! x=xl*(dble(ix)/dble(nx) -0.5d0)
        ! y=yl*(dble(iy)/dble(ny) -0.5d0)
        ! call bcval(bcx,bcy,x,y)
        ! write(*,*)usx(ix,iy,0),bcx

        ! compute time-stepper error for adaptive method:
        call maxdist(s1x,s1y,s1z,s2x,s2y,s2z,diff)
        call MPI_REDUCE(diff,diff0,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
             MPI_COMM_WORLD,mpicode)


        ! check how close the field is to obeying the boundary conditions
        call checkbc(lerror,usx,usy,usz)
        call MPI_REDUCE(lerror,error,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
             MPI_COMM_WORLD,mpicode)

        if (mpirank == 0) then
           diff0=diff0*pseudodt

           ! I really hate this part of Fortran.
20         format (es10.2,x,es10.2,x,es10.2) 
           write(*,20) pseudodt,diff0,error
           
           if(diff0 < pseudoerrmin) pseudodt=1.4d0*pseudodt ! step too small
           if(diff0 > pseudoerrmax) pseudodt=0.7d0*pseudodt ! step too large

           ! Loop exit conditions:
           if(myi > 10 .and. error < dsqrt(eps)) then
              write(*,*) "Converged in ",myi," iterations:"
              write(*,*) error," < dsqrt(",peps,")=",dsqrt(eps)

              keeponkeepingon= .false.
           endif
           if(myi > 1000000) then ! we've gone too far: abort
              write(*,*) myi," is too many iterations."
              keeponkeepingon= .false.
              call abort
           endif
           if(myi > 10 .and. error > 100.d0) then
              ! convergence isn't happening: abort
              write(*,*) "error is greater than 100; aborting due to instability"
              keeponkeepingon= .false.
              call abort
           endif
        endif
        call MPI_BCAST(pseudodt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)
        call MPI_BCAST(keeponkeepingon,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpicode)
     enddo ! keep on keeping on?

     if (mpirank == 0) write(*,*) "finished!"

     if (mpirank == 0) write(*,*) "setting velocity penalty field...."
     ! Velocity is no-slip:
     us(:,:,:,1)=0.d0
     us(:,:,:,2)=0.d0
     us(:,:,:,3)=0.d0

     if (mpirank == 0) write(*,*) "setting magnetic penalty field...."
     ! copy ust to appropriate field for time-stepping. (us 4 and us 5)
     do ix=ra(1),rb(1)  
        do iy=ra(2),rb(2)
           do iz=ra(3),rb(3)
              us(ix,iy,iz,4)=usx(ix,iy,iz)
              us(ix,iy,iz,5)=usy(ix,iy,iz)
           enddo
        enddo
     enddo
     ! the z-component of the magnetic field is penalized to B0
     us(:,:,:,6)=b0 
     
     if (mpirank == 0) write(*,*) "deallocating temporary buffers..."
     ! Deallocate temporary buffers
     deallocate(usx,usy,usz)
     deallocate(tusx,tusy,tusz)
     deallocate(s1x,s1y,s1z)
     deallocate(s2x,s2y,s2z)
     deallocate(uskx,usky,uskz)

!     if (mpirank == 0) write(*,*) "Testing: aborted. (FIXME!)"
!     call exit
     if (mpirank == 0) write(*,*) "Finished setting penalty fields."

  end if ! if(firstcall)
  
end subroutine smcnum_us_mhd
