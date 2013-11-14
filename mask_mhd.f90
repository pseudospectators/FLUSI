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
        call smc_mask_mhd()
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
  use mpi_header
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz

  mask=0.d0

  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)

        if(r >= R1) then
           do iz=ra(3),rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine smc_mask_mhd


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smc_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r,x,y
  integer :: ix,iy,iz
  real (kind=pr) :: a,b,c,d,k1,k2,h
  
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
        WRITE(*,*) "Fudge-factor not determined for given resolution."
     endif
     call exit
  endif

  on=.false.
  if(abs(x*x + y*y -r1*r1) < f*dx*dy) on=.true.
end subroutine bcpoint

subroutine bcval(bcx,bcy,x,y)
  use mpi_header
  use mhd_vars
  implicit none
  real(kind=pr), intent(out) :: bcx,bcy
  real(kind=pr), intent(in) :: x,y
  
  bcx=Bc*y/r1
  bcy=-Bc*x/r1
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
  integer a

  a=0
  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        
        call bcpoint(onboundary,x,y)
        if(onboundary) then
           a = a+1
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
  
  do ix=ra(1),rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        
        call bcpoint(onboundary,x,y)
        if(onboundary) then 
           call bcval(bcx,bcy,x,y)

           ! Un-comment to output points to check bcpoint validity.
           !write(*,*) x,y 

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
!  call abort ! Un-comment to stop after outputting points
end subroutine checkbc


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smcnum_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr) :: r,x,y,z,kx,ky,kz,k2,mydt,diff,diff0,globdiff,overthing,val
  real(kind=pr) :: vx,vy,vz
  integer :: ix,iy,iz,i,myi,mpicode

  logical, save :: FirstCall = .TRUE. 
  complex(kind=pr),dimension(:,:,:),allocatable :: pk1, pk2, pk3
  complex(kind=pr),dimension(:,:,:),allocatable :: ust1, ust2, ust3
  real(kind=pr),dimension(:,:,:),allocatable :: p1, p2, p3
  ! Local loop variables, which, in modern languages, are declared
  ! locally in the loop
  real(kind=pr) :: bcx,bcy
  complex(kind=pr) :: ux,uy,uz
  complex(kind=pr) :: pkx,pky,pkz
  real (kind=pr) :: peps
  logical keeponkeepingon

  if (FirstCall) then
     FirstCall = .FALSE.

     ! pseudo time-stepping parameters
     peps=pseudoeps
     mydt=pseudodt
     
     write(*,*) "Computing penalization field via pseudo time-stepping...."
     write(*,*) "pseudoeps=",pseudoeps
     write(*,*) "pseudodt=",pseudodt

     myi=0 ! iteration variable

     allocate(p1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
     allocate(p2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
     allocate(p3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

     allocate(pk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
     allocate(pk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
     allocate(pk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
     
     allocate(ust1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
     allocate(ust2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
     allocate(ust3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
    
     keeponkeepingon=.true.

     ! initialize penalization field to zero
     ust1=0.d0
     ust2=0.d0
     ust3=0.d0
     
     do while(keeponkeepingon) ! Solve for ust

        ! Compute the penalization term
        call setpen(p1,p2,p3)
        ! Transform the penalization term to Fourier space:
        call fft(pk1,p1)
        call fft(pk2,p2)
        call fft(pk3,p3)
        call dealias(pk1,pk2,pk3)

        ! compute the gradient and time-step
        do iz=ca(1),cb(1)
           kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           do ix=ca(2),cb(2)
              kx=scalex*dble(ix)
              do iy=ca(3),cb(3)
                 ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
                 
                 ux=ust1(iz,ix,iy)
                 uy=ust2(iz,ix,iy)
                 uz=ust3(iz,ix,iy)

                 k2=kx*kx +ky*ky +kz*kz

                 pkx=pk1(iz,ix,iy)
                 pky=pk2(iz,ix,iy)
                 pkz=pk3(iz,ix,iy)

                 ust1(iz,ix,iy)=ux +mydt*(-k2*ux +pkx/peps)
                 ust2(iz,ix,iy)=uy +mydt*(-k2*uy +pky/peps)
                 ust3(iz,ix,iy)=uz +mydt*(-k2*uz +pkz/peps)
              enddo
           enddo
        enddo

        call dealias(ust1,ust2,ust3)

        ! project onto the solenoidal manifold
        call div_field_nul(ust1,ust2,ust3)

        myi=myi+1

        ! transform ust to physical space
        call ifft(p1,ust1)
        call ifft(p2,ust2)
        call ifft(p3,ust3)
        ! check how close the field is to obeying the boundary conditions
        call checkbc(diff,p1,p2,p3)
        call MPI_REDUCE(diff,diff0,&
             1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
             MPI_COMM_WORLD,mpicode)

        ! output a sample bc point and what it should reach:
        ! ix=8
        ! iy=39
        ! x=xl*(dble(ix)/dble(nx) -0.5d0)
        ! y=yl*(dble(iy)/dble(ny) -0.5d0)
        ! call bcval(bcx,bcy,x,y)
        ! write(*,*)p1(ix,iy,0),bcx

        if (mpirank == 0) then
           write(*,*) "error=",diff0
           ! Loop exit conditions:
           if(myi > 10 .and. diff0 < dsqrt(eps)) then
              write(*,*) "finished: ",diff0," < ",peps
              keeponkeepingon= .FALSE.
           endif
           if(myi > 20000) then ! we've gone too far: abort
              if (mpirank == 0) WRITE(*,*) myi," is too many iterations."
              keeponkeepingon= .FALSE.
           endif
        endif
     enddo ! keep on keeping on?

     ! copy ust to appropriate field for time-stepping. (us 4 and us 5)
     do ix=ra(1),rb(1)  
        do iy=ra(2),rb(2)
           do iz=ra(3),rb(3)
              us(ix,iy,iz,4)=0.d0
              us(ix,iy,iz,5)=0.d0
           enddo
        enddo
     enddo
     
     deallocate(p1,p2,p3)
     deallocate(pk1,pk2,pk3)
     deallocate(ust1,ust2,ust3)
     
     if (mpirank == 0) WRITE(*,*) "Testing: aborted. (FIXME!)"
     call exit

  end if
  
  ! Velocity is no-slip:
  us(:,:,:,1)=0.d0
  us(:,:,:,2)=0.d0
  us(:,:,:,3)=0.d0

  ! the z-component of the magnetic field is penalized to B0
  us(:,:,:,6)=B0 
  
end subroutine smcnum_us_mhd
