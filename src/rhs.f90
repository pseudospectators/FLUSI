subroutine rhs_acm_spectral(u,rhs)
  use mpi 
  use fsi_vars
  use diff
  use p3dfft_wrapper
  implicit none
  real(kind=pr),intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  complex(kind=pr),dimension(:,:,:,:),allocatable :: nlk
  real(kind=pr),dimension(:,:,:,:),allocatable :: vort,u2
  complex(kind=pr),dimension(:,:,:,:),allocatable :: workc
  
  rhs=0.d0
  
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))   
  allocate(u2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))   
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) )
  
  call fft(inx=u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1),outk=uk(:,:,:,1))
  call fft(inx=u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2),outk=uk(:,:,:,2))
  call fft(inx=u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3),outk=uk(:,:,:,3))
  call fft(inx=u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4),outk=uk(:,:,:,4))
  

  call cal_nlk(0.d0,0,nlk,uk,u2,vort,workc)

  
  
  call ifft(ink=nlk(:,:,:,1),outx=rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1))
  call ifft(ink=nlk(:,:,:,2),outx=rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2))
  call ifft(ink=nlk(:,:,:,3),outx=rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3))
  call ifft(ink=nlk(:,:,:,4),outx=rhs(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4))
  
  deallocate(uk,nlk,vort,workc,u2)
end subroutine rhs_acm_spectral



! Wrapper for computing the nonlinear source term for Navier-Stokes/MHD
subroutine cal_nlk(time,it,nlk,uk,u,vort,workc)
  use fsi_vars
  use p3dfft_wrapper
  use basic_operators
  implicit none

  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) 
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(in) :: time
  real(kind=pr) :: t1,t0,kx,ky,kz
  integer, intent(in) :: it
  integer :: ix,iy,iz
  complex(kind=pr) :: imag,pk   ! imaginary unit
  imag = dcmplx(0.d0,1.d0)
  
  call cal_nlk_fsi(time,it,nlk,uk,u,vort,workc)
  ! add pressure gradient
  do iz=ca(1),cb(1)
    kz=wave_z(iz)
    do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
          kx=wave_x(ix) 
          
          pk = uk(iz,iy,ix,4)
          
          nlk(iz,iy,ix,1) = nlk(iz,iy,ix,1)-imag*kx*pk
          nlk(iz,iy,ix,2) = nlk(iz,iy,ix,2)-imag*ky*pk
          nlk(iz,iy,ix,3) = nlk(iz,iy,ix,3)-imag*kz*pk
        enddo
    enddo
  enddo
  
  ! compute pressure RHS
  call divergence( uk(:,:,:,1:3), workc(:,:,:,1) )
  nlk(:,:,:,4) = -c_0**2 * workc(:,:,:,1)
        
end subroutine cal_nlk


subroutine cal_nlk_fsi(time,it,nlk,uk,u,vort,workc)
  use mpi
  use p3dfft_wrapper
  use fsi_vars
  use basic_operators
  implicit none

  real(kind=pr),intent (in) :: time
  integer, intent(in) :: it
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(out)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq) 
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr) :: t0,t1,ux,uy,uz,vorx,vory,vorz,chi,usx,usy,usz
  real(kind=pr) :: fx,fy,fz,fx1,fy1,fz1
  real(kind=pr) :: soft_startup
  integer :: ix,iz,iy,mpicode

  !-----------------------------------------------------------------------------
  !-- Calculate velocity in physical space
  !-----------------------------------------------------------------------------
  call ifft3 (outx=u(:,:,:,1:3), ink=uk(:,:,:,1:3))
  
  !-----------------------------------------------------------------------------
  !-- Compute vorticity
  !-----------------------------------------------------------------------------
  call curl ( ink=uk(:,:,:,1:3), outk=nlk(:,:,:,1:3) )
  call ifft3 ( ink=nlk(:,:,:,1:3), outx=vort(:,:,:,1:3))  
  
    
  !-----------------------------------------------------------------------------
  !-- Non-Linear terms
  !-----------------------------------------------------------------------------
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)  
        ! local loop variables
        ux   = u(ix,iy,iz,1)
        uy   = u(ix,iy,iz,2)
        uz   = u(ix,iy,iz,3)
        vorx = vort(ix,iy,iz,1)
        vory = vort(ix,iy,iz,2)
        vorz = vort(ix,iy,iz,3)
        
        ! we overwrite the vorticity with the NL terms in phys space
        ! note this is indeed -(vor x u) (negative sign)
        vort(ix,iy,iz,1) = uy*vorz - uz*vory 
        vort(ix,iy,iz,2) = uz*vorx - ux*vorz 
        vort(ix,iy,iz,3) = ux*vory - uy*vorx 
      enddo
    enddo
  enddo
  ! to Fourier space
  call fft3( inx=vort(:,:,:,1:3),outk=nlk(:,:,:,1:3) )  
  
  workc(:,:,:,1:3)=uk(:,:,:,1:3)
  call laplacien_inplace(workc(:,:,:,1))
  call laplacien_inplace(workc(:,:,:,2))
  call laplacien_inplace(workc(:,:,:,3))
  
  nlk(:,:,:,1)=nlk(:,:,:,1)+nu*workc(:,:,:,1)
  nlk(:,:,:,2)=nlk(:,:,:,2)+nu*workc(:,:,:,2)
  nlk(:,:,:,3)=nlk(:,:,:,3)+nu*workc(:,:,:,3)
  
  
end subroutine cal_nlk_fsi



!-------------------------------------------------------------------------------
! Compute the pressure. It is given by the divergence of the non-linear
! terms (nlk: intent(in)) divided by k**2.
! so: p=(i*kx*sxk + i*ky*syk + i*kz*szk) / k**2 
! note: we use rotational formulation: p is NOT the physical pressure
! as the RHS in cal_nlk is on the left side, i.e. -(vor x u) -chi*(u-us)
! its sign is inversed when computing the pressure
!-------------------------------------------------------------------------------
subroutine pressure(nlk,pk)
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr),intent(out):: pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr) :: imag,nlkx,nlky,nlkz

  ! as the RHS in cal_nlk is on the left side, i.e. -(vor x u) -chi*(u-us)
  ! its sign is inversed when computing the pressure
  
  imag = dcmplx(0.d0,1.d0)

  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
          kx=wave_x(ix)
          k2=kx*kx + ky*ky + kz*kz
          if(k2 .ne. 0.0) then
            ! contains the pressure in Fourier space
            ! note "-" sign
            nlkx = nlk(iz,iy,ix,1)
            nlky = nlk(iz,iy,ix,2)
            nlkz = nlk(iz,iy,ix,3)
            pk(iz,iy,ix) = -imag*(kx*nlkx + ky*nlky + kz*nlkz) / k2
          else
            pk(iz,iy,ix) = dcmplx(0.d0,0.d0)
          endif
      enddo
    enddo
  enddo
end subroutine pressure


!-------------------------------------------------------------------------------
! Compute the pressure in an inefficient way given only the velocity
!-------------------------------------------------------------------------------
subroutine pressure_given_uk(time,u,uk,nlk,vort,work,workc,press)
  use mpi
  use p3dfft_wrapper
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
  
!   call cal_nlk_fsi (time,0,nlk,uk,u,vort,work,workc) 
!   call pressure(nlk,workc(:,:,:,1))
!   call ifft(ink=workc(:,:,:,1), outx=press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  
end subroutine 



!-------------------------------------------------------------------------------
! Add the gradient of the pressure to the nonlinear term, which is the actual
! projection scheme used in this code. The non-linear term comes in with NL and
! penalization and leaves divergence free
!-------------------------------------------------------------------------------
subroutine add_grad_pressure(nlk1,nlk2,nlk3)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr),intent(inout):: nlk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: nlk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: nlk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr) :: qk, nlx,nly,nlz
  complex(kind=pr) :: imag   ! imaginary unit

!   imag = dcmplx(0.d0,1.d0)
!   
!   do iz=ca(1),cb(1)
!      kz=wave_z(iz)
!      do iy=ca(2),cb(2)
!         ky=wave_y(iy)
!         do ix=ca(3),cb(3)
!            kx=wave_x(ix)
!            
!            k2=kx*kx + ky*ky + kz*kz
! 
!            if (k2 .ne. 0.0) then
!               nlx=nlk1(iz,iy,ix)
!               nly=nlk2(iz,iy,ix)
!               nlz=nlk3(iz,iy,ix)
! 
!               ! qk is the Fourier coefficient of thr pressure
!               qk=(kx*nlx + ky*nly + kz*nlz)/k2
!               ! add the gradient to the non-linear terms
!               nlk1(iz,iy,ix)=nlx - kx*qk
!               nlk2(iz,iy,ix)=nly - ky*qk
!               nlk3(iz,iy,ix)=nlz - kz*qk
!            endif
!         enddo
!      enddo
!   enddo
end subroutine add_grad_pressure


! Compute the nonlinear source term of the mhd equations,
! including penality term, in Fourier space.

! Input: ubk, which is a 4D array containing the Fourier-space version
! of the velocity and magnetic field.

! Output: nlk is the (Fourier-space) nonlinear source term for both
! the velocity and magnetic field. ub is the physical-space version of
! the input velocity and magnetic field.

! Work memory: wj is a 4D array which is used to compute the voriticy
! and current density and other quantities, but doesn't contain anything
! useful at the end of the subroutine.

! Other (global) arrays: mask is 3D physical-space array containing
! the mask function for both the velocity and magnetic field. us is a
! 4D array containing the imposed velocity and magnetic field in
! phsyical space.
subroutine cal_nlk_mhd(nlk,ubk,ub,wj)
  use mpi
  use fsi_vars
  use p3dfft_wrapper
  use basic_operators
  implicit none

  complex(kind=pr),intent(inout) ::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
end subroutine cal_nlk_mhd


! Render the input field divergence-free via a Helmholtz
! decomposition. The zero-mode is left untouched.
subroutine div_field_nul(fx,fy,fz)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr), intent(inout) :: fx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fy(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: ix, iy, iz
  real(kind=pr) :: kx, ky, kz, k2
  complex(kind=pr) :: val, vx,vy,vz

  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3), cb(3)
           kx=wave_x(ix)
           
           k2=kx*kx +ky*ky +kz*kz

           if(k2 /= 0.d0) then
              ! val = (k \cdot{} f) / k^2
              vx=fx(iz,iy,ix)
              vy=fy(iz,iy,ix)
              vz=fz(iz,iy,ix)

              val=(kx*vx + ky*vy + kz*vz)/k2

              ! f <- f - k \cdot{} val
              fx(iz,iy,ix)=vx -kx*val
              fy(iz,iy,ix)=vy -ky*val
              fz(iz,iy,ix)=vz -kz*val
           endif
        enddo
     enddo
  enddo
  
end subroutine div_field_nul



subroutine rhs_acm(u,rhs)
  use mpi 
  use fsi_vars
  use diff
  implicit none
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  
  integer::ix,iy,iz
  real(kind=pr)::ux,uy,uz,vorx,vory,vorz,uxdx,uxdy,uxdz,uydx,uydy,uydz,&
  uzdx,uzdy,uzdz,uxdxdx,uxdydy,uxdzdz,uydxdx,uydydy,uydzdz,uzdxdx,uzdydy,uzdzdz,&
  dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv,pdx,pdy,pdz
  
  call synchronize_ghosts_FD (u)
  
  dxinv = 1.d0/(2.d0*dx)
  dyinv = 1.d0/(2.d0*dy)
  dzinv = 1.d0/(2.d0*dz)
  
  dx2inv = 1.d0/(dx**2)
  dy2inv = 1.d0/(dy**2)
  dz2inv = 1.d0/(dz**2)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ux = u(ix,iy,iz,1)
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)
        
        uxdx = (u(ix+1,iy,iz,1) - u(ix-1,iy,iz,1))*dxinv
        uxdy = (u(ix,iy+1,iz,1) - u(ix,iy-1,iz,1))*dyinv
        uxdz = (u(ix,iy,iz+1,1) - u(ix,iy,iz-1,1))*dzinv
        
        uydx = (u(ix+1,iy,iz,2) - u(ix-1,iy,iz,2))*dxinv
        uydy = (u(ix,iy+1,iz,2) - u(ix,iy-1,iz,2))*dyinv
        uydz = (u(ix,iy,iz+1,2) - u(ix,iy,iz-1,2))*dzinv
        
        uzdx = (u(ix+1,iy,iz,3) - u(ix-1,iy,iz,3))*dxinv
        uzdy = (u(ix,iy+1,iz,3) - u(ix,iy-1,iz,3))*dyinv
        uzdz = (u(ix,iy,iz+1,3) - u(ix,iy,iz-1,3))*dzinv
        
        pdx = (u(ix+1,iy,iz,4) - u(ix-1,iy,iz,4))*dxinv
        pdy = (u(ix,iy+1,iz,4) - u(ix,iy-1,iz,4))*dyinv
        pdz = (u(ix,iy,iz+1,4) - u(ix,iy,iz-1,4))*dzinv        
        
        vorx = uzdy - uydz
        vory = uxdz - uzdx
        vorz = uydx - uxdy
        
        uxdxdx = (u(ix-1,iy,iz,1)-2.d0*u(ix,iy,iz,1)+u(ix+1,iy,iz,1))*dx2inv 
        uxdydy = (u(ix,iy-1,iz,1)-2.d0*u(ix,iy,iz,1)+u(ix,iy+1,iz,1))*dy2inv 
        uxdzdz = (u(ix,iy,iz-1,1)-2.d0*u(ix,iy,iz,1)+u(ix,iy,iz+1,1))*dz2inv 
        
        uydxdx = (u(ix-1,iy,iz,2)-2.d0*u(ix,iy,iz,2)+u(ix+1,iy,iz,2))*dx2inv 
        uydydy = (u(ix,iy-1,iz,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy+1,iz,2))*dy2inv 
        uydzdz = (u(ix,iy,iz-1,2)-2.d0*u(ix,iy,iz,2)+u(ix,iy,iz+1,2))*dz2inv 
        
        uzdxdx = (u(ix-1,iy,iz,3)-2.d0*u(ix,iy,iz,3)+u(ix+1,iy,iz,3))*dx2inv 
        uzdydy = (u(ix,iy-1,iz,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy+1,iz,3))*dy2inv 
        uzdzdz = (u(ix,iy,iz-1,3)-2.d0*u(ix,iy,iz,3)+u(ix,iy,iz+1,3))*dz2inv 
        
        rhs(ix,iy,iz,1) = uy*vorz -uz*vory - pdx + nu*(uxdxdx+uxdydy+uxdzdz)
        rhs(ix,iy,iz,2) = uz*vorx -ux*vorz - pdy + nu*(uydxdx+uydydy+uydzdz)
        rhs(ix,iy,iz,3) = ux*vory -uy*vorx - pdz + nu*(uzdxdx+uzdydy+uzdzdz)
        rhs(ix,iy,iz,4) = -(c_0**2)*(uxdx+uydy+uzdz)
      enddo
    enddo
  enddo     
end subroutine




subroutine rhs_acm_4th(u,rhs)
  use mpi 
  use fsi_vars
  use diff
  implicit none
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  
  integer::ix,iy,iz
  real(kind=pr)::ux,uy,uz,vorx,vory,vorz,uxdx,uxdy,uxdz,uydx,uydy,uydz,&
  uzdx,uzdy,uzdz,uxdxdx,uxdydy,uxdzdz,uydxdx,uydydy,uydzdz,uzdxdx,uzdydy,uzdzdz,&
  dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv,pdx,pdy,pdz,a1,a2,a4,a5,&
  b1,b2,b3,b4,b5
  
  
  call synchronize_ghosts_FD (u)
  
  a1 = 1.d0/12.d0
  a2 =-2.d0/3.d0
  a4 = 2.d0/3.d0
  a5 = -1.d0/12.d0
    
  
  b1=-1.d0/12.d0
  b2=4.d0/3.d0
  b3=-5.d0/2.d0
  b4=4.d0/3.d0
  b5=-1.d0/12.d0
  
  dxinv = 1.d0/dx
  dyinv = 1.d0/dy
  dzinv = 1.d0/dz
  
  dx2inv = 1.d0/(dx**2)
  dy2inv = 1.d0/(dy**2)
  dz2inv = 1.d0/(dz**2)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ux = u(ix,iy,iz,1)
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)
        
        uxdx = (a1*u(ix-2,iy,iz,1)+a2*u(ix-1,iy,iz,1)+a4*u(ix+1,iy,iz,1)+a5*u(ix+2,iy,iz,1))*dxinv        
        uxdy = (a1*u(ix,iy-2,iz,1)+a2*u(ix,iy-1,iz,1)+a4*u(ix,iy+1,iz,1)+a5*u(ix,iy+2,iz,1))*dyinv
        uxdz = (a1*u(ix,iy,iz-2,1)+a2*u(ix,iy,iz-1,1)+a4*u(ix,iy,iz+1,1)+a5*u(ix,iy,iz+2,1))*dzinv
        
        uydx = (a1*u(ix-2,iy,iz,2)+a2*u(ix-1,iy,iz,2)+a4*u(ix+1,iy,iz,2)+a5*u(ix+2,iy,iz,2))*dxinv        
        uydy = (a1*u(ix,iy-2,iz,2)+a2*u(ix,iy-1,iz,2)+a4*u(ix,iy+1,iz,2)+a5*u(ix,iy+2,iz,2))*dyinv
        uydz = (a1*u(ix,iy,iz-2,2)+a2*u(ix,iy,iz-1,2)+a4*u(ix,iy,iz+1,2)+a5*u(ix,iy,iz+2,2))*dzinv
        
        uzdx = (a1*u(ix-2,iy,iz,3)+a2*u(ix-1,iy,iz,3)+a4*u(ix+1,iy,iz,3)+a5*u(ix+2,iy,iz,3))*dxinv        
        uzdy = (a1*u(ix,iy-2,iz,3)+a2*u(ix,iy-1,iz,3)+a4*u(ix,iy+1,iz,3)+a5*u(ix,iy+2,iz,3))*dyinv
        uzdz = (a1*u(ix,iy,iz-2,3)+a2*u(ix,iy,iz-1,3)+a4*u(ix,iy,iz+1,3)+a5*u(ix,iy,iz+2,3))*dzinv
        
        pdx = (a1*u(ix-2,iy,iz,4)+a2*u(ix-1,iy,iz,4)+a4*u(ix+1,iy,iz,4)+a5*u(ix+2,iy,iz,4))*dxinv        
        pdy = (a1*u(ix,iy-2,iz,4)+a2*u(ix,iy-1,iz,4)+a4*u(ix,iy+1,iz,4)+a5*u(ix,iy+2,iz,4))*dyinv
        pdz = (a1*u(ix,iy,iz-2,4)+a2*u(ix,iy,iz-1,4)+a4*u(ix,iy,iz+1,4)+a5*u(ix,iy,iz+2,4))*dzinv       
        
        vorx = uzdy - uydz
        vory = uxdz - uzdx
        vorz = uydx - uxdy
        
        uxdxdx = (b1*u(ix-2,iy,iz,1)+b2*u(ix-1,iy,iz,1)+b3*u(ix,iy,iz,1)+b4*u(ix+1,iy,iz,1)+b5*u(ix+2,iy,iz,1))*dx2inv   
        uxdydy = (b1*u(ix,iy-2,iz,1)+b2*u(ix,iy-1,iz,1)+b3*u(ix,iy,iz,1)+b4*u(ix,iy+1,iz,1)+b5*u(ix,iy+2,iz,1))*dy2inv   
        uxdzdz = (b1*u(ix,iy,iz-2,1)+b2*u(ix,iy,iz-1,1)+b3*u(ix,iy,iz,1)+b4*u(ix,iy,iz+1,1)+b5*u(ix,iy,iz+2,1))*dz2inv   
        
        uydxdx = (b1*u(ix-2,iy,iz,2)+b2*u(ix-1,iy,iz,2)+b3*u(ix,iy,iz,2)+b4*u(ix+1,iy,iz,2)+b5*u(ix+2,iy,iz,2))*dx2inv   
        uydydy = (b1*u(ix,iy-2,iz,2)+b2*u(ix,iy-1,iz,2)+b3*u(ix,iy,iz,2)+b4*u(ix,iy+1,iz,2)+b5*u(ix,iy+2,iz,2))*dy2inv   
        uydzdz = (b1*u(ix,iy,iz-2,2)+b2*u(ix,iy,iz-1,2)+b3*u(ix,iy,iz,2)+b4*u(ix,iy,iz+1,2)+b5*u(ix,iy,iz+2,2))*dz2inv   
        
        uzdxdx = (b1*u(ix-2,iy,iz,3)+b2*u(ix-1,iy,iz,3)+b3*u(ix,iy,iz,3)+b4*u(ix+1,iy,iz,3)+b5*u(ix+2,iy,iz,3))*dx2inv   
        uzdydy = (b1*u(ix,iy-2,iz,3)+b2*u(ix,iy-1,iz,3)+b3*u(ix,iy,iz,3)+b4*u(ix,iy+1,iz,3)+b5*u(ix,iy+2,iz,3))*dy2inv   
        uzdzdz = (b1*u(ix,iy,iz-2,3)+b2*u(ix,iy,iz-1,3)+b3*u(ix,iy,iz,3)+b4*u(ix,iy,iz+1,3)+b5*u(ix,iy,iz+2,3))*dz2inv   
        
        rhs(ix,iy,iz,1) = uy*vorz -uz*vory - pdx + nu*(uxdxdx+uxdydy+uxdzdz)
        rhs(ix,iy,iz,2) = uz*vorx -ux*vorz - pdy + nu*(uydxdx+uydydy+uydzdz)
        rhs(ix,iy,iz,3) = ux*vory -uy*vorx - pdz + nu*(uzdxdx+uzdydy+uzdzdz)
        rhs(ix,iy,iz,4) = -(c_0**2)*(uxdx+uydy+uzdz)
      enddo
    enddo
  enddo     
end subroutine




subroutine rhs_acm_4th2(u,rhs)
  use mpi 
  use fsi_vars
  use diff
  implicit none
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  
  integer::ix,iy,iz
  real(kind=pr)::ux,uy,uz,vorx,vory,vorz,uxdx,uxdy,uxdz,uydx,uydy,uydz,&
  uzdx,uzdy,uzdz,uxdxdx,uxdydy,uxdzdz,uydxdx,uydydy,uydzdz,uzdxdx,uzdydy,uzdzdz,&
  dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv,pdx,pdy,pdz,a1,a2,a4,a5,&
  b1,b2,b3,b4,b5
  
  
  call synchronize_ghosts_FD (u)
  
  a1 = 1.d0/12.d0
  a2 =-2.d0/3.d0
  a4 = 2.d0/3.d0
  a5 = -1.d0/12.d0
    
  
  b1=-1.d0/12.d0
  b2=4.d0/3.d0
  b3=-5.d0/2.d0
  b4=4.d0/3.d0
  b5=-1.d0/12.d0
  
  dxinv = 1.d0/dx
  dyinv = 1.d0/dy
  dzinv = 1.d0/dz
  
  dx2inv = 1.d0/(dx**2)
  dy2inv = 1.d0/(dy**2)
  dz2inv = 1.d0/(dz**2)
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ux = u(ix,iy,iz,1)
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)
        
        uxdy = (a1*u(ix,iy-2,iz,1)+a2*u(ix,iy-1,iz,1)+a4*u(ix,iy+1,iz,1)+a5*u(ix,iy+2,iz,1))*dyinv
        uxdz = (a1*u(ix,iy,iz-2,1)+a2*u(ix,iy,iz-1,1)+a4*u(ix,iy,iz+1,1)+a5*u(ix,iy,iz+2,1))*dzinv
        
        uydx = (a1*u(ix-2,iy,iz,2)+a2*u(ix-1,iy,iz,2)+a4*u(ix+1,iy,iz,2)+a5*u(ix+2,iy,iz,2))*dxinv        
        uydz = (a1*u(ix,iy,iz-2,2)+a2*u(ix,iy,iz-1,2)+a4*u(ix,iy,iz+1,2)+a5*u(ix,iy,iz+2,2))*dzinv
        
        uzdx = (a1*u(ix-2,iy,iz,3)+a2*u(ix-1,iy,iz,3)+a4*u(ix+1,iy,iz,3)+a5*u(ix+2,iy,iz,3))*dxinv        
        uzdy = (a1*u(ix,iy-2,iz,3)+a2*u(ix,iy-1,iz,3)+a4*u(ix,iy+1,iz,3)+a5*u(ix,iy+2,iz,3))*dyinv

        vorx = uzdy - uydz
        vory = uxdz - uzdx
        vorz = uydx - uxdy        
        
        rhs(ix,iy,iz,1) = uy*vorz -uz*vory
        rhs(ix,iy,iz,2) = uz*vorx -ux*vorz
        rhs(ix,iy,iz,3) = ux*vory -uy*vorx
        
      enddo
    enddo
  enddo     
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)        
        uxdx = (a1*u(ix-2,iy,iz,1)+a2*u(ix-1,iy,iz,1)+a4*u(ix+1,iy,iz,1)+a5*u(ix+2,iy,iz,1))*dxinv  
        uydy = (a1*u(ix,iy-2,iz,2)+a2*u(ix,iy-1,iz,2)+a4*u(ix,iy+1,iz,2)+a5*u(ix,iy+2,iz,2))*dyinv
        uzdz = (a1*u(ix,iy,iz-2,3)+a2*u(ix,iy,iz-1,3)+a4*u(ix,iy,iz+1,3)+a5*u(ix,iy,iz+2,3))*dzinv
        rhs(ix,iy,iz,4) = -(c_0**2)*(uxdx+uydy+uzdz)
      enddo
    enddo
  enddo  
  

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)        
        uxdxdx = (b1*u(ix-2,iy,iz,1)+b2*u(ix-1,iy,iz,1)+b3*u(ix,iy,iz,1)+b4*u(ix+1,iy,iz,1)+b5*u(ix+2,iy,iz,1))*dx2inv   
        uxdydy = (b1*u(ix,iy-2,iz,1)+b2*u(ix,iy-1,iz,1)+b3*u(ix,iy,iz,1)+b4*u(ix,iy+1,iz,1)+b5*u(ix,iy+2,iz,1))*dy2inv   
        uxdzdz = (b1*u(ix,iy,iz-2,1)+b2*u(ix,iy,iz-1,1)+b3*u(ix,iy,iz,1)+b4*u(ix,iy,iz+1,1)+b5*u(ix,iy,iz+2,1))*dz2inv   
        
        uydxdx = (b1*u(ix-2,iy,iz,2)+b2*u(ix-1,iy,iz,2)+b3*u(ix,iy,iz,2)+b4*u(ix+1,iy,iz,2)+b5*u(ix+2,iy,iz,2))*dx2inv   
        uydydy = (b1*u(ix,iy-2,iz,2)+b2*u(ix,iy-1,iz,2)+b3*u(ix,iy,iz,2)+b4*u(ix,iy+1,iz,2)+b5*u(ix,iy+2,iz,2))*dy2inv   
        uydzdz = (b1*u(ix,iy,iz-2,2)+b2*u(ix,iy,iz-1,2)+b3*u(ix,iy,iz,2)+b4*u(ix,iy,iz+1,2)+b5*u(ix,iy,iz+2,2))*dz2inv   
        
        uzdxdx = (b1*u(ix-2,iy,iz,3)+b2*u(ix-1,iy,iz,3)+b3*u(ix,iy,iz,3)+b4*u(ix+1,iy,iz,3)+b5*u(ix+2,iy,iz,3))*dx2inv   
        uzdydy = (b1*u(ix,iy-2,iz,3)+b2*u(ix,iy-1,iz,3)+b3*u(ix,iy,iz,3)+b4*u(ix,iy+1,iz,3)+b5*u(ix,iy+2,iz,3))*dy2inv   
        uzdzdz = (b1*u(ix,iy,iz-2,3)+b2*u(ix,iy,iz-1,3)+b3*u(ix,iy,iz,3)+b4*u(ix,iy,iz+1,3)+b5*u(ix,iy,iz+2,3))*dz2inv   
      
        rhs(ix,iy,iz,1) = rhs(ix,iy,iz,1)+nu*(uxdxdx+uxdydy+uxdzdz)
        rhs(ix,iy,iz,2) = rhs(ix,iy,iz,2)+nu*(uydxdx+uydydy+uydzdz)
        rhs(ix,iy,iz,3) = rhs(ix,iy,iz,3)+nu*(uzdxdx+uzdydy+uzdzdz)
      enddo
    enddo
  enddo 
  
  
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)        
        rhs(ix,iy,iz,1) = rhs(ix,iy,iz,1)-(a1*u(ix-2,iy,iz,4)+a2*u(ix-1,iy,iz,4)+a4*u(ix+1,iy,iz,4)+a5*u(ix+2,iy,iz,4))*dxinv        
        rhs(ix,iy,iz,2) = rhs(ix,iy,iz,2)-(a1*u(ix,iy-2,iz,4)+a2*u(ix,iy-1,iz,4)+a4*u(ix,iy+1,iz,4)+a5*u(ix,iy+2,iz,4))*dyinv
        rhs(ix,iy,iz,3) = rhs(ix,iy,iz,3)-(a1*u(ix,iy,iz-2,4)+a2*u(ix,iy,iz-1,4)+a4*u(ix,iy,iz+1,4)+a5*u(ix,iy,iz+2,4))*dzinv       
      enddo
    enddo
  enddo    
  
end subroutine