! Wrapper for computing the nonlinear source term for Navier-Stokes/MHD
subroutine cal_nlk(time,it,nlk,uk,u,vort,work)
  use vars
  implicit none

  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(in) :: time
  integer, intent(in) :: it

  select case(method(1:3))
  case("fsi") 
     call cal_nlk_fsi(time,it,nlk,uk,u,vort,work)
  case("mhd") 
     call cal_nlk_mhd(nlk,uk,u,vort)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in cal_nlk"
     call abort
  end select
end subroutine cal_nlk


! Compute the nonlinear source term of the Navier-Stokes equation,
! including penality term, in Fourier space. Seven real-valued
! arrays are required for working memory.
! FIXME: this does other things as well, like computing energy
! dissipation.
! FIXME: add documentation: which arguments are used for what?
subroutine cal_nlk_fsi(time,it,nlk,uk,u,vort,work)
  use mpi_header
  use fsi_vars
  implicit none

  complex(kind=pr),intent(in):: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent (in) :: time
  real(kind=pr) :: t1,t0
  integer, intent(in) :: it
  logical :: TimeForDrag ! FIXME: move to time_step routine?
  integer i

  ! FIXME: move into write_integrals_fsi
  ! is it time for save global quantities?
  TimeForDrag=.false.
  ! note we do this every itdrag time steps
  if (modulo(it,itdrag)==0) TimeForDrag=.true. ! yes, indeed
  
  ! performance measurement in global variables
  t0=MPI_wtime()
  time_fft2 =0.0 ! time_fft2 is the time spend on ffts during cal_nlk only
  time_ifft2=0.0 ! time_ifft2 is the time spend on iffts during cal_nlk only

  !-----------------------------------------------
  !-- Calculate ux and uy in physical space
  !-----------------------------------------------
  t1=MPI_wtime()
  do i=1,nd
     call ifft(u(:,:,:,i),uk(:,:,:,i))
  enddo
  time_u=time_u + MPI_wtime() - t1

  ! FIXME: move into write_integrals_fsi
  !------------------------------------------------
  ! TEMP: compute divergence
  !-----------------------------------------------
  ! if (TimeForDrag) call compute_divergence(FIXME)
  
  !-----------------------------------------------
  !-- Compute vorticity
  !-----------------------------------------------
  t1=MPI_wtime()
  ! NB: nlk is temporarily used for vortk
  call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
             uk(:,:,:,1), uk(:,:,:,2), uk(:,:,:,3)) 

  do i=1,3
     call ifft(vort(:,:,:,i),nlk(:,:,:,i))
  enddo

  ! Timing statistics
  time_vor=time_vor + MPI_wtime() - t1

  ! FIXME: move into write_integrals_fsi
  !-----------------------------------------------
  !-- Compute kinetic energy and dissipation rate + mask volume
  !-----------------------------------------------  
  if((TimeForDrag) .and. (iKinDiss==1)) then
     call Energy_Dissipation (u,vort)
  endif
  
  !-------------------------------------------------------------
  !-- Calculate omega x u (cross-product)
  !-- add penalization term
  !-- and transform the result into Fourier space 
  !-------------------------------------------------------------
  t1=MPI_wtime()
  if (iPenalization==1) then
     call omegacrossu_penalize(work,u,vort,TimeForDrag,nlk)
  else ! no penalization
     call omegacrossu_nopen(work,u,vort,nlk)
  endif
  ! timing statistics
  time_curl=time_curl + MPI_wtime() - t1

  t1=MPI_wtime()
  call add_grad_pressure(nlk(:,:,:,1),nlk(:,:,:,3),nlk(:,:,:,3))
  time_p=time_p + MPI_wtime() - t1

  ! this is for the timing statistics.
  ! how much time was spend on ffts in cal_nlk?
  time_nlk_fft=time_nlk_fft + time_fft2 + time_ifft2
  ! how much time was spend on cal_nlk
  time_nlk=time_nlk + MPI_wtime() - t0
end subroutine cal_nlk_fsi


! This subroutine takes one component of the penalization term (work)
! computes the integral over it, which is the hydrodynamic force in
! the direction iDirection. The force is stored in the GlobalIntegrals
! structure
subroutine IntegralForce(work,iDirection) 
  use fsi_vars
  use mpi_header
  implicit none
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
       intent (in) :: work
  integer, intent (in) :: iDirection
  integer :: mpicode
  real (kind=pr) :: Force_local

  Force_local =dx*dy*dz*sum( work )  
  call MPI_REDUCE (Force_local,GlobalIntegrals%Force(iDirection),1,mpireal,&
       MPI_SUM,0,MPI_COMM_WORLD,mpicode)  ! max at 0th process  
end subroutine IntegralForce


! Compute the kinetic energy, dissipation rate and mask volume.  Store
! all these in the structure GlobalIntegrals (definition see
! fsi_vars).
subroutine Energy_Dissipation(u,vort)
  use fsi_vars
  use mpi_header
  implicit none

  real (kind=pr),intent (in) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent (in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: Ekin_local, Dissip_local, Volume_local
  integer :: mpicode

  Ekin_local =0.5d0*dx*dy*dz*sum( u*u )
  Dissip_local=-nu*dx*dy*dz*sum( vort*vort )

  if (iPenalization==1) then
     Volume_local=dx*dy*dz*sum( mask )*eps
  else
     Volume_local=0.d0
  endif

   ! sum at 0th process
  call MPI_REDUCE(Ekin_local,GlobalIntegrals%Ekin,1,mpireal,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode) 
  call MPI_REDUCE(Dissip_local,GlobalIntegrals%Dissip,1,mpireal,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(Volume_local,GlobalIntegrals%Volume,1,mpireal,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)
end subroutine Energy_Dissipation


! Compute the pressure. It is given by the divergence of the non-linear
! terms (nlk: intent(in)) divided by k**2.
! so: p=(i*kx*sxk + i*ky*syk + i*kz*szk) / k**2 
! note: we use rotational formulation: p is NOT the physical pressure
subroutine compute_pressure(pk,nlk)
  use mpi_header
  use vars
  implicit none

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr),intent(out):: pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1
     kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
     do ix=ca(2),cb(2) ! kx : 0..nx/2
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)  ! ky : 0..ny/2-1 ,then, -ny/2..-1
           ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)

           k2=kx*kx + ky*ky + kz*kz
           if(k2 .ne. 0.0) then
              ! contains the pressure in Fourier space
              pk(iz,ix,iy)=imag*(&
                   kx*nlk(iz,ix,iy,1)&
                   +ky*nlk(iz,ix,iy,2)&
                   +kz*nlk(iz,ix,iy,3)&
                   )/k2
           endif
        enddo
     enddo
  enddo
end subroutine compute_pressure


! Add the gradient of the pressure to the nonlinear term, which is the actual
! projection scheme used in this code. The non-linear term comes in with NL and
! penalization and leaves divergence free
subroutine add_grad_pressure(nlk1,nlk2,nlk3)
  use mpi_header
  use vars
  implicit none

  complex(kind=pr),intent(inout):: nlk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: nlk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: nlk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr) :: qk, nlx,nly,nlz
  complex(kind=pr) :: imag   ! imaginary unit

  imag = dcmplx(0.d0,1.d0)
  
  do iz=ca(1),cb(1)
     kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
           
           k2=kx*kx + ky*ky + kz*kz

           if (k2 .ne. 0.0) then
              nlx=nlk1(iz,ix,iy)
              nly=nlk2(iz,ix,iy)
              nlz=nlk3(iz,ix,iy)

              ! qk is the Fourier coefficient of thr pressure
              qk=(kx*nlx + ky*nly + kz*nlz)/k2
              ! add the gradient to the non-linear terms
              nlk1(iz,ix,iy)=nlx - kx*qk
              nlk2(iz,ix,iy)=nly - ky*qk
              nlk3(iz,ix,iy)=nlz - kz*qk
           endif
        enddo
     enddo
  enddo
end subroutine add_grad_pressure


!FIXME: temp code, removed from main cal_nlk
subroutine compute_divergence()
  use mpi_header
  use vars
  implicit none
!!$  
!!$  ! compute max val of {|div(.)|/|.|} over entire domain
!!$  do iz=ca(1),cb(1)
!!$     kz=scalez*(modulo(iz+nz/2,nz) -nz/2)
!!$     do iy=ca(3),cb(3)
!!$        ky=scaley*(modulo(iy+ny/2,ny) -ny/2)
!!$        do ix=ca(2),cb(2)
!!$           kx=scalex*ix
!!$           ! divergence of velocity field
!!$           nlk(iz,ix,iy,1)=dcmplx(0.d0,1.d0)*(kx*uk(iz,ix,iy,1)+ky*uk(iz,ix,iy,2)+kz*uk(iz,ix,iy,3))
!!$        enddo
!!$     enddo
!!$  enddo
!!$  ! now nlk(:,:,:,1) contains divergence field
!!$  call ifft(nlk(:,:,:,1),work)
!!$  call MPI_REDUCE (maxval(work), GlobalIntegrals%Divergence, 1, mpireal, &
!!$       MPI_MAX, 0, MPI_COMM_WORLD, mpicode)  ! max at 0th process  
!!$
!!$  call Energy_Dissipation ( GlobalIntegrals, u, vort )
!!$
!!$  if (mpirank ==0) then
!!$     write (*,'("max{div(u)}=",es15.8,"max{div(u)}/||u||=",es15.8)') &
!!$          GlobalIntegrals%Divergence, GlobalIntegrals%Divergence/(2.d0*GlobalIntegrals%Ekin)
!!$  endif
end subroutine compute_divergence


! Compute non-linear transport term (omega cross u) and transform it to Fourier
! space. This is the case without penalization.
subroutine omegacrossu_nopen(work,u,vort,nlk)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(out) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  work=u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2)
  call fft(nlk(:,:,:,1),work)
  work=u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3)
  call fft(nlk(:,:,:,2),work)
  work=u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1)
  call fft(nlk(:,:,:,3),work)
endsubroutine omegacrossu_nopen


! Compute non-linear transport term (omega cross u) and transform it
! to Fourier space. This is the case with penalization. Therefore we
! compute the penalty term mask*(u-us) as well. This gives the
! occasion to compute the drag forces, if it is time to do so
! (TimeForDrag=.true.). The drag is returned in GlobalIntegrals.
subroutine omegacrossu_penalize(work,u,vort,TimeForDrag,nlk)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  logical,intent(in) :: TimeForDrag
  
  ! x component
  call Penalize(work,u,1,TimeForDrag)
  ! FIXME: move into write_integrals_fsi?
  ! if its time, compute drag forces
  if ((TimeForDrag).and.(iDrag==1)) then
     call IntegralForce (work,1) 
  endif
  work=work + u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2)
  call fft(nlk(:,:,:,1),work)
  
  ! y component
  call Penalize(work,u,2,TimeForDrag)
  ! FIXME: move into write_integrals_fsi?
  ! if its time, compute drag forces
  if ((TimeForDrag).and.(iDrag==1)) then
     call IntegralForce (work,2) 
  endif
  work=work + u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3)
  call fft(nlk(:,:,:,2),work)
  
  ! z component
  call Penalize(work,u,3,TimeForDrag)
  ! FIXME: move into write_integrals_fsi?
  ! if its time, compute drag forces
  if ((TimeForDrag).and.(iDrag==1)) then
     call IntegralForce (work,3) 
  endif
  work=work + u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1)
  call fft(nlk(:,:,:,3),work)
end subroutine omegacrossu_penalize


! The actual penalization (even though its a fairly simple process) is
! done by this routine instead of in the computation of the nonlinear
! term in order to remove some lines in the actual cal_nlk also.
subroutine Penalize(work,u,iDir,TimeForDrag)
  use fsi_vars
  implicit none

  integer, intent(in) :: iDir
  logical, intent(in) :: TimeForDrag
  real (kind=pr), intent(out) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr), intent(in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  ! compute penalization term
  if (iMoving == 1) then
     work=-mask*(u(:,:,:,iDir) - us(:,:,:,iDir))  
  else
     work=-mask*(u(:,:,:,iDir))
  endif
end subroutine Penalize


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
  use mpi_header
  use fsi_vars
  implicit none

  complex(kind=pr),intent(inout) ::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
   !  real(kind=pr) :: t1,t0
  integer :: i,ix,iy,iz
  real(kind=pr) :: w1,w2,w3,j1,j2,j3
  real(kind=pr) :: u1,u2,u3,b1,b2,b3
  real(kind=pr) :: m,us1,us2,us3

  ! Transform u and B into physical space:
  do i=1,nd
     call ifft(ub(:,:,:,i),ubk(:,:,:,i))
  enddo

  ! Compute us, the imposed penalty field:
  call update_us(ub) ! TODO: only update when necessary

  ! Compute the vorticity and store the result in the first three 3D
  ! arrays of nlk.
  call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
       ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3))

  ! Compute the current density and store the result in the last three
  ! 3D arrays of nlk.
  call curl(nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6),&
       ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6))

  ! Transform vorcitity and current density to physical space, store
  ! in wj
  do i=1,nd
     call ifft(wj(:,:,:,i),nlk(:,:,:,i))
  enddo

  ! Put the x-space version of the nonlinear source term in wj.
  do ix=ra(1),rb(1)
     do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
           ! Loop-local variables for velocity and magnetic field:
           u1=ub(ix,iy,iz,1)
           u2=ub(ix,iy,iz,2)
           u3=ub(ix,iy,iz,3)
           b1=ub(ix,iy,iz,4)
           b2=ub(ix,iy,iz,5)
           b3=ub(ix,iy,iz,6)

           ! Loop-local variables for vorticity and current density:
           w1=wj(ix,iy,iz,1)
           w2=wj(ix,iy,iz,2)
           w3=wj(ix,iy,iz,3)
           j1=wj(ix,iy,iz,4)
           j2=wj(ix,iy,iz,5)
           j3=wj(ix,iy,iz,6)

           ! Loop-local variables for mask and imposed velocity field:
           m=mask(ix,iy,iz)
           us1=us(ix,iy,iz,1)
           us2=us(ix,iy,iz,2)
           us3=us(ix,iy,iz,3)

            ! Nonlinear source term for fluid, including penalization:
            wj(ix,iy,iz,1)=u2*w3 - u3*w2 + j2*b3 - j3*b2 -m*(u1-us1)
            wj(ix,iy,iz,2)=u3*w1 - u1*w3 + j3*b1 - j1*b3 -m*(u2-us2)
            wj(ix,iy,iz,3)=u1*w2 - u2*w1 + j1*b2 - j2*b1 -m*(u3-us3)

            ! Nonlinear source term for magnetic field (missing the
            ! curl and without penalization):
            wj(ix,iy,iz,4)=u2*b3 - u3*b2
            wj(ix,iy,iz,5)=u3*b1 - u1*b3
            wj(ix,iy,iz,6)=u1*b2 - u2*b1
        enddo
     enddo
  enddo

  ! Transform B to Fourier space.  Keep the first three fields free so
  ! that we can use it to store the penalization for the B field.
  do i=4,nd
     call fft(nlk(:,:,:,i),wj(:,:,:,i))
  enddo
  ! NB: the last three sub-arrays of wj and the first three sub-arrays
  ! of nlk are free.

  ! Add the curl to the magnetic source term:
  call curl_inplace(nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6))

  ! Penalization for B-field:
  if(iPenalization == 1) then
     do i=4,nd
        wj(:,:,:,4)=-mask*(ub(:,:,:,i) - us(:,:,:,i))
        call fft(nlk(:,:,:,1),wj(:,:,:,4))
        nlk(:,:,:,i)=nlk(:,:,:,i) + nlk(:,:,:,1)
     enddo
  endif
  
  ! Transform u source-term to Fourier space:
  do i=1,3
     call fft(nlk(:,:,:,i),wj(:,:,:,i))
  enddo

  ! NB: wj is now completely free, and contains nothing useful.

  ! Add the gradient of the pseudo-pressure to the source term of the
  ! fluid.
  call add_grad_pressure(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3))

  ! Make the source term for the magnetic field divergence-free via a
  ! Helmholtz decomposition.
  call div_field_nul(nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6))
end subroutine cal_nlk_mhd


! Given three components of an input fields in Fourier space, compute
! the curl in physical space.  Arrays are 3-dimensional.
subroutine curl(out1,out2,out3,in1,in2,in3)
  use mpi_header
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
     kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)

           out1(iz,ix,iy)=imag*(ky*in3(iz,ix,iy) -kz*in2(iz,ix,iy))
           out2(iz,ix,iy)=imag*(kz*in1(iz,ix,iy) -kx*in3(iz,ix,iy))
           out3(iz,ix,iy)=imag*(kx*in2(iz,ix,iy) -ky*in1(iz,ix,iy))
        enddo
     enddo
  enddo
end subroutine curl


! Given three components of a fields in Fourier space, compute the
! curl in physical space.  Arrays are 3-dimensional.
subroutine curl_inplace(fx,fy,fz)
  use mpi_header
  use vars
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
     kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
           
           t1=fx(iz,ix,iy)
           t2=fy(iz,ix,iy)
           t3=fz(iz,ix,iy)

           fx(iz,ix,iy)=imag*(ky*t3 - kz*t2)
           fy(iz,ix,iy)=imag*(kz*t1 - kx*t3)
           fz(iz,ix,iy)=imag*(kx*t2 - ky*t1)
        enddo
     enddo
  enddo
end subroutine curl_inplace


! Render the input field divergence-free via a Helmholtz
! decomposition. The zero-mode is left untouched.
subroutine div_field_nul(fx,fy,fz)
  use mpi_header
  use vars
  implicit none

  complex(kind=pr), intent(inout) :: fx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fy(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  integer :: ix, iy, iz
  real(kind=pr) :: kx, ky, kz, k2
  complex(kind=pr) :: val, vx,vy,vz

  do iz=ca(1),cb(1)
     kz=scalez*dble(modulo(iz+nz/2,nz) -nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3), cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny) -ny/2)
           
           k2=kx*kx +ky*ky +kz*kz

           if(k2 /= 0.d0) then
              ! val = (k \cdot{} f) / k^2
              vx=fx(iz,ix,iy)
              vy=fy(iz,ix,iy)
              vz=fz(iz,ix,iy)

              val=(kx*vx + ky*vy + kz*vz)/k2

              ! f <- f - k \cdot{} val
              fx(iz,ix,iy)=vx -kx*val
              fy(iz,ix,iy)=vy -ky*val
              fz(iz,ix,iy)=vz -kz*val
           endif
        enddo
     enddo
  enddo
  
end subroutine div_field_nul
