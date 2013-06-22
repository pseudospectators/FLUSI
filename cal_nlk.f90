! Wrapper for computing the nonlinear source term for Navier-Stokes/MHD
subroutine cal_nlk(time,it,nlk,uk,work_u,work_vort,work)
  use vars
  implicit none

  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  complex(kind=pr),intent(out)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::&
       work_vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent(inout)::work_u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent (in) :: time
  integer, intent(in) :: it

  select case(method(1:3))
  case("fsi") 
     call cal_nlk_fsi(time,it,nlk,uk,work_u,work_vort,work)
  case("mhd") 
     call cal_nlk_mhd(time,it,nlk,uk,work_u,work_vort,work)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in get_params"
     call abort
  end select
end subroutine cal_nlk


! Compute the nonlinear source term of the Navier-Stokes equation,
! including penality term, in Fourier space. Seven real-valued
! arrays are required for working memory.
! FIXME: this does other things as well, like computing energy
! dissipation.
! FIXME: add documentation: which arguments are used for what?
subroutine cal_nlk_fsi(time,it,nlk,uk,work_u,work_vort,work)
  use mpi_header
  use fsi_vars
  implicit none

  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  complex(kind=pr),intent(out)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::&
       work_vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent(inout)::work_u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent (in) :: time
  real(kind=pr) :: t1,t0
  integer, intent(in) :: it
  logical :: TimeForDrag
  integer i

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
  call cofitxyz (uk(:,:,:,1), work_u(:,:,:,1))
  call cofitxyz (uk(:,:,:,2), work_u(:,:,:,2))
  call cofitxyz (uk(:,:,:,3), work_u(:,:,:,3))
  time_u=time_u + MPI_wtime() - t1

  !------------------------------------------------
  ! TEMP: compute divergence
  !-----------------------------------------------
  ! if (TimeForDrag) call compute_divergence(FIXME)
  
  !-----------------------------------------------
  !-- Compute vorticity
  !-----------------------------------------------
  t1=MPI_wtime()
  !nlk is temporarily used for vortk
  call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
       uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3)) 
  do i=1,3
     call cofitxyz(work_vort(:,:,:,i),nlk(:,:,:,i))
  enddo
  ! timing statistics
  time_vor=time_vor + MPI_wtime() - t1

  !-----------------------------------------------
  !-- Compute kinetic energy and dissipation rate + mask volume
  !-----------------------------------------------  
  if((TimeForDrag) .and. (iKinDiss==1)) then
     call Energy_Dissipation (work_u, work_vort )
  endif

  !-------------------------------------------------------------
  !-- Calculate omega x u (cross-product)
  !-- add penalization term
  !-- and transform the result into Fourier space 
  !-------------------------------------------------------------
  t1=MPI_wtime()
  if (iPenalization==1) then
     call omegacrossu_penalize(nlk,work,work_u,work_vort,TimeForDrag)
  else ! no penalization
     call omegacrossu_nopen(nlk,work,work_u,work_vort)
  endif
  ! timing statistics
  time_curl=time_curl + MPI_wtime() - t1

  t1=MPI_wtime()
  call add_grad_pressure(nlk)
  time_p=time_p + MPI_wtime() - t1

  ! this is for the timing statistics.
  ! how much time was spend on ffts in cal_nlk?
  time_nlk_fft=time_nlk_fft + time_fft2 + time_ifft2
  ! how much time was spend on cal_nlk
  time_nlk=time_nlk + MPI_wtime() - t0
end subroutine cal_nlk_fsi


! This subroutine takes one component of the penalization term (work)
! computes the integral over it, which is the hydrodynamic force in
! the direction iDirection. The force is stored in the GlobIntegrals
! structure
subroutine IntegralForce (work, iDirection ) 
  use fsi_vars
  use mpi_header
  implicit none
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
       intent (in) :: work
  integer, intent (in) :: iDirection
  integer :: mpicode
  real (kind=pr) :: Force_local

  Force_local =dx*dy*dz*sum( work )  
  call MPI_REDUCE (Force_local,GlobIntegrals%Force(iDirection),1,mpireal,&
       MPI_SUM,0,MPI_COMM_WORLD,mpicode)  ! max at 0th process  
end subroutine IntegralForce


! Compute the kinetic energy, dissipation rate and mask volume.  Store
! all these in the structure GlobIntegrals (definition see
! fsi_vars).
subroutine Energy_Dissipation (u, vort )
  use fsi_vars
  use mpi_header
  implicit none

  real (kind=pr),intent (in) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real (kind=pr),intent (in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real (kind=pr) :: E_kin_local, Dissip_local, Volume_local
  integer :: mpicode

  E_kin_local =0.5d0*dx*dy*dz*sum( u*u )
  Dissip_local=-nu*dx*dy*dz*sum( vort*vort )

  if (iPenalization==1) then
     Volume_local=dx*dy*dz*sum( mask )*eps
  else
     Volume_local=0.d0
  endif

  call MPI_REDUCE(E_kin_local,GlobIntegrals%E_kin,1,mpireal,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)  ! sum at 0th process
  call MPI_REDUCE(Dissip_local,GlobIntegrals%Dissip,1,mpireal,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)  ! sum at 0th process
  call MPI_REDUCE(Volume_local,GlobIntegrals%Volume,1,mpireal,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)  ! sum at 0th process
end subroutine Energy_Dissipation


! Compute the pressure. It is given by the divergence of the non-linear
! terms (nlk: intent(in)) divided by k**2.
! so: p=(i*kx*sxk + i*ky*syk + i*kz*szk) / k**2 
! note: we use rotational formulation: p is NOT the physical pressure

! FIXME: I think there is the imaginary unit missing!! 

subroutine compute_pressure(pk,nlk)
  use mpi_header
  use vars
  implicit none

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr),intent(out):: pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  do iy=ca(3),cb(3)  ! ky : 0..ny/2-1 ,then, -ny/2..-1
     ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
     do ix=ca(2),cb(2) ! kx : 0..nx/2
        kx=scalex*dble(ix)
        do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1
           kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           k2=kx*kx+ky*ky+kz*kz
           if(k2 .ne. 0.0) then
              ! contains the pressure in Fourier space
              pk(iz,ix,iy)=(kx*nlk(iz,ix,iy,1)&
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
subroutine add_grad_pressure(nlk)
  use mpi_header
  use vars
  implicit none

  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr) :: qk
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)

  do iy=ca(3),cb(3) ! ky : 0..ny/2-1 ,then, -ny/2..-1     
     ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
     do ix=ca(2),cb(2) ! kx : 0..nx/2
        kx=scalex*dble(ix)
        do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1         
           kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           k2=kx*kx+ky*ky+kz*kz
           if (k2 .ne. 0.0) then  
              ! qk is the Fourier coefficient of thr pressure
              qk=(kx*nlk(iz,ix,iy,1)+ky*nlk(iz,ix,iy,2)+kz*nlk(iz,ix,iy,3))/k2
              ! add the gradient to the non-linear terms
              nlk(iz,ix,iy,1)=nlk(iz,ix,iy,1) - kx*qk  
              nlk(iz,ix,iy,2)=nlk(iz,ix,iy,2) - ky*qk
              nlk(iz,ix,iy,3)=nlk(iz,ix,iy,3) - kz*qk
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
!!$  call cofitxyz(nlk(:,:,:,1),work)
!!$  call MPI_REDUCE (maxval(work), GlobIntegrals%Divergence, 1, mpireal, &
!!$       MPI_MAX, 0, MPI_COMM_WORLD, mpicode)  ! max at 0th process  
!!$
!!$  call Energy_Dissipation ( GlobIntegrals, work_u, work_vort )
!!$
!!$  if (mpirank ==0) then
!!$     write (*,'("max{div(u)}=",es15.8,"max{div(u)}/||u||=",es15.8)') &
!!$          GlobIntegrals%Divergence, GlobIntegrals%Divergence/(2.d0*GlobIntegrals%E_kin)
!!$  endif
end subroutine compute_divergence


! Compute non-linear transport term (omega cross u) and transform it to Fourier
! space. This is the case without penalization.
subroutine omegacrossu_nopen(nlk,work,u,vort)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  
  work=u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2)
  call coftxyz(work,nlk(:,:,:,1))
  work=u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3)
  call coftxyz(work,nlk(:,:,:,2))
  work=u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1)
  call coftxyz(work,nlk(:,:,:,3))
endsubroutine omegacrossu_nopen


! Compute non-linear transport term (omega cross u) and transform it to Fourier
! space. This is the case with penalization. Therefore we compute the penalty
! term mask*(u-us) as well. This gives the occasion to compute the drag forces,
! if it is time to do so (TimeForDrag=.true.). The drag is returned in 
! GlobIntegrals.
subroutine omegacrossu_penalize(nlk,work,u,vort,TimeForDrag)
  use mpi_header
  use fsi_vars
  implicit none

  real(kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  logical,intent(in) :: TimeForDrag
  
  ! x component
  call Penalize ( work, u, 1, TimeForDrag)
  work=u(:,:,:,2)*vort(:,:,:,3) - u(:,:,:,3)*vort(:,:,:,2) + work
  call coftxyz (work, nlk(:,:,:,1))
  
  ! y component
  call Penalize ( work, u, 2, TimeForDrag)
  work=u(:,:,:,3)*vort(:,:,:,1) - u(:,:,:,1)*vort(:,:,:,3) + work
  call coftxyz (work, nlk(:,:,:,2))
  
  ! z component
  call Penalize ( work, u, 3, TimeForDrag)
  work=u(:,:,:,1)*vort(:,:,:,2) - u(:,:,:,2)*vort(:,:,:,1) + work
  call coftxyz (work, nlk(:,:,:,3))     
end subroutine omegacrossu_penalize


! we outsource the actual penalization (even though its a fairly
! simple process) to remove some lines in the actual cal_nlk also, at
! this occasion, we directly compute the integral hydrodynamic forces,
! if its time to do so (TimeForDrag=.true.)
subroutine Penalize(work,work_u,iDir,TimeForDrag)    
  use fsi_vars
  implicit none

  integer, intent (in) :: iDir
  logical, intent (in) :: TimeForDrag
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)),&
       intent (out) :: work
  real (kind=pr), dimension (ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf),&
       intent (in) :: work_u


  ! compute penalization term
  if (iMoving == 1) then
     work=-mask*(work_u(:,:,:,iDir) - us(:,:,:,iDir))  
  else
     work=-mask*(work_u(:,:,:,iDir))
  endif

  ! if its time, compute drag forces
  if ((TimeForDrag).and.(iDrag==1)) then
     call IntegralForce (work, iDir ) 
  endif
end subroutine Penalize



! Compute the nonlinear source term of the mhd equations,
! including penality term, in Fourier space. 
! FIXME: add documentation: which arguments are used for what?
subroutine cal_nlk_mhd(time,it,nlk,uk,work_u,work_vort,work)
  use mpi_header
  use fsi_vars
  implicit none

  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  complex(kind=pr),intent(out)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::&
       work_vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent(inout)::work_u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nf)
  real(kind=pr),intent (in) :: time
  real(kind=pr) :: t1,t0
  integer, intent(in) :: it
  logical :: TimeForDrag

  
  ! FIXME: actually write the nonlinear term here, yo.
  
  ! compute u and B in physical space

  ! compute vorticity and current density

  ! u x omega + j x B

!  call curl(vort(:,:,:,1:3),vortk(:,:,:,1:3),uk(:,:,:,1:3))
!  call curl(vort(:,:,:,1),vortk(:,:,:,1),uk(:,:,:,1))

  ! ik x (u x B)
  
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
  
  ! comput vorticity in Fourier space:
  do iy=ca(3),cb(3)    ! ky : 0..ny/2-1 ,then,-ny/2..-1
     ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)
     do ix=ca(2),cb(2)  ! kx : 0..nx/2
        kx=scalex*dble(ix)
        do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then,-nz/2..-1
           kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
           out1(iz,ix,iy)=imag*(ky*in3(iz,ix,iy)-kz*in2(iz,ix,iy))
           out2(iz,ix,iy)=imag*(kz*in1(iz,ix,iy)-kx*in3(iz,ix,iy))
           out3(iz,ix,iy)=imag*(kx*in2(iz,ix,iy)-ky*in1(iz,ix,iy))
        enddo
     enddo
  enddo

  ! Transform to physical space
end subroutine curl
