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


!-------------------------------------------------------------------------------
! Compute the nonlinear source term of the Navier-Stokes equation,
! including penalty term, in Fourier space. Seven real-valued
! arrays are required for working memory.
! Input:
!       time: guess what!
!       uk: 4D field of velocity in Fourier space
! Output:
!       nlk:  The right hand side of penalized Navier-Stokes in Fourier space
!       vort: work array, contains NL+penal term in phys space
!             this is reused in save_fields to compute the pressure snapshot
!       u:    work array, contains the velocity in phys space
!             this is reused in the caller FluidTimestep to adjust dt
!       work: work array, currently unused (will be used for pressure)
! Side Effects:
!      * computes the hydrodynamic forces+torques at every time step and stores 
!        them in a global struct on root rank. we will need this for free flight
!      * if present, a sponge is applied to remove incoming vorticity
!
! To Do:
!       * for "true" FSI, we'll need to return the pressure field in phys space
!-------------------------------------------------------------------------------
subroutine cal_nlk_fsi(time,it,nlk,uk,u,vort,work)
  use mpi_header
  use fsi_vars
  implicit none

  complex(kind=pr),intent(in)::   uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(out):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout):: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent (in) :: time
  real(kind=pr) :: t1,t0,ux,uy,uz,vorx,vory,vorz,chi,usx,usy,usz
  real(kind=pr) :: penalx,penaly,penalz,forcex,forcey,forcez,xlev,ylev,zlev
  real(kind=pr) :: torquex,torquey,torquez
  integer, intent(in) :: it
  integer :: ix,iz,iy,mpicode, i


  
  ! performance measurement in global variables
  t0         = MPI_wtime()
  time_fft2  = 0.d0 ! time_fft2 is the time spend on ffts during cal_nlk only
  time_ifft2 = 0.d0 ! time_ifft2 is the time spend on iffts during cal_nlk only
  penalx     = 0.d0
  penaly     = 0.d0
  penalz     = 0.d0
  forcex     = 0.d0
  forcey     = 0.d0
  forcez     = 0.d0  
  torquex    = 0.d0
  torquey    = 0.d0
  torquez    = 0.d0    
  usx     = 0.d0
  usy     = 0.d0
  usz     = 0.d0  

  !-----------------------------------------------
  !-- Calculate velocity in physical space
  !-----------------------------------------------
  t1 = MPI_wtime()
  call ifft3 (u, uk)
  time_u = time_u + MPI_wtime() - t1
  
  !-----------------------------------------------
  !-- Compute vorticity
  !-----------------------------------------------
  t1 = MPI_wtime()
  ! nlk is temporarily used for vortk
  call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
             uk(:,:,:,1), uk(:,:,:,2), uk(:,:,:,3)) 
  ! transform it to physical space
  call ifft3 (vort, nlk)  
  time_vor = time_vor + MPI_wtime() - t1
  
  !-----------------------------------------------
  !-- vorticity sponge term
  !-----------------------------------------------
  t1 = MPI_wtime()
  if (iVorticitySponge == "yes") then  
    ! loop over components
    do i=1,3  
      ! penalize vorticity in the work array
      call penalize_vort ( work, vort(:,:,:,i) )
      ! then transform it to fourier space and store it in the global
      ! array sponge
      call fft ( sponge(:,:,:,i), work )  
    enddo  
    ! transform vorticity sponge term to velocity
    call vort_sponge()  
  endif
  time_sponge = time_sponge + MPI_wtime() - t1
  
  
  !-------------------------------------------------------------
  !-- Non-Linear terms
  !-------------------------------------------------------------
  t1 = MPI_wtime()
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
        
        ! local variables for penalization
        if (iPenalization==1) then
          ! for torque moment
          xlev = dble(ix)*dx - x0
          ylev = dble(iy)*dy - y0
          zlev = dble(iz)*dz - z0

          chi  = mask(ix,iy,iz)
          if (iMoving==1) then
            usx  = us(ix,iy,iz,1)
            usy  = us(ix,iy,iz,2)
            usz  = us(ix,iy,iz,3)
          endif
          ! actual penalization term
          penalx = -chi*(ux-usx)
          penaly = -chi*(uy-usy)
          penalz = -chi*(uz-usz)
          
          ! integrate forces + torques (note sign inversion!)
          forcex = forcex - penalx
          forcey = forcey - penaly
          forcez = forcez - penalz              
          torquex = torquex - (ylev*penalz - zlev*penaly)
          torquey = torquey - (zlev*penalx - xlev*penalz)
          torquez = torquez - (xlev*penaly - ylev*penalx)
        endif
        
        ! we overwrite the vorticity with the NL terms in phys space
        ! note this is indeed -(vor x u) (negative sign)
        vort(ix,iy,iz,1) = uy*vorz - uz*vory + penalx
        vort(ix,iy,iz,2) = uz*vorx - ux*vorz + penaly
        vort(ix,iy,iz,3) = ux*vory - uy*vorx + penalz
      enddo
    enddo
  enddo
  ! to Fourier space
  call fft3(nlk,vort)  
  time_curl = time_curl + MPI_wtime() - t1
  
  !-----------------------------------------------
  ! add sponge term
  !-----------------------------------------------
  t1 = MPI_wtime()
  if (iVorticitySponge == "yes") then
    nlk(:,:,:,1) = nlk(:,:,:,1) + sponge(:,:,:,1)
    nlk(:,:,:,2) = nlk(:,:,:,2) + sponge(:,:,:,2)
    nlk(:,:,:,3) = nlk(:,:,:,3) + sponge(:,:,:,3)
  endif
  time_sponge = time_sponge + MPI_wtime() - t1
  
  !-----------------------------------------------
  !-- add pressure gradient
  !-----------------------------------------------  
  t1 = MPI_wtime()
  call add_grad_pressure(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3))
  time_p = time_p + MPI_wtime() - t1

  !-----------------------------------------------
  !-- done. nlk is the right hand side of naiver-stokes in Fourier space
  !-----------------------------------------------  
  
  ! save global forces
  forcex = forcex*dx*dy*dz
  forcey = forcey*dx*dy*dz
  forcez = forcez*dx*dy*dz  
  call MPI_REDUCE (forcex,GlobalIntegrals%Force(1),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)  
  call MPI_REDUCE (forcey,GlobalIntegrals%Force(2),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode) 
  call MPI_REDUCE (forcez,GlobalIntegrals%Force(3),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)    
                   
  torquex = torquex*dx*dy*dz
  torquey = torquey*dx*dy*dz
  torquez = torquez*dx*dy*dz  
  call MPI_REDUCE (torquex,GlobalIntegrals%Torque(1),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)  
  call MPI_REDUCE (torquey,GlobalIntegrals%Torque(2),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode) 
  call MPI_REDUCE (torquez,GlobalIntegrals%Torque(3),1,mpireal,MPI_SUM,0,&
                   MPI_COMM_WORLD,mpicode)                    
                   
  ! this is for the timing statistics.
  ! how much time was spend on ffts in cal_nlk?
  time_nlk_fft=time_nlk_fft + time_fft2 + time_ifft2
  ! how much time was spend on cal_nlk
  time_nlk=time_nlk + MPI_wtime() - t0
end subroutine cal_nlk_fsi



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

  ! important remark:
  ! as the RHS in cal_nlk is on the left side, i.e. -(vor x u) -chi*(u-us)
  ! its sign is inversed when computing the pressure. we therefore inverse the 
  ! sign here, which was wrong in old versions.
  
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
          pk(iz,ix,iy) = - imag*(kx*nlk(iz,ix,iy,1)+ky*nlk(iz,ix,iy,2)+&
                                 kz*nlk(iz,ix,iy,3) )/k2
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



! ------------------------------------------------------------------------------
! Vorticity sponge technology
!
! adds a damping term to the vorticity and thus removes partly the periodicity
!
! input is penalized vorticity in fourier space, output is the penalization 
! term to be added to the Navier--Stokes eqn.
!
! the vorticity sponge term is chi_sp * (vort) / eta_sponge
! but this is in the vorticity formulation
! so we compute the streamfunctions vort = - LAPLACE(psi)
! and then take the curl sponge = nabla \crossproduct psi
!
! currently, the sponge array is global so no arguments
! ------------------------------------------------------------------------------
subroutine vort_sponge()
  use mpi_header
  use fsi_vars
  complex (kind=pr) :: im, spx,spy,spz
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  
  ! imaginary unit
  im=dcmplx(0.d0,1.d0)  
  
  do iy=ca(3), cb(3)    ! ky : 0..ny/2-1 ,then, -ny/2..-1     
    ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
    ky2=ky*ky
    do ix=ca(2), cb(2)  ! kx : 0..nx/2
      kx=scalex*dble(ix)                
      kx2=kx*kx
      do iz=ca(1),cb(1)  ! kz : 0..nz/2-1 ,then, -nz/2..-1           
        kz     =scalez*dble(modulo(iz+nz/2,nz)-nz/2)
        kz2    =kz*kz
        k_abs_2=kx2+ky2+kz2
        if (abs(k_abs_2) .ne. 0.0) then  
          ! we first "solve" the poisson eqn 
          ! which gives us the streamfunction components
          spx = sponge(iz,ix,iy,1) / k_abs_2
          spy = sponge(iz,ix,iy,2) / k_abs_2
          spz = sponge(iz,ix,iy,3) / k_abs_2
          
          ! we then take the curl of the streamfunction
          sponge(iz,ix,iy,1)=im*(ky*spz - kz*spy)
          sponge(iz,ix,iy,2)=im*(kz*spx - kx*spz)
          sponge(iz,ix,iy,3)=im*(kx*spy - ky*spx)          
          
        else
          sponge(iz,ix,iy,1)=dcmplx(0.d0,0.d0)
          sponge(iz,ix,iy,2)=dcmplx(0.d0,0.d0)
          sponge(iz,ix,iy,3)=dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo
  
end subroutine vort_sponge


!-------------------------------------------------------------------------------
! vorticity penalization
!
! computes chi_sponge * (vort-vort0) / eta_sponge in physical space
! we currently do not allocate a mask_sponge array
!
! in this version, the mask is applied in a layer on top (in Z-direction)
! of the domain.
!-------------------------------------------------------------------------------
subroutine penalize_vort ( vort_penalized, vort )
  use mpi_header
  use fsi_vars
  
  ! input: vorticity in phys space
  real(kind=pr),intent(in):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  ! output: penalized vorticity in phys space
  real(kind=pr),intent(out):: vort_penalized(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: iz
  
  vort_penalized = 0.d0
  
  do iz=ra(3),rb(3) 
    ! note we currenly do not allocate a mask for this
    if (iz>nz-32) then
      vort_penalized(:,:,iz) = - vort(:,:,iz) / eps_sponge
    endif
  enddo    
end subroutine penalize_vort
  