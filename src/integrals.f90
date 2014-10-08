!-------------------------------------------------------------------------------
! Wrapper for writing integral quantities to file
! Input:
!       uk: the neq-component vector of the unknowns in F-space
! Input/Output:
!       u, vort: two 3D-work arrays (free on entry and exit)
!       nlk: 3D complex work array (free on entry and exit)
! Output:
!       all output is done directly to hard disk in the *.t files
!-------------------------------------------------------------------------------
subroutine write_integrals(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  real(kind=pr) :: t1
  
  t1=MPI_wtime()
  call write_integrals_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)  
  time_integrals = time_integrals + MPI_wtime()-t1
end subroutine 


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use basic_operators
  use solid_model
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(in)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect 
  
  real(kind=pr) :: kx, ky, kz, maxdiv,maxdiv_fluid, maxdiv_loc,volume, t3
  real(kind=pr) :: concentration, conc
  real(kind=pr) :: ekinf, ekinxf, ekinyf, ekinzf
  real(kind=pr) :: ekin, ekinx, ekiny, ekinz
  real(kind=pr) :: diss, dissx, dissy, dissz
  real(kind=pr) :: dissf, dissxf, dissyf, disszf
  integer :: ix,iy,iz,mpicode
  
!   !-----------------------------------------------------------------------------
!   ! hydrodynamic forces (except for AB2_rigid_solid time stepper)
!   !-----------------------------------------------------------------------------
!   ! the stepper AB2_rigid_solid has to compute the drag at every time step, so
!   ! we can skip the separate computation in INTEGRALS
!   if (compute_forces==1 .and. iTimeMethodFluid/="AB2_rigid_solid" ) then
!     t3 = MPI_wtime()    
!     ! to compute the forces, we need the mask at time t. not we cannot suppose
!     ! that mask after fluidtimestep is at time t, it is rather at t-dt, thus we
!     ! have to reconstruct the mask now. solids are also at time t
!     if(iMoving==1) call create_mask(time, Insect, beams)
!     call cal_drag (time, u, Insect)
!     time_drag = time_drag + MPI_wtime() - t3
!   endif
  
  !-----------------------------------------------------------------------------
  ! divergence of velocity field (in the entire domain and in the fluid domain)
  !-----------------------------------------------------------------------------
  call divergence( u(:,:,:,1:3), work(:,:,:,1) )

  maxdiv = fieldmax( work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
  work(:,:,:,1) = work(:,:,:,1)*(1.d0-mask*eps)
  maxdiv_fluid = fieldmax( work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
  
  if(mpirank == 0) then
     open(14,file='divu.t',status='unknown',position='append')
     write (14,'(4(es15.8,1x))') time%time,maxdiv,maxdiv_fluid
     close(14)
  endif
  
!   !-----------------------------------------------------------------------------
!   ! fluid energy and dissipation
!   !-----------------------------------------------------------------------------
!   ! total kinetic energy (including solid)
!   call compute_energies(u(:,:,:,1:3),ekin,ekinx,ekiny,ekinz)
!   ! fluid kinetic energy (excluding solid)
!   u(:,:,:,1)=u(:,:,:,1)*(1.d0-mask*eps)
!   u(:,:,:,2)=u(:,:,:,2)*(1.d0-mask*eps)
!   u(:,:,:,3)=u(:,:,:,3)*(1.d0-mask*eps)
!   call compute_energies(u(:,:,:,1:3),ekinf,ekinxf,ekinyf,ekinzf)
! 
!   ! compute dissipation rate
!   call curl( uk, work3c )
!   call ifft3( ink=work3c, outx=work3r )
!   ! dissipation in the whole domain:
!   call compute_energies(work3r(:,:,:,1:3),diss,dissx,dissy,dissz)
!   ! again consider only fluid domain
!   work3r(:,:,:,1)=work3r(:,:,:,1)*(1.d0-mask*eps)
!   work3r(:,:,:,2)=work3r(:,:,:,2)*(1.d0-mask*eps)
!   work3r(:,:,:,3)=work3r(:,:,:,3)*(1.d0-mask*eps)
!   call compute_energies(work3r(:,:,:,1:3),dissf,dissxf,dissyf,disszf)
!        
!   ! add missing factor (from enstrophy to dissipation rate)
!   diss  = 2.d0*nu*diss  ! note hidden factor of 2 in compute_energies
!   dissx = 2.d0*nu*dissx
!   dissy = 2.d0*nu*dissy
!   dissz = 2.d0*nu*dissz
!   dissf  = 2.d0*nu*dissf
!   dissxf = 2.d0*nu*dissxf
!   dissyf = 2.d0*nu*dissyf
!   disszf = 2.d0*nu*disszf
!    
!   ! dump to disk     
!   if(mpirank == 0) then
!      open(14,file='energy.t',status='unknown',position='append')
!      write (14,'(21(es15.8,1x))') time,&
!        ekinf,ekinxf,ekinyf,ekinzf,&
!        dissf,dissxf,dissyf,disszf,&
!        ekin,ekinx,ekiny,ekinz,&
!        diss,dissx,dissy,dissz,&
!        GlobalIntegrals%penalization_power,& ! note this is computed in drag.f90
!        GlobalIntegrals%penalization_power_x,&
!        GlobalIntegrals%penalization_power_y,&
!        GlobalIntegrals%penalization_power_z
!      close(14)
!   endif
!   
!   !-----------------------------------------------------------------------------
!   ! Save mean flow values
!   !-----------------------------------------------------------------------------
!   if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
!      ! This is done only by one CPU (which is not nessesarily the root rank)
!      open  (14,file='meanflow.t',status='unknown',position='append')
!      write (14,'(4(es15.8,1x))') time, dreal(uk(0,0,0,1:3))
!      close (14)
!   endif
!   
!   !-----------------------------------------------------------------------------
!   ! integral of scalar concentration
!   !-----------------------------------------------------------------------------
!   if ((use_passive_scalar==1).and.(compute_scalar)) then
!     call ifft( ink=uk(:,:,:,4), outx=work1)
!     work1 = work1*(1.d0-mask*eps)
!     conc = sum(work1)*dx*dy*dz
!     call MPI_REDUCE(conc,concentration,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!          MPI_COMM_WORLD,mpicode)
!          
!     if (mpirank == 0) then
!       open  (14,file='scalar.t',status='unknown',position='append')
!       write (14,'(2(es15.8,1x))') time, concentration
!       close (14)
!     endif
!   endif
!   
!   !-----------------------------------------------------------------------------
!   ! mask volume
!   !-----------------------------------------------------------------------------
!   mask = mask*eps
!   call compute_mask_volume(volume)
!   mask = mask/eps
!   if(mpirank == 0) then
!     open(14,file='mask_volume.t',status='unknown',position='append')
!     write (14,'(2(es15.8,1x))') time,volume
!     close(14)
!   endif
  
end subroutine 
! 
! 

! ! 
! ! 
! ! Compute the average total energy and energy in each direction for a
! ! physical-space vector fields with components f1, f2, f3, leaving the
! ! input vector field untouched. ACTS ON FLUID DOMAIN ONLY
! subroutine compute_energies_f(E,Ex,Ey,Ez,f1,f2,f3)
!   use mpi
!   use vars
!   implicit none
!   
!   real(kind=pr),intent(out) :: E,Ex,Ey,Ez
!   real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr) :: LE,LEx,LEy,LEz ! local quantities
!   real(kind=pr) :: v1,v2,v3
!   integer :: ix,iy,iz,mpicode
! 
!   ! initialize local variables
!   LE=0.d0
!   LEx=0.d0
!   LEy=0.d0
!   LEz=0.d0
! 
!   ! Add contributions in physical space
!   do ix=ra(1),rb(1)
!      do iy=ra(2),rb(2)
!         do iz=ra(3),rb(3)
!            if(mask(ix,iy,iz) == 0.d0) then
!               
!               v1=f1(ix,iy,iz)
!               v2=f2(ix,iy,iz)
!               v3=f3(ix,iy,iz)
!               
!               LE=Le + v1*v1 + v2*v2 + v3*v3
!               LEx=LEx + v1*v1
!               LEy=LEy + v2*v2
!               LEz=LEz + v3*v3
!            endif
!         enddo
!      enddo
!   enddo
! 
!   LE=0.5*dx*dy*dz*LE
!   LEx=0.5*dx*dy*dz*LEx
!   LEy=0.5*dx*dy*dz*LEy
!   LEz=0.5*dx*dy*dz*LEz
! 
!   ! Sum over all MPI processes
!   call MPI_ALLREDUCE(LE,E,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEx,Ex,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEy,Ey,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEz,Ez,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
! end subroutine compute_energies_f
! 
! 
! !-------------------------------------------------------------------------------
! ! Compute the average total energy and energy in each direction for a
! ! physical-space vector fields u(:,:,:,1:3), leaving the
! ! input vector field untouched.
! !-------------------------------------------------------------------------------
! subroutine compute_energies(u,E,Ex,Ey,Ez)
!   use mpi
!   use vars
!   implicit none
!   
!   real(kind=pr),intent(in):: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
!   real(kind=pr),intent(out) :: E,Ex,Ey,Ez
!   real(kind=pr) :: LE,LEx,LEy,LEz ! local quantities
!   real(kind=pr) :: ux,uy,uz
!   integer :: ix,iy,iz,mpicode
! 
!   ! initialize local variables
!   LE=0.d0
!   LEx=0.d0
!   LEy=0.d0
!   LEz=0.d0
! 
!   ! Add contributions in physical space
!   do ix=ra(1),rb(1)
!      do iy=ra(2),rb(2)
!         do iz=ra(3),rb(3)
!           ux=u(ix,iy,iz,1)
!           uy=u(ix,iy,iz,2)
!           uz=u(ix,iy,iz,3)
!           
!           LEx=LEx + ux*ux
!           LEy=LEy + uy*uy
!           LEz=LEz + uz*uz
!         enddo
!      enddo
!   enddo
! 
!   LE= 0.5*dx*dy*dz*(LEx+LEy+LEz)
!   LEx=0.5*dx*dy*dz*LEx
!   LEy=0.5*dx*dy*dz*LEy
!   LEz=0.5*dx*dy*dz*LEz
! 
!   ! Sum over all MPI processes
!   call MPI_ALLREDUCE(LE,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEx,Ex,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEy,Ey,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_ALLREDUCE(LEz,Ez,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
!        MPI_COMM_WORLD,mpicode)
! end subroutine compute_energies
! 
! 
! ! Compute the average average component in each direction for a
! ! physical-space vector fields with components f1, f2, f3, leaving the
! ! input vector field untouched.
! subroutine compute_components(mask,Cx,Cy,Cz,f1,f2,f3)
!   use mpi
!   use vars
!   implicit none
!   
!   real(kind=pr),intent(out) :: Cx,Cy,Cz
!   real(kind=pr),intent(inout):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(inout):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(inout):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr) :: LCx,LCy,LCz ! local quantities
!   real(kind=pr) :: v1,v2,v3
!   integer :: ix,iy,iz,mpicode
! 
!   ! initialize local variables
!   LCx=0.d0
!   LCy=0.d0
!   LCz=0.d0
! 
!   ! Add contributions in physical space
!   do ix=ra(1),rb(1)
!      do iy=ra(2),rb(2)
!         do iz=ra(3),rb(3)
!            if(mask(ix,iy,iz) == 0.d0) then
!               v1=f1(ix,iy,iz)
!               v2=f2(ix,iy,iz)
!               v3=f3(ix,iy,iz)
!               
!               LCx=LCx + v1
!               LCy=LCy + v2
!               LCz=LCz + v3
!            endif
!         enddo
!      enddo
!   enddo
! 
!   LCx=LCx*dx*dy*dz
!   LCy=LCy*dx*dy*dz
!   LCz=LCz*dx*dy*dz
!   
!   ! Sum over all MPI processes
!   call MPI_REDUCE(LCx,Cx,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_REDUCE(LCy,Cy,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_REDUCE(LCz,Cz,&
!        1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
!        MPI_COMM_WORLD,mpicode)
! end subroutine compute_components

! 
! 
! ! Compute the maximum non-normalized divergence of the given 3D field
! ! fk1, fk2, fk3, 
! subroutine compute_max_div(maxdiv,fk1,fk2,fk3,f1,f2,f3,div,divk)
!   use mpi
!   use vars
!   use p3dfft_wrapper
!   implicit none
! 
!   real(kind=pr),intent(out) :: maxdiv  
!   real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   complex(kind=pr),intent(in) ::fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
!   complex(kind=pr),intent(in) ::fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
!   complex(kind=pr),intent(in) ::fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
!   complex(kind=pr),intent(inout) ::divk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
!   real(kind=pr),intent(inout):: div(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   integer :: ix,iy,iz,mpicode
!   real(kind=pr) :: kx, ky, kz, locmax,  v1,v2,v3,d
!   ! real(kind=pr) fnorm ! Only used for normalized version.
!   complex(kind=pr) :: imag ! imaginary unit
! 
!   imag = dcmplx(0.d0,1.d0)
! 
!   ! Compute the divergence in Fourier space, store in divk
!   do iz=ca(1),cb(1)          
!     !-- wavenumber in z-direction
!     kz = wave_z(iz)       
!     do iy=ca(2), cb(2)
!       !-- wavenumber in y-direction
!       ky = wave_y(iy)      
!       do ix=ca(3), cb(3)
!         !-- wavenumber in x-direction
!         kx = wave_x(ix)
! 
!         divk(iz,iy,ix)=imag*&
!             (kx*fk1(iz,iy,ix)&
!             +ky*fk2(iz,iy,ix)&
!             +kz*fk3(iz,iy,ix))
!       enddo
!     enddo
!   enddo
! 
!   call ifft(div,divk)
!   
!   ! Find the local max
!   
!   ! FIXME: at least in the present version, this can be simplified to
!   ! locmax = maxval(abs(div))
!   ! without loss of functionality or performance.
!   
!   locmax=0.d0
!   do ix=ra(1),rb(1)
!      do iy=ra(2),rb(2)
!         do iz=ra(3),rb(3)
!            if(mask(ix,iy,iz) == 0.d0) then
!               
!               v1=f1(ix,iy,iz)
!               v2=f2(ix,iy,iz)
!               v3=f3(ix,iy,iz)
!               
!               ! Normalized version:
!               ! fnorm=v1*v2 + v2*v2 + v3*v3 + 1d-8 ! avoid division by zero
!               ! d=abs(div(ix,iy,iz))/fnorm
!               
!               ! Non-normalized version:
!               d=abs(div(ix,iy,iz))
!               
!               if(d > locmax) then
!                  locmax=d
!               endif
!            endif
!         enddo
!      enddo
!   enddo
! 
!   ! Find the global max
!   call MPI_REDUCE(locmax,maxdiv,&
!        1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
!        MPI_COMM_WORLD,mpicode)
! end subroutine compute_max_div
! 
! 
! ! Compute the maximum components of the given 3D field with
! ! componennts f1, f2, f3.
! subroutine compute_max(vmax,xmax,ymax,zmax,f1,f2,f3)
!   use mpi
!   use vars
!   implicit none
! 
!   real(kind=pr),intent(out) :: vmax,xmax,ymax,zmax
!   real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
!   integer :: ix,iy,iz,mpicode
!   real(kind=pr) :: v1,v2,v3
!   real(kind=pr) :: Lmax,Lxmax,Lymax,Lzmax
! 
!   Lmax=0.d0
!   Lxmax=0.d0
!   Lymax=0.d0
!   Lzmax=0.d0
! 
!   ! Find the (per-process) max norm and max components in physical
!   ! space
!   do ix=ra(1),rb(1)
!      do iy=ra(2),rb(2)
!         do iz=ra(3),rb(3)
!            if(mask(ix,iy,iz) == 0.d0) then
!               v1=f1(ix,iy,iz)
!               v2=f2(ix,iy,iz)
!               v3=f3(ix,iy,iz)
!               Lmax=max(Lmax,dsqrt(v1*v1 + v2*v2 + v3*v3))
!               Lxmax=max(Lxmax,v1)
!               Lymax=max(Lymax,v2)
!               Lzmax=max(Lzmax,v3)
!            endif
!         enddo
!      enddo
!   enddo
! 
!   ! Determine the global max
!   call MPI_REDUCE(Lmax,vmax,&
!        1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_REDUCE(Lxmax,xmax,&
!        1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_REDUCE(Lymax,ymax,&
!        1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
!        MPI_COMM_WORLD,mpicode)
!   call MPI_REDUCE(Lzmax,zmax,&
!        1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
!        MPI_COMM_WORLD,mpicode)
! end subroutine compute_max
! 
! 
 
 
! Compute the fluid volume.
subroutine compute_fluid_volume(mask,volume)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(in)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out) :: volume
  integer :: ix,iy,iz,mpicode
  real(kind=pr) :: Lvolume ! Process-local volume
  real(kind=pr) :: dxyz ! Volume of pixel

  Lvolume=0.d0
  dxyz=dx*dy*dz

  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
          if(mask(ix,iy,iz) == 0.d0) then
              Lvolume=Lvolume +dxyz
          endif
        enddo
    enddo
  enddo
  
  call MPI_REDUCE(Lvolume,volume,&
      1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
      MPI_COMM_WORLD,mpicode)
end subroutine compute_fluid_volume


! Compute the mask volume.
subroutine compute_mask_volume(mask,volume)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(in)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out) :: volume
  integer :: mpicode
  real(kind=pr) :: Lvolume ! Process-local volume

  Lvolume=sum(mask)*dx*dy*dz
     
  call MPI_REDUCE(Lvolume,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,mpicode)
          
end subroutine compute_mask_volume