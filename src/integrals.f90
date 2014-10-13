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
  use basic_operators
  use insect_module
  implicit none

  type(timetype), intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
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
  real(kind=pr) :: t1
  t1=MPI_wtime()
  
  !-----------------------------------------------------------------------------
  ! hydrodynamic forces (except for AB2_rigid_solid time stepper)
  !-----------------------------------------------------------------------------
  ! the stepper AB2_rigid_solid has to compute the drag at every time step, so
  ! we can skip the separate computation in INTEGRALS
  if (compute_forces==1 .and. iTimeMethodFluid/="AB2_rigid_solid" ) then
    t3 = MPI_wtime()    
    ! to compute the forces, we need the mask at time t. not we cannot suppose
    ! that mask after fluidtimestep is at time t, it is rather at t-dt, thus we
    ! have to reconstruct the mask now. solids are also at time t
    if(iMoving==1) call create_mask( time%time, mask, mask_color, us, Insect, beams )
    call cal_drag ( time, u, mask, mask_color, us, Insect )
    time_drag = time_drag + MPI_wtime() - t3
  endif
  
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
  
  !-----------------------------------------------------------------------------
  ! fluid energy and dissipation
  !-----------------------------------------------------------------------------
  ! total kinetic energy (including solid)
  call compute_energies(u(:,:,:,1:3),ekin,ekinx,ekiny,ekinz)
  ! fluid kinetic energy (excluding solid)
  nlk(:,:,:,1,1)=u(:,:,:,1)*(1.d0-mask*eps)
  nlk(:,:,:,2,1)=u(:,:,:,2)*(1.d0-mask*eps)
  nlk(:,:,:,3,1)=u(:,:,:,3)*(1.d0-mask*eps)
  call compute_energies(nlk(:,:,:,1:3,1),ekinf,ekinxf,ekinyf,ekinzf)

  ! compute dissipation rate
  call curl( u(:,:,:,1:3), nlk(:,:,:,1:3,1) )
  ! dissipation in the whole domain:
  call compute_energies(nlk(:,:,:,1:3,1),diss,dissx,dissy,dissz)
  ! again consider only fluid domain
  nlk(:,:,:,1,1)=nlk(:,:,:,1,1)*(1.d0-mask*eps)
  nlk(:,:,:,2,1)=nlk(:,:,:,2,1)*(1.d0-mask*eps)
  nlk(:,:,:,3,1)=nlk(:,:,:,3,1)*(1.d0-mask*eps)
  call compute_energies(nlk(:,:,:,1:3,1),dissf,dissxf,dissyf,disszf)
  
  if (iTimeMethodFluid=="AB2") then
    write(*,*) "ATTENTION write_integrals is NOT YET READY for AB2, it overwrites&
    & the nlk(:,:,:,1:3,1) vector which is only ok for RK schemes"
  endif
       
  ! add missing factor (from enstrophy to dissipation rate)
  diss  = 2.d0*nu*diss  ! note hidden factor of 2 in compute_energies
  dissx = 2.d0*nu*dissx
  dissy = 2.d0*nu*dissy
  dissz = 2.d0*nu*dissz
  dissf  = 2.d0*nu*dissf
  dissxf = 2.d0*nu*dissxf
  dissyf = 2.d0*nu*dissyf
  disszf = 2.d0*nu*disszf
   
  ! dump to disk     
  if(mpirank == 0) then
     open(14,file='energy.t',status='unknown',position='append')
     write (14,'(21(es15.8,1x))') time%time,&
       ekinf,ekinxf,ekinyf,ekinzf,&
       dissf,dissxf,dissyf,disszf,&
       ekin,ekinx,ekiny,ekinz,&
       diss,dissx,dissy,dissz,&
       GlobalIntegrals%penalization_power,& ! note this is computed in drag.f90
       GlobalIntegrals%penalization_power_x,&
       GlobalIntegrals%penalization_power_y,&
       GlobalIntegrals%penalization_power_z
     close(14)
  endif
  
  !-----------------------------------------------------------------------------
  ! Save mean flow values
  !-----------------------------------------------------------------------------
  uxmean = volume_integral(u(:,:,:,1)) / (xl*yl*zl)
  uymean = volume_integral(u(:,:,:,2)) / (xl*yl*zl)
  uzmean = volume_integral(u(:,:,:,3)) / (xl*yl*zl)
  if (mpirank==0) then
     open  (14,file='meanflow.t',status='unknown',position='append')
     write (14,'(4(es15.8,1x))') time%time, uxmean, uymean, uzmean
     close (14)
  endif
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
  !-----------------------------------------------------------------------------
  ! mask volume
  !-----------------------------------------------------------------------------
  mask = mask*eps
  call compute_mask_volume(mask,volume)
  mask = mask/eps
  if(mpirank == 0) then
    open(14,file='mask_volume.t',status='unknown',position='append')
    write (14,'(2(es15.8,1x))') time%time, volume
    close(14)
  endif
  
  
  time_integrals = time_integrals + MPI_wtime()-t1
end subroutine 



! Compute the average total energy and energy in each direction for a
! physical-space vector fields with components f1, f2, f3, leaving the
! input vector field untouched. ACTS ON FLUID DOMAIN ONLY
subroutine compute_energies_f(E,Ex,Ey,Ez,f1,f2,f3,mask)
  use mpi
  use vars
  implicit none
  
  real(kind=pr),intent(out) :: E,Ex,Ey,Ez
  real(kind=pr),intent(in):: f1(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(in):: f2(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(in):: f3(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(in):: mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr) :: LE,LEx,LEy,LEz ! local quantities
  real(kind=pr) :: v1,v2,v3
  integer :: ix,iy,iz,mpicode

  ! initialize local variables
  LE=0.d0
  LEx=0.d0
  LEy=0.d0
  LEz=0.d0

  ! Add contributions in physical space
  do iz=ra(3),rb(3)
     do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
           if(mask(ix,iy,iz) == 0.d0) then
              
              v1=f1(ix,iy,iz)
              v2=f2(ix,iy,iz)
              v3=f3(ix,iy,iz)
              
              LE=Le + v1*v1 + v2*v2 + v3*v3
              LEx=LEx + v1*v1
              LEy=LEy + v2*v2
              LEz=LEz + v3*v3
           endif
        enddo
     enddo
  enddo

  LE=0.5*dx*dy*dz*LE
  LEx=0.5*dx*dy*dz*LEx
  LEy=0.5*dx*dy*dz*LEy
  LEz=0.5*dx*dy*dz*LEz

  ! Sum over all MPI processes
  call MPI_ALLREDUCE(LE,E,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEx,Ex,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEy,Ey,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEz,Ez,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
end subroutine compute_energies_f


!-------------------------------------------------------------------------------
! Compute the average total energy and energy in each direction for a
! physical-space vector fields u(:,:,:,1:3), leaving the
! input vector field untouched.
!-------------------------------------------------------------------------------
subroutine compute_energies(u,E,Ex,Ey,Ez)
  use mpi
  use vars
  implicit none
  
  real(kind=pr),intent(in):: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:3)
  real(kind=pr),intent(out) :: E,Ex,Ey,Ez
  real(kind=pr) :: LE,LEx,LEy,LEz ! local quantities
  real(kind=pr) :: ux,uy,uz
  integer :: ix,iy,iz,mpicode

  ! initialize local variables
  LE=0.d0
  LEx=0.d0
  LEy=0.d0
  LEz=0.d0

  ! Add contributions in physical space
  do iz=ra(3),rb(3)
     do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          ux=u(ix,iy,iz,1)
          uy=u(ix,iy,iz,2)
          uz=u(ix,iy,iz,3)
          
          LEx=LEx + ux*ux
          LEy=LEy + uy*uy
          LEz=LEz + uz*uz
        enddo
     enddo
  enddo

  LE= 0.5*dx*dy*dz*(LEx+LEy+LEz)
  LEx=0.5*dx*dy*dz*LEx
  LEy=0.5*dx*dy*dz*LEy
  LEz=0.5*dx*dy*dz*LEz

  ! Sum over all MPI processes
  call MPI_ALLREDUCE(LE,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEx,Ex,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEy,Ey,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LEz,Ez,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       MPI_COMM_WORLD,mpicode)
end subroutine compute_energies

 
 
! Compute the fluid volume.
subroutine compute_fluid_volume(mask,volume)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(in)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
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

  real(kind=pr),intent(in)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(out) :: volume
  integer :: mpicode
  real(kind=pr) :: Lvolume ! Process-local volume

  Lvolume=sum(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))*dx*dy*dz
     
  call MPI_REDUCE(Lvolume,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,mpicode)
          
end subroutine compute_mask_volume