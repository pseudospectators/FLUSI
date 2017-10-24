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
subroutine write_integrals(time,uk,u,vort,nlk,work,scalars,Insect,beams)
  use mpi
  use vars
  use solid_model
  use insect_module
  implicit none

  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(in):: time
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: t1

  t1=MPI_wtime()

  select case(method)
  case("fsi")
    call write_integrals_fsi(time,uk,u,vort,nlk,work(:,:,:,1),scalars,Insect,beams)
  case("mhd")
    call write_integrals_mhd(time,uk,u,vort,nlk,work(:,:,:,1))
  case default
    call abort(1, "Error! Unknown method in write_integrals")
  end select

  time_integrals = time_integrals + MPI_wtime()-t1
end subroutine write_integrals


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time,uk,u,work3r,work3c,work1,scalars,Insect,beams)
  use mpi
  use vars
  use p3dfft_wrapper
  use basic_operators
  use solid_model
  use insect_module
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in)::time
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::work3r(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::work1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::work3c(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect

  real(kind=pr) :: kx, ky, kz, maxdiv,maxdiv_fluid, maxdiv_loc,volume, t3
  real(kind=pr) :: concentration, conc, maxi
  real(kind=pr) :: ekinf, ekinxf, ekinyf, ekinzf
  real(kind=pr) :: ekin, ekinx, ekiny, ekinz
  real(kind=pr) :: diss, dissx, dissy, dissz
  real(kind=pr) :: dissf, dissxf, dissyf, disszf
  real(kind=pr), dimension(0:nx-1) :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec
  integer :: ix,iy,iz,mpicode,j
  character(len=9) scalarfile

  ! fetch u in x-space at time t. note the output u of fluidtimestep is not
  ! nessesarily what we need, so we must ensure u=ifft(uk) here
  call ifft3 (ink=uk, outx=u)

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
    if(iMoving==1) call create_mask(time, Insect, beams)
    call cal_drag (time, u, Insect)
    time_drag = time_drag + MPI_wtime() - t3
  endif

  !-----------------------------------------------------------------------------
  ! divergence of velocity field (in the entire domain and in the fluid domain)
  !-----------------------------------------------------------------------------
  call divergence( ink=uk, outk=work3c(:,:,:,1) )
  call ifft( ink=work3c(:,:,:,1), outx=work1 ) ! work1 is now div in phys space

  maxdiv = fieldmax(work1)
  maxdiv_fluid = fieldmax(work1*(1.d0-mask*eps))
  if(mpirank == 0) then
    open(14,file='divu.t',status='unknown',position='append')
    write (14,'(4(es15.8,1x))') time,maxdiv,maxdiv_fluid
    close(14)
  endif

  !-----------------------------------------------------------------------------
  ! fluid energy and dissipation
  !-----------------------------------------------------------------------------
  ! total kinetic energy (including solid)
  call compute_energies(u(:,:,:,1:3),ekin,ekinx,ekiny,ekinz)
  ! fluid kinetic energy (excluding solid)
  u(:,:,:,1)=u(:,:,:,1)*(1.d0-mask*eps)
  u(:,:,:,2)=u(:,:,:,2)*(1.d0-mask*eps)
  u(:,:,:,3)=u(:,:,:,3)*(1.d0-mask*eps)
  call compute_energies(u(:,:,:,1:3),ekinf,ekinxf,ekinyf,ekinzf)

  ! compute dissipation rate
  call curl( uk, work3c )
  call ifft3( ink=work3c, outx=work3r )
  ! dissipation in the whole domain:
  call compute_energies(work3r(:,:,:,1:3),diss,dissx,dissy,dissz)
  ! again consider only fluid domain
  work3r(:,:,:,1)=work3r(:,:,:,1)*(1.d0-mask*eps)
  work3r(:,:,:,2)=work3r(:,:,:,2)*(1.d0-mask*eps)
  work3r(:,:,:,3)=work3r(:,:,:,3)*(1.d0-mask*eps)
  call compute_energies(work3r(:,:,:,1:3),dissf,dissxf,dissyf,disszf)

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
    write (14,'(21(es15.8,1x))') time,&
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
  if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
    ! This is done only by one CPU (which is not nessesarily the root rank)
    open  (14,file='meanflow.t',status='unknown',position='append')
    write (14,'(4(es15.8,1x))') time, dreal(uk(0,0,0,1:3))
    close (14)
  endif

  !-----------------------------------------------------------------------------
  ! integral of scalar concentration
  !-----------------------------------------------------------------------------
  if ((use_passive_scalar==1).and.(compute_scalar)) then
    do j=1,n_scalars
      ! exclude solid regions
      work1 = scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j)*(1.d0-mask*eps)
      conc = sum(work1)*dx*dy*dz
      call MPI_REDUCE(conc,concentration,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
      MPI_COMM_WORLD,mpicode)

      maxi = fieldmax(scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j))
      write(scalarfile,'("scalar",i1,".t")') j
      if (mpirank == 0) then
        open  (14,file=scalarfile,status='unknown',position='append')
        write (14,'(3(es15.8,1x))') time, concentration, maxi
        close (14)
      endif
    enddo
  endif

  !-----------------------------------------------------------------------------
  ! mask volume
  !-----------------------------------------------------------------------------
  mask = mask*eps
  ! the mask volume is useful for debugging, since we woul see some changes if occasionally
  ! the code does not set some values at some grid points
  call compute_mask_volume(volume)
  mask = mask/eps
  ! the mask is not everything: the solid velocity has to be correct too, so here we
  ! compute another integral quantity:
  kx = mpisum( sum(mask*us(:,:,:,1))*dx*dy*dz*eps )
  ky = mpisum( sum(mask*us(:,:,:,2))*dx*dy*dz*eps )
  kz = mpisum( sum(mask*us(:,:,:,3))*dx*dy*dz*eps )

  if(mpirank == 0) then
    open(14,file='mask_volume.t',status='unknown',position='append')
    write (14,'(5(es15.8,1x))') time, volume, kx, ky, kz
    close(14)
  endif

  !-----------------------------------------------------------------------------
  ! Spectrum
  !-----------------------------------------------------------------------------
  if (iSaveSpectrae=="yes") then
    call compute_spectrum(time,kvec,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
    if(mpirank == 0) then
      open(14,file='spectrum_tot.t',status='unknown',position='append')
      write (14,'(4096(es15.8,1x))') time,S_Ekin
      close(14)
      open(14,file='spectrum_x.t',status='unknown',position='append')
      write (14,'(4096(es15.8,1x))') time,S_Ekinx
      close(14)
      open(14,file='spectrum_y.t',status='unknown',position='append')
      write (14,'(4096(es15.8,1x))') time,S_Ekiny
      close(14)
      open(14,file='spectrum_z.t',status='unknown',position='append')
      write (14,'(4096(es15.8,1x))') time,S_Ekinz
      close(14)
      ! the wavenumber vector does not change, so just keep one line
      open(14,file='spectrum_k.t', status='replace')
      write (14,'(4096(es15.8,1x))') kvec
      close(14)
    endif
  endif
end subroutine write_integrals_fsi


! The mhd version of writing integral quantities to disk.
! In order to make the asy files useful for both hd and mhd codes,
! please output velocity and magnetic fields quantities in separate
! files, or (if there aren't too many columns) put all the
! velocity-only quantities first.
subroutine write_integrals_mhd(time,ubk,ub,wj,nlk,work)
  use mpi
  use vars
  use penalization ! mask array etc
  use p3dfft_wrapper
  use basic_operators
  implicit none

  complex (kind=pr),intent(in)::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr), intent(in) :: time
  integer :: i
  ! Local loop variables
  real(kind=pr) :: Ekin,Ekinx,Ekiny,Ekinz
  real(kind=pr) :: Emag,Emagx,Emagy,Emagz
  real(kind=pr) :: meanjx,meanjy,meanjz
  real(kind=pr) :: jmax,jxmax,jymax,jzmax
  real(kind=pr) :: divu,divb
  real(kind=pr) :: dissu,dissb
  real(kind=pr) :: fluid_volume

  ! Make sure that we have the fields that we need in the space we need:

  ! Compute u and B to physical space
  do i=1,nd
    call ifft(ub(:,:,:,i),ubk(:,:,:,i))
  enddo

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

  ! Compute the integral quantities and output to disk:

  ! Compute the fluid volume.
  call compute_fluid_volume(fluid_volume)

  ! Compute kinetic energies.
  call compute_energies_f(Ekin,Ekinx,Ekiny,Ekinz,&
  ub(:,:,:,1),ub(:,:,:,2),ub(:,:,:,3))
  if(mpirank == 0) then
    open(14,file='ek.t',status='unknown',position='append')
    ! 9 outputs, including tabs
    write(14,'(e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16)') &
    time,tab,Ekin,tab,Ekinx,tab,Ekiny,tab,Ekinz
    close(14)
  endif

  ! Comptue magnetic energies.
  call compute_energies_f(Emag,Emagx,Emagy,Emagz,&
  ub(:,:,:,4),ub(:,:,:,5),ub(:,:,:,6))
  if(mpirank == 0) then
    open(14,file='eb.t',status='unknown',position='append')
    ! 9 outputs, including tabs
    write(14,'(e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16)') &
    time,tab,Emag,tab,Emagx,tab,Emagy,tab,Emagz
    close(14)
  endif

  ! Compute current density values.
  call compute_components(meanjx,meanjy,meanjz,&
  wj(:,:,:,4),wj(:,:,:,5),wj(:,:,:,6))
  call compute_max(jmax,jxmax,jymax,jzmax,wj(:,:,:,4),wj(:,:,:,5),wj(:,:,:,6))
  if(mpirank == 0) then
    open(14,file='j.t',status='unknown',position='append')
    ! 15 outputs, including tabs
    write(14,&
    '(e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16,A,e25.16)')&
    time,tab,meanjx,tab,meanjy,tab,meanjz,tab,jmax,tab,jxmax,tab,&
    jymax,tab,jzmax
    close(14)
  endif

  ! Compute kinetic and magnetic energy dissipation
  ! Kinetic energy dissipation is nu*< |vorticity| >
  call compute_mean_norm(dissu,wj(:,:,:,1),wj(:,:,:,2),wj(:,:,:,3))
  dissu=nu*dissu
  ! Magnetic energy dissipation is eta*< |current density| >
  call compute_mean_norm(dissb,wj(:,:,:,4),wj(:,:,:,5),wj(:,:,:,6))
  dissb=eta*dissb
  if(mpirank == 0) then
    open(14,file='diss.t',status='unknown',position='append')
    ! 3 outputs
    write(14,'(e25.16,A,e25.16,A,e25.16)') time,tab,dissu,tab,dissb
    close(14)
  endif

  ! Compute max divergence.
  call compute_max_div(divu,&
  ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3),&
  ub(:,:,:,1),ub(:,:,:,2),ub(:,:,:,3),&
  work,nlk(:,:,:,1))
  call compute_max_div(divb,&
  ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6),&
  ub(:,:,:,4),ub(:,:,:,5),ub(:,:,:,6),&
  work,nlk(:,:,:,1))
  if(mpirank == 0) then
    open(14,file='d.t',status='unknown',position='append')
    ! 3 outputs
    write(14,'(e25.16,A,e25.16,A,e25.16)') time,tab,divu,tab,divb
    close(14)
  endif
end subroutine write_integrals_mhd


! Compute the average total energy and energy in each direction for a
! physical-space vector fields with components f1, f2, f3, leaving the
! input vector field untouched. ACTS ON FLUID DOMAIN ONLY
subroutine compute_energies_f(E,Ex,Ey,Ez,f1,f2,f3)
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: E,Ex,Ey,Ez
  real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
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

  LE=0.50d0*dx*dy*dz*LE
  LEx=0.50d0*dx*dy*dz*LEx
  LEy=0.50d0*dx*dy*dz*LEy
  LEz=0.50d0*dx*dy*dz*LEz

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

  real(kind=pr),intent(in):: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
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

  LE= 0.50d0*dx*dy*dz*(LEx+LEy+LEz)
  LEx=0.50d0*dx*dy*dz*LEx
  LEy=0.50d0*dx*dy*dz*LEy
  LEz=0.50d0*dx*dy*dz*LEz

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


!-------------------------------------------------------------------------------
! Compute the average total energy for a physical-space SCALAR fiels u, leaving the
! input vector field untouched.
!-------------------------------------------------------------------------------
subroutine compute_energies1(u,E)
  use mpi
  use vars
  implicit none

  real(kind=pr),intent(in):: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(out) :: E
  real(kind=pr) :: LE ! local quantities
  real(kind=pr) :: ux
  integer :: ix,iy,iz,mpicode

  ! initialize local variables
  LE=0.d0

  ! Add contributions in physical space
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ux=u(ix,iy,iz)
        LE = LE + ux*ux
      enddo
    enddo
  enddo

  LE = 0.50d0*dx*dy*dz*LE

  ! Sum over all MPI processes
  call MPI_ALLREDUCE(LE,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
end subroutine compute_energies1


!-------------------------------------------------------------------------------
! Compute the integral total energy for a input field given in Fourier space
!-------------------------------------------------------------------------------
subroutine compute_energies_k(uk,E)
    use vars
    use helpers
    implicit none
    complex(kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
    integer :: ix,iy,iz
    real(kind=pr) :: E,e1

    ! compute total kinetic energy (including the solid domain) in Fourier space
    ! using Parseval's identity
    e1 = 0.d0
    do iz=ca(1),cb(1)
      do iy=ca(2),cb(2)
        do ix=ca(3),cb(3)
          if ( ix==0 .or. ix==nx/2 ) then
            ! note zero mode and highest wavenumber is special
            e1=e1+dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2)/2. &
                 +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2)/2. &
                 +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)/2.
          else
            e1=e1+dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2) &
                 +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2) &
                 +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)
          endif
        enddo
      enddo
    enddo

    ! note the scaling factor for parsevals identity
    E = mpisum(e1)*xl*yl*zl
end subroutine compute_energies_k

!-------------------------------------------------------------------------------
! Compute the integral total energy for a input field given in Fourier space
!-------------------------------------------------------------------------------
subroutine compute_energies1_k(uk,E)
    use vars
    use helpers
    implicit none
    complex(kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
    integer :: ix,iy,iz
    real(kind=pr) :: E,e1

    ! compute total kinetic energy (including the solid domain) in Fourier space
    ! using Parseval's identity
    e1 = 0.d0
    do iz=ca(1),cb(1)
      do iy=ca(2),cb(2)
        do ix=ca(3),cb(3)
          if ( ix==0 .or. ix==nx/2 ) then
            ! note zero mode and highest wavenumber is special
            e1=e1+dble(real(uk(iz,iy,ix))**2+aimag(uk(iz,iy,ix))**2)/2.
          else
            e1=e1+dble(real(uk(iz,iy,ix))**2+aimag(uk(iz,iy,ix))**2)
          endif
        enddo
      enddo
    enddo

    ! note the scaling factor for parsevals identity
    E = mpisum(e1)*xl*yl*zl
end subroutine compute_energies1_k

! Compute the average average component in each direction for a
! physical-space vector fields with components f1, f2, f3, leaving the
! input vector field untouched.
subroutine compute_components(Cx,Cy,Cz,f1,f2,f3)
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: Cx,Cy,Cz
  real(kind=pr),intent(inout):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: LCx,LCy,LCz ! local quantities
  real(kind=pr) :: v1,v2,v3
  integer :: ix,iy,iz,mpicode

  ! initialize local variables
  LCx=0.d0
  LCy=0.d0
  LCz=0.d0

  ! Add contributions in physical space
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        if(mask(ix,iy,iz) == 0.d0) then
          v1=f1(ix,iy,iz)
          v2=f2(ix,iy,iz)
          v3=f3(ix,iy,iz)

          LCx=LCx + v1
          LCy=LCy + v2
          LCz=LCz + v3
        endif
      enddo
    enddo
  enddo

  LCx=LCx*dx*dy*dz
  LCy=LCy*dx*dy*dz
  LCz=LCz*dx*dy*dz

  ! Sum over all MPI processes
  call MPI_ALLREDUCE(LCx,Cx,&
  1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LCy,Cy,&
  1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(LCz,Cz,&
  1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
end subroutine compute_components


! Compute the maximum non-normalized divergence of the given 3D field
! fk1, fk2, fk3,
subroutine compute_max_div(maxdiv,fk1,fk2,fk3,f1,f2,f3,div,divk)
  use mpi
  use vars
  use p3dfft_wrapper
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: maxdiv
  real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(in) ::fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in) ::fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(in) ::fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout) ::divk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr),intent(inout):: div(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz,mpicode
  real(kind=pr) :: kx, ky, kz, locmax,  v1,v2,v3,d
  ! real(kind=pr) fnorm ! Only used for normalized version.
  complex(kind=pr) :: imag ! imaginary unit

  imag = dcmplx(0.d0,1.d0)

  ! Compute the divergence in Fourier space, store in divk
  do iz=ca(1),cb(1)
    !-- wavenumber in z-direction
    kz = wave_z(iz)
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)

        divk(iz,iy,ix)=imag*&
        (kx*fk1(iz,iy,ix)&
        +ky*fk2(iz,iy,ix)&
        +kz*fk3(iz,iy,ix))
      enddo
    enddo
  enddo

  call ifft(div,divk)

  ! Find the local max

  ! FIXME: at least in the present version, this can be simplified to
  ! locmax = maxval(abs(div))
  ! without loss of functionality or performance.

  locmax=0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        if(mask(ix,iy,iz) == 0.d0) then

          v1=f1(ix,iy,iz)
          v2=f2(ix,iy,iz)
          v3=f3(ix,iy,iz)

          ! Normalized version:
          ! fnorm=v1*v2 + v2*v2 + v3*v3 + 1d-8 ! avoid division by zero
          ! d=abs(div(ix,iy,iz))/fnorm

          ! Non-normalized version:
          d=abs(div(ix,iy,iz))

          if(d > locmax) then
            locmax=d
          endif
        endif
      enddo
    enddo
  enddo

  ! Find the global max
  call MPI_ALLREDUCE(locmax,maxdiv,&
  1,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)
end subroutine compute_max_div


! Compute the maximum components of the given 3D field with
! componennts f1, f2, f3.
subroutine compute_max(vmax,xmax,ymax,zmax,f1,f2,f3)
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: vmax,xmax,ymax,zmax
  real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz,mpicode
  real(kind=pr) :: v1,v2,v3
  real(kind=pr) :: Lmax,Lxmax,Lymax,Lzmax

  Lmax=0.d0
  Lxmax=0.d0
  Lymax=0.d0
  Lzmax=0.d0

  ! Find the (per-process) max norm and max components in physical
  ! space
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        if(mask(ix,iy,iz) == 0.d0) then
          v1=f1(ix,iy,iz)
          v2=f2(ix,iy,iz)
          v3=f3(ix,iy,iz)
          Lmax=max(Lmax,dsqrt(v1*v1 + v2*v2 + v3*v3))
          Lxmax=max(Lxmax,v1)
          Lymax=max(Lymax,v2)
          Lzmax=max(Lzmax,v3)
        endif
      enddo
    enddo
  enddo

  ! Determine the global max
  call MPI_ALLREDUCE(Lmax,vmax,&
  1,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(Lxmax,xmax,&
  1,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(Lymax,ymax,&
  1,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(Lzmax,zmax,&
  1,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)
end subroutine compute_max


! Compute the meannorm of the given field with x-space components f1, f2, f3.
subroutine compute_mean_norm(mean,f1,f2,f3)
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: mean
  real(kind=pr),intent(in):: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer :: ix,iy,iz,mpicode
  real(kind=pr) :: v1,v2,v3
  real(kind=pr) :: Lmean ! Process-local mean

  Lmean = 0.d0
  mean = 0.d0

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        if(mask(ix,iy,iz) == 0.d0) then
          v1=f1(ix,iy,iz)
          v2=f2(ix,iy,iz)
          v3=f3(ix,iy,iz)

          Lmean=Lmean + v1*v1 + v2*v2 + v3*v3
        endif
      enddo
    enddo
  enddo

  Lmean=Lmean*dx*dy*dz

  call MPI_ALLREDUCE(Lmean,mean,&
  1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
end subroutine compute_mean_norm


! Compute the fluid volume.
! NB: mask is a global!
subroutine compute_fluid_volume(volume)
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(out) :: volume
  integer :: ix,iy,iz,mpicode
  real(kind=pr) :: Lvolume ! Process-local volume
  real(kind=pr) :: dxyz ! Volume of pixel

  if(iPenalization == 0) then
    volume=xl*yl*zl
  else
    Lvolume=0.d0
    dxyz=dx*dy*dz

    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          if(mask(ix,iy,iz) == 0.d0) then
            Lvolume=Lvolume +dxyz
          endif
        enddo
      enddo
    enddo

    call MPI_ALLREDUCE(Lvolume,volume,&
    1,MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,mpicode)
  endif
end subroutine compute_fluid_volume



! Compute the fluid volume.
! NB: mask is a global!
subroutine compute_mask_volume(volume)
  use mpi
  use penalization
  use vars
  implicit none

  real(kind=pr),intent(out) :: volume
  integer :: mpicode
  real(kind=pr) :: Lvolume ! Process-local volume

  Lvolume=sum(mask)*dx*dy*dz

  call MPI_ALLREDUCE(Lvolume,volume,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

end subroutine compute_mask_volume


! ----------------------------------------------------------------------------
! Compute the energy spectrum of a velocity field in in Fourier space
! It is assumed that the x-direction is not split among processes, i.e. we
! deal with a 2D data decomposition.
!
! This code also works for non-cubic data
!
! The returned array "kvec" contains the "scaled wavenumber" (on the actual,
! physical domain and not the 2*pi interval)
!
! The wavenumber ordering is given by the wrapper is p3dfft_wrapper
! Assumes using STRIDE-1 ordering: (ix,iy,iz) but (kz,ky,kx) and NOT (kz,kx,ky)
! ----------------------------------------------------------------------------
subroutine compute_spectrum(time,kvec,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  real(kind=pr), intent(in) :: time
  complex(kind=pr),intent(inout) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  ! NOTE: we assume the x-direction to be contiguous (i.e. it is NOT split among CPU)
  real(kind=pr), dimension(0:nx-1), intent(inout) :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin,kvec

  real(kind=pr), dimension(0:nx-1) :: S_Ekinx_loc,S_Ekiny_loc,S_Ekinz_loc,S_Ekin_loc
  real(kind=pr) :: kx, ky, kz, kreal, kmax, dk
  real(kind=pr) :: sum_ux, sum_uy, sum_uz, sum_u
  integer :: ix, iy, iz, ik, mpicode, nk


  ! initalize local and global arrays. note all CPU hold an (0:nx-1) array (since
  ! in flusi, the x-direction is always local)
  S_Ekinx=0.0d0
  S_Ekiny=0.0d0
  S_Ekinz=0.0d0
  S_Ekin =0.0d0

  kvec = 0.d0

  S_Ekinx_loc=0.0d0
  S_Ekiny_loc=0.0d0
  S_Ekinz_loc=0.0d0
  S_Ekin_loc =0.0d0

  ! there is 2 kinds of wavenumbers. The FFT assumes a domain length of 2*pi, then the k
  ! are integers. How should the FFT know otherwise? it just gets a number of datapoints
  ! To get the "real" wavenumber, for instance to compute a derivative, the k has to be scaled
  ! by scalex = 2.d0*pi/xl

  ! the very largest wavenumber occuring in the field is given by
  ! kmax = norm2( (/scalex*nx/2,scaley*ny/2,scalez*nz/2/) )
  ! and corresponds to the corner in wavenumber space. However, for a spectrum, we integrate
  ! over wavenumber shells: therefore, the largest COMPLETE shell is given by the smallest maximum
  ! wavenumber in one direction.
  ! On the other hand, if we do that, \int E(k) does not equal the energy of the field..
!  kmax = minval( (/scalex*dble(nx/2),scaley*dble(ny/2),scalez*dble(nz/2)/) )
  kmax = norm2( (/scalex*dble(nx/2),scaley*dble(ny/2),scalez*dble(nz/2)/) )

  ! spacing in wavenumber space oddly is the scale factors itself
  dk = minval( (/scalex,scaley,scalez/) )

  ! the number of bins for the spectrum: (note that due to symmetry, only half the
  ! range is actually used) (we also return a spectrum array of size (0:nx-1), which contains
  ! zeros out of the valid range)
  nk = nint(kmax / dk)

  if (nx <= nk) then
    call abort(12,"for some reason, we have too many wavenumbers in spectrum..")
  endif

  do iz=ca(1),cb(1)
    kz = wave_z(iz)
    do iy=ca(2),cb(2)
      ky = wave_y(iy)
      do ix=ca(3),cb(3)
        kx = wave_x(ix)
        ! compute 2-norm of wavenumber, scaled to our actual domain size
        kreal = dsqrt( (kx*kx)+(ky*ky)+(kz*kz))
        ! note spectrum omits parts of the data (circle in square problem)
        if (kreal <= kmax) then

          ! then round it so that we can fit it in the corresponding bin, i.e.
          ! for example 0.5 <= k <= 1.5 is K=1 (this is a spherical wavenumber
          ! shell in k-space). this operation can easily be achieved by rounding
          ! kabs to the nearest integer
          ik = nint( kreal/dk )

          ! NOTE what we are computing here is the power spectrum, so a sine-wave 1*sin(x)
          ! has the power 0.25 (or generall A*sin(x) has (A/2)^2 energy)
          if ( ix==0 .or. ix==nx/2 ) then
            sum_ux = dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2)/2.d0
            sum_uy = dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2)/2.d0
            sum_uz = dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)/2.d0
          else
            sum_ux = dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2)
            sum_uy = dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2)
            sum_uz = dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)
          endif

          sum_u = sum_ux + sum_uy + sum_uz
          S_Ekinx_loc(ik) = S_Ekinx_loc(ik) + sum_ux
          S_Ekiny_loc(ik) = S_Ekiny_loc(ik) + sum_uy
          S_Ekinz_loc(ik) = S_Ekinz_loc(ik) + sum_uz
          S_Ekin_loc(ik)  = S_Ekin_loc(ik)  + sum_u
        endif
      enddo
    enddo
  enddo

  ! the returned wavenumber is the middle of the bin: (and not kreal!)
  do ik = 0, nk
    kvec(ik) = dble(ik)*dk
  enddo

  call MPI_ALLREDUCE(S_Ekinx_loc,S_Ekinx,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(S_Ekiny_loc,S_Ekiny,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(S_Ekinz_loc,S_Ekinz,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE(S_Ekin_loc ,S_Ekin ,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

end subroutine compute_spectrum
