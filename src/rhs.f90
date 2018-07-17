! Wrapper for computing the nonlinear source term for Navier-Stokes/MHD
subroutine cal_nlk(time,it,nlk,uk,u,vort,work,workc,press,scalars,scalars_rhs,Insect,beams)
  use vars
  use p3dfft_wrapper
  use solid_model
  use module_insects
  use passive_scalar_module
  implicit none

  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(in) :: time
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
  real(kind=pr) :: t1,t0
  integer, intent(in) :: it
  t0 = MPI_wtime()

  !-----------------------------------------------------------------------------
  ! If the mask is time-dependend, we create it here
  ! 26/01/2015 Moved the mask routine here to RHS, since this is where it con-
  ! ceptually belongs. RK2 and RK4 need to create the mask several times per time
  ! step.
  !-----------------------------------------------------------------------------
  if ((iMoving==1).and.(iPenalization==1).and.(iTimeMethodFluid/="FSI_AB2_iteration")) then
    ! for the iterative FSI schemes, the mask is created in fluidtimestep
    ! (so iTimeMethodFluid==FSI_AB2_iteration skips mask generation)
    call create_mask( time, Insect, beams )
  endif

  select case(method)
  case("fsi")
    !---------------------------------------------------------------------------
    ! FSI case. note projection is outsourced and performed here.
    !---------------------------------------------------------------------------

    select case (equation)
    case ("navier-stokes")
      ! compute source-terms, *not* divergence-free
      t1 = MPI_wtime()
      call cal_nlk_fsi(time,it,nlk,uk,u,vort,work,workc,Insect)
      time_nlk2 = time_nlk2 + MPI_wtime() - t1

      t1 = MPI_wtime()
      ! if we compute active FSI (with flexible obstacles), we need the pressure
      if (use_solid_model=="yes") then
        call pressure( nlk,workc(:,:,:,1) )
        ! transform it to phys space (note "press" has ghostpoints, cut them here)
        call ifft( ink=workc(:,:,:,1), outx=press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
      endif
      ! project the right hand side to the incompressible manifold
      call add_grad_pressure(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3))
      ! for global performance measurement
      time_p = time_p + MPI_wtime() - t1

      !---------------------------------------------------------------------------
      ! passive scalar. the new module uses finite differences and evolves up to 9
      ! different colors. for FD, no work arrays are required.
      !---------------------------------------------------------------------------
      t1 = MPI_wtime()
      if ((use_passive_scalar==1).and.(compute_scalar)) then
        call cal_nlk_scalar(time,it,u,scalars,scalars_rhs)
      endif
      time_nlk_scalar = time_nlk_scalar + MPI_wtime() - t1

    case ("artificial-compressibility")
      t1 = MPI_wtime()
      call rhs_acm(time,it,nlk,uk,u,vort,work,workc,Insect)
      time_nlk2 = time_nlk2 + MPI_wtime() - t1

    case default
      call abort(129388,"For some reason, the equation is unkown. Set navier-stokes or artificial-compressibility")

    end select

  case("mhd")
     !--------------------------------------------------------------------------
     ! MHD case
     !--------------------------------------------------------------------------
     call cal_nlk_mhd(nlk,uk,u,vort)

  case default
     call abort(1,"Error! Unknown method in cal_nlk")
  end select

  time_rhs = time_rhs + MPI_wtime() - t0
end subroutine cal_nlk


!-------------------------------------------------------------------------------
! Compute the nonlinear source term of the Navier-Stokes equation,
! including penalty term, in Fourier space. Seven real-valued
! arrays are required for working memory. The term in NLK reads
! nlk = omega x u - chi/eta * (u-us) - sponge
! Input:
!       time: guess what!
!       uk: vector field of velocity in Fourier space, holding the 3 velocity
!           components, if present
! Output:
!       nlk:  The right hand side of penalized Navier-Stokes in Fourier space,
!             ie the NL term, penalty term, sponge term (NOT THE PRESSURE)
!       vort: work array, can be reused immediatly (is free after this routine)
!       u:    work array, contains the velocity in phys space this is reused in
!             the caller FluidTimestep to adjust dt
!       work: work array (real)
!       workc: work array (cmplx), for sponge and/or passive scalar
!
! NOTES:
!       16 Jul 2014: This routine excludes the pressure. the field NLK is
!          NOT DIVERGENCE FREE. The projection takes place in the
!          caller "cal_nlk"
!       09 Aug 2014: enable dynamic mean flow forcing. this solves the eqn
!          d(u_mean)/dt = F / m_fluid
!          This is useful if one wants the mean flow to be independend of the
!          domain size. It accelerates a given mass of fluid, which is fixed.
!          Solving this eqn requires to compute the sum of the penalty term
!          in this routine, which is dead weight when not using this function.
!-------------------------------------------------------------------------------
subroutine cal_nlk_fsi(time,it,nlk,uk,u,vort,work,workc, Insect)
  use mpi
  use p3dfft_wrapper
  use vars
  use module_insects
  use basic_operators
  use penalization ! mask array etc

  implicit none

  real(kind=pr),intent (in) :: time
  integer, intent(in) :: it
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  type(diptera),intent(inout) :: Insect
  real(kind=pr) :: t0,t1,ux,uy,uz,vorx,vory,vorz,chi,usx,usy,usz
  real(kind=pr) :: fx,fy,fz,fx1,fy1,fz1
  real(kind=pr) :: soft_startup, dt, eps_inv
  integer :: ix,iz,iy,mpicode
  t0 = MPI_wtime()
  fx=0.d0; fy=0.d0; fz=0.d0

  if (eps<1.0d-11) call abort(77363800, "Value of eps is very small, maybe even zero.")
  eps_inv = 1.0_pr / eps


  !-----------------------------------------------------------------------------
  !-- Calculate velocity in physical space
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  call ifft3 (outx=u, ink=uk)
  time_u = time_u + MPI_wtime() - t1

  !-----------------------------------------------------------------------------
  !-- Compute vorticity
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  ! nlk is temporarily used for vortk
  call curl ( ink=uk, outk=nlk )
  ! transform it to physical space
  call ifft3 ( ink=nlk, outx=vort )
  time_vor = time_vor + MPI_wtime() - t1

  !-----------------------------------------------------------------------------
  ! compute time-averaged enstrophy Z_avg
  !-----------------------------------------------------------------------------
  ! we do this here since we happen to have the voriticy and want to save FFTs
  if ((time_avg=="yes").and.(enstrophy_avg=="yes").and.(time>=tstart_avg)) then
    ! we need to know here what the time step will be
    call adjust_dt(time,u,dt)
    ! compute incremental avg
    Z_avg = ( (vort(:,:,:,1)**2 + vort(:,:,:,2)**2 + vort(:,:,:,3)**2)*dt  &
    + (time-tstart_avg)*Z_avg ) / ( (time-tstart_avg)+dt )
  endif

  !-----------------------------------------------------------------------------
  !-- vorticity sponge term
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  if (iVorticitySponge == "yes") then
    call vorticity_sponge( vort, work(:,:,:,1), workc, Insect )
  endif
  time_sponge = time_sponge + MPI_wtime() - t1

  !-----------------------------------------------------------------------------
  !-- Non-Linear terms
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ! local loop variables
        ux   = u(ix,iy,iz,1)
        uy   = u(ix,iy,iz,2)
        uz   = u(ix,iy,iz,3)
        vorx = vort(ix,iy,iz,1)
        vory = vort(ix,iy,iz,2)
        vorz = vort(ix,iy,iz,3)

        ! local variables for penalization. NEW: since 07/2018, we divide the mask
        ! by eps here AND ONLY HERE! it is much easier to understand then. there may be
        ! a very slight performance penalty if the obstacle does not move.
        chi = mask(ix,iy,iz) * eps_inv
        usx = us(ix,iy,iz,1)
        usy = us(ix,iy,iz,2)
        usz = us(ix,iy,iz,3)

        ! compute sum of penalty term while computing it
        fx = fx + chi*(ux-usx)
        fy = fy + chi*(uy-usy)
        fz = fz + chi*(uz-usz)

        ! we overwrite the vorticity with the NL terms in phys space
        ! note this is indeed -(vor x u) (negative sign)
        vort(ix,iy,iz,1) = uy*vorz - uz*vory -chi*(ux-usx)
        vort(ix,iy,iz,2) = uz*vorx - ux*vorz -chi*(uy-usy)
        vort(ix,iy,iz,3) = ux*vory - uy*vorx -chi*(uz-usz)
      enddo
    enddo
  enddo
  ! to Fourier space
  call fft3( inx=vort,outk=nlk )
  time_curl = time_curl + MPI_wtime() - t1

  !-----------------------------------------------------------------------------
  ! add sponge term
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  if (iVorticitySponge == "yes") then
    nlk(:,:,:,1:3) = nlk(:,:,:,1:3) + workc(:,:,:,1:3)
  endif
  time_sponge = time_sponge + MPI_wtime() - t1

  !-----------------------------------------------------------------------------
  ! dynamic mean flow forcing (fix fluid mass manually, domain-independent)
  ! The Mean Flow is governed by the zeroth Fourier mode of the RHS. Note this
  ! is independent of the pressure (since it has vannishing spatial avg).
  ! If the meanflow at t=0 is not 0, the startup singularity in the forces
  ! causes problems. for this case, we use the startup conditioner to keep
  ! the meanflow const until T_release_meanflow, then gently turning it on
  ! during the time tau_meanflow
  !-----------------------------------------------------------------------------
  if(iMeanFlow_x=="dynamic".or.iMeanFlow_y=="dynamic".or.iMeanFlow_z=="dynamic") then
    ! integral forces:
    fx = fx*dx*dy*dz
    fy = fy*dx*dy*dz
    fz = fz*dx*dy*dz
    call MPI_ALLREDUCE ( fx,fx1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
    call MPI_ALLREDUCE ( fy,fy1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
    call MPI_ALLREDUCE ( fz,fz1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)

    ! startup conditioner (See vars.f90)
    if (iMeanFlowStartupConditioner=="yes") then
      soft_startup = startup_conditioner(time,T_release_meanflow,tau_meanflow)
      fx1 = fx1*soft_startup
      fy1 = fy1*soft_startup
      fz1 = fz1*soft_startup
    endif

    ! fixing the fluid mass means modifying the zero mode of RHS term
    if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
      if(iMeanFlow_x=="dynamic") nlk(0,0,0,1) = -fx1 / m_fluid
      if(iMeanFlow_y=="dynamic") nlk(0,0,0,2) = -fy1 / m_fluid
      if(iMeanFlow_z=="dynamic") nlk(0,0,0,3) = -fz1 / m_fluid
    endif
  endif

  !-----------------------------------------------------------------------------
  ! forcing that accelerates mean flow from rest to unity and keeps it there
  !-----------------------------------------------------------------------------
  if(iMeanFlow_x=="accelerate_to_unity".or.iMeanFlow_y=="accelerate_to_unity".or.iMeanFlow_z=="accelerate_to_unity") then
    if (ca(1) == 0 .and. ca(2) == 0 .and. ca(3) == 0) then
      fx = max(0.d0,1.d0-dreal(uk(0,0,0,1)))
      fy = max(0.d0,1.d0-dreal(uk(0,0,0,2)))
      fz = max(0.d0,1.d0-dreal(uk(0,0,0,3)))
      if(iMeanFlow_x=="accelerate_to_unity") nlk(0,0,0,1) = nlk(0,0,0,1) + fx
      if(iMeanFlow_y=="accelerate_to_unity") nlk(0,0,0,2) = nlk(0,0,0,2) + fy
      if(iMeanFlow_z=="accelerate_to_unity") nlk(0,0,0,3) = nlk(0,0,0,3) + fz
    endif
  endif

  !---------------------------------------------------------------------------
  ! Add explicit diffusion term here, if RK4 is used (other time steppers have
  ! integrating factors and thus implicit diffusion)
  !---------------------------------------------------------------------------
  if (iTimeMethodFluid=="RK4") then
    call add_explicit_diffusion(uk,nlk)
  endif

  !-----------------------------------------------------------------------------
  ! Add forcing term for isotropic turbulence, if used
  !-----------------------------------------------------------------------------
  if (forcing_type/="none") then
    call add_forcing_term(time,uk,nlk,vort)
  endif


  time_nlk = time_nlk + MPI_wtime() - t0
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

  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)
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
! Compute the pressure given the velocity field.
!-------------------------------------------------------------------------------
! Note the pressure is computed from the divergence of the RHS of Navier-Stokes.
! This implies that you actually need to compute the RHS in order to compute
! the pressure from it - thus, this is (almost) the price of a Navier-Stokes
! step.
! The good way to compute the pressure is to do it when you compute the RHS anyways
! but in some situations, this is not possible:
!   - in the semi-implicit FSI scheme (advance fluid, then get new pressure, then advance solid)
!   - in postprocessing
! For these reasons, this routine exists.
! IMPORTANT: NOTE: the routine does *NOT* create mask and us fields. Why? Because in the FSI
! semiimplicit scheme, we compute the new pressure using the OLD mask (since the beam is not
! yet advanced to the new time level).
! In postprocessing, be sure to create the mask before calling this routine!!
!-------------------------------------------------------------------------------
subroutine pressure_from_uk_use_existing_mask(time,u,uk,nlk,vort,work,workc,press,Insect)
  use mpi
  use p3dfft_wrapper
  use vars
  use module_insects
  implicit none

  real(kind=pr),intent(in) :: time
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  type(diptera),intent(inout)::Insect

  ! first, compute the source terms: non-linear operator and penalization term.
  ! NOTE: we do not create the mask here (nor in cal_nlk_fsi!) see note above!
  ! This implies that if you do not have the mask at the right time or not set at
  ! all, the result is wrong
  call cal_nlk_fsi(time,0,nlk,uk,u,vort,work,workc,Insect)

  ! compute the pressure from source-terms
  call pressure(nlk,workc(:,:,:,1))

  ! transform pressure back to phys. space (NOTE: press is an ghost-point array,
  ! therefore you need to copy only the ra:rb parts)
  call ifft(ink=workc(:,:,:,1), outx=press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  ! as we employ the rotational formulation for the nonlinear term, the pressure
  ! is PI and not p. Return p by sustracting kinetic energy:
  press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) = press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) &
  - 0.5d0*(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)
end subroutine


!-------------------------------------------------------------------------------
! Add the explicit diffusion term to the right hand side. This is only necessary
! for RK4 scheme currently (01/2015) since all other time steppers have integrating
! factors.
!-------------------------------------------------------------------------------
subroutine add_explicit_diffusion(uk,nlk)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout):: uk (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz,k2

  do iz=ca(1),cb(1)
    !-- wavenumber in z-direction
    kz = wave_z(iz)
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)
        k2 = kx*kx + ky*ky + kz*kz
        nlk(iz,iy,ix,1) = nlk(iz,iy,ix,1) - nu * k2 * uk(iz,iy,ix,1)
        nlk(iz,iy,ix,2) = nlk(iz,iy,ix,2) - nu * k2 * uk(iz,iy,ix,2)
        nlk(iz,iy,ix,3) = nlk(iz,iy,ix,3) - nu * k2 * uk(iz,iy,ix,3)
      enddo
    enddo
  enddo
end subroutine add_explicit_diffusion


! forcing terms for Homogeneous Isotropic Turbulence (HIT), added directly to NLK
! The forcing injects energy into some wavenumbers, trying to keep the total
! energy constant by compensating for viscous losses.
subroutine add_forcing_term(time,uk,nlk,work)
  use vars
  use p3dfft_wrapper
  use basic_operators
  implicit none

  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(inout):: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout):: uk (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  real(kind=pr), dimension(0:nx-1) :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin, kvec
  real(kind=pr) :: kx,ky,kz,kreal,factor,epsilon,epsilon_loc,k2
  real(kind=pr) :: E
  integer :: ix,iy,iz, k, mpicode

  select case(forcing_type)
  case ("machiels")
    ! Forcing by Machiels PRL1997 ("Predictability of Small-scale Motion in Isotropic...")
    ! This forcing aims at specifying the dissipation rate epsilon, which is
    ! directly related to the kolmogorov scales. The desired value is set in
    ! eps_forcing, which is read from the parameter file.

    ! get spectrum (on all procs, MPI_ALLREDUCE)
    call compute_spectrum(time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
    factor = eps_forcing / (2.d0*(S_Ekin( int(kf) )+S_Ekin( int(kf)+1 )))

    ! force wavenumber shells
    do ix=ca(3),cb(3)
      kx=wave_x(ix)
      do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
          kz=wave_z(iz)
          ! compute 2-norm of wavenumber
          kreal = dsqrt( (kx*kx)+(ky*ky)+(kz*kz) )

          ! forcing lives around this shell
          if ( (kreal>=kf-0.5d0) .and. (kreal<=kf+1.5d0) ) then
            nlk(iz,iy,ix,1) = nlk(iz,iy,ix,1) + uk(iz,iy,ix,1)*factor
            nlk(iz,iy,ix,2) = nlk(iz,iy,ix,2) + uk(iz,iy,ix,2)*factor
            nlk(iz,iy,ix,3) = nlk(iz,iy,ix,3) + uk(iz,iy,ix,3)*factor
          endif

        enddo
      enddo
    enddo

  case ("kaneda")
    ! get spectrum (on all procs, MPI_ALLREDUCE)
    call compute_spectrum(time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)

    ! compute current dissipation rate. There are several ways to do this; one
    ! is computing the enstrophy in x-space, using the vorticity, but we can also
    ! stay with the velocity in k-space and integrate k^2 * E(k)
    ! However, this is exact only if E(k) is not averaged over the wavenumber shell
    ! as it is done when computing the spectrum.
    call dissipation_rate(uk, epsilon)

    ! The actual forcing term follows. We now have the current dissipation rate
    ! epsilon, which we would like to be balanced by the energy input through
    ! the forcing. The forcing is c*uk (where c="factor" here). Note the dissipation
    ! rate is divided by the energy of the first two wavenumber shells
    factor = epsilon / (2.d0*(S_Ekin( int(kf) )+S_Ekin( int(kf)+1 )))

    ! force wavenumber shells
    do ix=ca(3),cb(3)
      kx=wave_x(ix)
      do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
          kz=wave_z(iz)
          ! compute 2-norm of wavenumber
          kreal = dsqrt( (kx*kx)+(ky*ky)+(kz*kz) )

          ! forcing lives around this shell
          if ( (kreal>=kf-0.5d0) .and. (kreal<=kf+1.5d0) ) then
            nlk(iz,iy,ix,1) = nlk(iz,iy,ix,1) + uk(iz,iy,ix,1)*factor
            nlk(iz,iy,ix,2) = nlk(iz,iy,ix,2) + uk(iz,iy,ix,2)*factor
            nlk(iz,iy,ix,3) = nlk(iz,iy,ix,3) + uk(iz,iy,ix,3)*factor
          endif

        enddo
      enddo
    enddo

  case default
    call abort(99929,"unknown HIT forcing method..")
  end select

  ! enforce hermitian symmetry the lazy ways
  call ifft3( ink=nlk, outx=work )
  call fft3( inx=work, outk=nlk  )

end subroutine add_forcing_term

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

  imag = dcmplx(0.d0,1.d0)

  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)
           k2=kx*kx + ky*ky + kz*kz

           if (k2 .ne. 0.0d0) then
              nlx=nlk1(iz,iy,ix)
              nly=nlk2(iz,iy,ix)
              nlz=nlk3(iz,iy,ix)

              ! qk is the Fourier coefficient of thr pressure
              qk=(kx*nlx + ky*nly + kz*nlz)/k2
              ! add the gradient to the non-linear terms
              nlk1(iz,iy,ix)=nlx - kx*qk
              nlk2(iz,iy,ix)=nly - ky*qk
              nlk3(iz,iy,ix)=nlz - kz*qk
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
  use mpi
  use vars
  use p3dfft_wrapper
  use basic_operators
  use penalization ! mask array etc
  implicit none

  complex(kind=pr),intent(inout) ::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
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
  do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
          do ix=ra(1),rb(1)
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
              m=mask(ix,iy,iz) / eps
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
  call curl(nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6))

  ! Penalization for B-field:
  if(iPenalization == 1) then
     do i=4,nd
        wj(:,:,:,4)=-(mask/eps)*(ub(:,:,:,i) - us(:,:,:,i))
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


! Render the input field divergence-free via a Helmholtz
! decomposition. The zero-mode is left untouched. Not this is exactly what is done
! in add_grad_pressure, so this is just an alias.
subroutine div_field_nul(fx,fy,fz)
  use vars
  implicit none

  complex(kind=pr), intent(inout) :: fx(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fy(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr), intent(inout) :: fz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  call add_grad_pressure(fx,fy,fz)

end subroutine div_field_nul



subroutine rhs_acm(time, it, nlk, uk, u, vort, work, workc, Insect)
  use mpi
  use p3dfft_wrapper
  use vars
  use module_insects
  use basic_operators
  use penalization ! mask array etc

  implicit none

  real(kind=pr),intent (in) :: time
  integer, intent(in) :: it
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  type(diptera),intent(inout) :: Insect
  real(kind=pr) :: t0,t1,uy,uz,chi,usy,usz, uy_dy, uy_dz, uz_dy, uz_dz, p
  real(kind=pr) :: kx, ky, kz, k2, eps_inv
  complex(kind=pr) :: imag   ! imaginary unit
  integer :: ix,iz,iy,mpicode
  t0 = MPI_wtime()

  nlk = 0.0_pr
  ! reset ux (we assume you deal with 2d flows...)
  uk(:,:,:,1) = 0.0_pr

  if (nx/=1) call abort(33427, "artificial-compressibility only for 2d implemented")
  ! no integrating factor
  if (iTimeMethodFluid/="RK4")  call abort(33427, "artificial-compressibility only for RK4 implemented")
  if (ncw<4) call abort(77283,"acm: not enough cmplx work arrays")
  if (nrw<4) call abort(77283,"acm: not enough real work arrays")

  eps_inv = 1.0_pr / eps

  !-----------------------------------------------------------------------------
  !-- Calculate velocity in physical space
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  call ifft(outx=u(:,:,:,2), ink=uk(:,:,:,2))
  call ifft(outx=u(:,:,:,3), ink=uk(:,:,:,3))

  ! derivatives (requires for non-linear term in convection formulation)
  imag = dcmplx(0.d0,1.d0)
  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)
           workc(iz,iy,ix,1) = imag*ky*uk(iz,iy,ix,2) ! work1 = uy_dy
           workc(iz,iy,ix,2) = imag*kz*uk(iz,iy,ix,2) ! work2 = uy_dz
           workc(iz,iy,ix,3) = imag*ky*uk(iz,iy,ix,3) ! work3 = uz_dy
           workc(iz,iy,ix,4) = imag*kz*uk(iz,iy,ix,3) ! work4 = uz_dz
       enddo
    enddo
  enddo

  call ifft(ink=workc(:,:,:,1), outx=work(:,:,:,1)) ! work1 = uy_dy
  call ifft(ink=workc(:,:,:,2), outx=work(:,:,:,2)) ! work2 = uy_dz
  call ifft(ink=workc(:,:,:,3), outx=work(:,:,:,3)) ! work3 = uz_dy
  call ifft(ink=workc(:,:,:,4), outx=work(:,:,:,4)) ! work4 = uz_dz

  time_u = time_u + MPI_wtime() - t1
  !-----------------------------------------------------------------------------
  !-- Non-Linear terms
  !-----------------------------------------------------------------------------
  t1 = MPI_wtime()
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ! local loop variables
        uy = u(ix,iy,iz,2)
        uz = u(ix,iy,iz,3)

        ! derivatives
        uy_dy = work(ix,iy,iz,1)
        uy_dz = work(ix,iy,iz,2)
        uz_dy = work(ix,iy,iz,3)
        uz_dz = work(ix,iy,iz,4)

        ! local variables for penalization
        chi = mask(ix,iy,iz) * eps_inv
        usy = us(ix,iy,iz,2)
        usz = us(ix,iy,iz,3)

        ! we overwrite the vorticity with the NL terms in phys space
        ! note this is indeed -(vor x u) (negative sign)
        vort(ix,iy,iz,2) = -uy*uy_dy -uz*uy_dz -chi*(uy-usy)
        vort(ix,iy,iz,3) = -uy*uz_dy -uz*uz_dz -chi*(uz-usz)
      enddo
    enddo
  enddo

  if (acm_sponge==1) then
    ! get pressure in x space
    call ifft(ink=uk(:,:,:,4), outx=work(:,:,:,1) )

    ! penalization term for pressure in sponge areas
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          ! local loop variables
          p = work(ix,iy,iz,1)
          chi = mask(ix,iy,iz) * eps_inv
          ! color 0 is sponges
          if (mask_color(ix,iy,iz)==0) work(ix,iy,iz,2) = -chi*p
        enddo
      enddo
    enddo

    ! to k-space
    call fft( inx=work(:,:,:,2),outk=nlk(:,:,:,4) )
  endif

  ! to Fourier space
  call fft( inx=vort(:,:,:,2),outk=nlk(:,:,:,2) )
  call fft( inx=vort(:,:,:,3),outk=nlk(:,:,:,3) )
  time_curl = time_curl + MPI_wtime() - t1

  ! add remaining terms in fourier space (divergence, grad, laplace)
  imag = dcmplx(0.d0, 1.d0)
  do ix=ca(3),cb(3)
     kx=wave_x(ix)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do iz=ca(1),cb(1)
           kz=wave_z(iz)
           k2 = kx*kx + ky*ky + kz*kz

           ! add pressure gradient and laplace
           nlk(iz,iy,ix,2) = nlk(iz,iy,ix,2) - imag*ky*uk(iz,iy,ix,4) - nu*k2*uk(iz,iy,ix,2)
           nlk(iz,iy,ix,3) = nlk(iz,iy,ix,3) - imag*kz*uk(iz,iy,ix,4) - nu*k2*uk(iz,iy,ix,3)

           ! rhs for pressure (divergence)
           nlk(iz,iy,ix,4) = nlk(iz,iy,ix,4) -(c_0**2)*imag*(ky*uk(iz,iy,ix,2) + kz*uk(iz,iy,ix,3)) -gamma_p*uk(iz,iy,ix,4)
       enddo
    enddo
  enddo

  ! rhs(ix,iy,1) = -u(ix,iy)*u_dx - v(ix,iy)*u_dy - p_dx + nu*(u_dxdx + u_dydy) + penalx
  ! rhs(ix,iy,2) = -u(ix,iy)*v_dx - v(ix,iy)*v_dy - p_dy + nu*(v_dxdx + v_dydy) + penaly
  ! rhs(ix,iy,3) = -(c_0**2)*div_U - gamma*p(ix,iy)

  time_nlk = time_nlk + MPI_wtime() - t0
end subroutine rhs_acm
