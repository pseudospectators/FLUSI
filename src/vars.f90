! the mask and us (the solid Velocity) are for historic reasons global arrays
! that do not appear in the argument lists. to gain flexibility, we move them out
! of "vars". Then newer routines from newer codes, that do not have mask as a global
! can be adapted easier.
module penalization
  implicit none
  ! The mask array.  TODO: move out of shave_vars?
  real(kind=8),dimension (:,:,:),allocatable,save :: mask ! mask function
  integer(kind=2),dimension (:,:,:),allocatable,save :: mask_color ! mask color function
  ! Velocity field inside the solid.  TODO: move out of shave_vars?
  real(kind=8),allocatable,save :: us(:,:,:,:)  ! Velocity in solid
end module penalization


! Variables for pseudospectral simnulations
module vars
  use mpi

  ! include the new timing module code-wide so it can be used everywhere
  use module_timing

  implicit none

  character(len=1),save:: tab ! Fortran lacks a native tab, so we set one up.
  ! Used in params.f90
  integer,parameter :: strlen=120   ! standard string length


 ! Precision of doubles
  integer,parameter :: pr = 8
  integer,parameter :: i8 = 8

  ! Method variables set in the program file:
  character(len=strlen),save :: method ! mhd  or fsi
  character(len=strlen),save :: equation ! navier-stokes; artificial-compressibility
  integer,save :: nf  ! number of linear exponential fields (1 for HYD, 2 for MHD)
  integer,save :: nd  ! number of fields (3 for NS, 6 for MHD)
  integer,save :: neq ! number of fields in u-vector (3 for HYD, 6 for MHD)
  integer,save :: nrw ! number of real work arrays in work
  integer,save :: ncw ! number of complex work arrays in workc
  integer,save :: ng  ! number of ghostpoints (if used)
  integer,save :: nrhs ! number of registers for right-hand sides

  ! MPI and p3dfft variables and parameters
  integer,save :: mpisize, mpirank
  ! Local array bounds
  integer,dimension (1:3),save :: ra,rb,rs,ca,cb,cs
  ! Local array bounds with ghost points
  integer,dimension (1:3),save :: ga,gb
  ! Local array bounds for real arrays for all MPI processes
  integer, dimension (:,:), allocatable, save :: ra_table, rb_table, yz_plane_ranks
  integer, dimension (:,:), allocatable, save :: ca_table, cb_table
  ! for simplicity, store what decomposition we use
  character(len=strlen), save :: decomposition

  ! p3dfft domain decomposition parameters and communicators
  integer,save :: mpicommcart,mpicommy,mpicommz,mpitaskid,mpitasks
  integer,dimension(2),save :: mpidims,mpicoords,mpicommslab
  ! only root rank has this true:
  logical, save :: root=.false.

  real(kind=pr),parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164d0

  real(kind=pr),dimension(:),allocatable,save :: lin ! contains nu and eta

  real(kind=pr), save :: tstart=0.0d0, time_total

  ! Variables set via the parameters file
  real(kind=pr),save :: length, alpha_generic
  real(kind=pr),dimension(1:3),save :: us_fixed

  ! krylov
  integer :: M_max = 20
  real(kind=pr) :: krylov_err_threshold = 1.0e-4

  ! Domain size variables:
  integer,save :: nx,ny,nz
  real(kind=pr),save :: xl,yl,zl,dx,dy,dz,scalex,scaley,scalez

  ! dealiase (done in cal_vis.f90)
  integer,save :: iDealias

  ! Isotropic Turbulence Forcing
  character(len=strlen), save :: forcing_type="none"
  real(kind=pr), save :: kf, eps_forcing

  ! Parameters to set which files are saved and how often:
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask,iSaveMagVorticity
  integer,save :: iSaveMagneticField,iSaveCurrent,iSaveSolidVelocity
  integer,save :: striding=1,idobackup
  character(len=strlen),save :: backup_type = "one-file-backup"
  real(kind=pr),save :: tintegral ! Time between output of integral quantities
  real(kind=pr),save :: tsave ! Time between outpout of entire fields.
  real(kind=pr),save :: tsave_first ! don't save before this time
  ! compute drag force every itdrag time steps and compute unst corrections if
  ! you've told to do so.
  integer,save :: itdrag, unst_corrections
  ! save beam every itbeam time steps
  integer,save :: itbeam
  real(kind=pr),save :: truntime, truntimenext ! Number of hours bet
  real(kind=pr),save :: wtimemax ! Stop after a certain number of hours of wall.
  ! for periodically repeating flows, it may be better to always have only
  ! one set of files on the disk
  character(len=strlen),save :: save_only_one_period, field_precision, naming
  real(kind=pr),save :: tsave_period ! then this is period time
  character(len=strlen),save :: iSaveSpectrae

  ! Time-stepping parameters
  real(kind=pr),save :: tmax
  real(kind=pr),save :: dt_fixed
  real(kind=pr),save :: dt_max=0.d0
  real(kind=pr),save :: cfl, CFL_eta=0.99d0
  integer,save :: nt
  character(len=strlen),save :: iTimeMethodFluid, intelligent_dt

  ! viscosity:
  real(kind=pr),save :: nu

  ! artificial-compressibility
  real(kind=pr),save :: c_0, gamma_p
  integer, save :: acm_sponge
  character(len=strlen) :: acm_inipressure

  ! Initial conditions:
  character(len=strlen),save :: inicond, file_ux, file_uy, file_uz, file_p
  character(len=strlen),save :: file_bx, file_by, file_bz, inicond_spectrum_file
  real(kind=pr),save :: omega1, nu_smoothing

  ! Boundary conditions:
  character(len=strlen),save :: iMask   
  integer(kind=2),save :: endcolor ! the highest color value used
  integer,save :: iMoving,iPenalization
  real(kind=pr),save :: eps
  real(kind=pr),save :: r1,r2,r3 ! Parameters for boundary conditions
  real(kind=pr),save :: pseudoeps, pseudodt, pseudoerrmin, pseudoerrmax
  logical, save :: periodic = .false.

  ! turbulent inlet:
  character(len=strlen),save :: use_turbulent_inlet
  real(kind=pr),save :: rescale
  integer, save :: inlet_thickness

  ! averaging
  character(len=strlen),save :: time_avg, vel_avg, ekin_avg, save_one_only
  character(len=strlen),save :: enstrophy_avg
  ! this is a hack, currently the avg velocity field is global
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk_avg
  ! kinetic energy:
  real(kind=pr),dimension(:,:,:),allocatable :: e_avg
  ! enstrophy:
  real(kind=pr),dimension(:,:,:),allocatable :: Z_avg
  real(kind=pr)::tstart_avg

  ! saving of slices
  character(len=strlen),save :: use_slicing
  integer,save :: itslice, ncache_slices
  integer,save :: slices_to_save(1:4)
  real(kind=pr) :: tslice,tslice_first

  ! switch for vorticity sponge:
  character(len=strlen), save :: iVorticitySponge, iSpongeType
  real(kind=pr), save :: eps_sponge
  integer, save :: sponge_thickness

  ! cavity mask:
  character(len=strlen), save :: iCavity, iChannel
  integer, save :: cavity_size
  ! wall thickness
  real(kind=pr),save :: thick_wall
  ! wall position (solid from pos_wall to pos_wall+thick_wall)
  real(kind=pr),save :: pos_wall

  ! save forces and use unsteady corrections?
  integer, save :: compute_forces

  ! for the iterative FSI coupling, we need these two work
  ! arrays (THIS IS A HACK - TO BE REMOVED)
  complex(kind=pr),allocatable,dimension(:,:,:,:)::nlk_tmp
  complex(kind=pr),dimension(:,:,:,:),allocatable:: uk_old ! TODO: allocate only once

  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle
  real(kind=pr),save :: origin(1:3) = 0.0d0 ! origin of grid (only used in postprocessing, usually the grid starts at 0,0,0)

  ! mean flow control
  real(kind=pr),save :: Uxmean,Uymean,Uzmean, m_fluid, umean_amplitude(1:3)
  real(kind=pr),save :: umean_freq
  character(len=strlen),save :: iMeanFlow_x,iMeanFlow_y,iMeanFlow_z
  ! mean flow startup conditioner (if "dynamic" and mean flow at t=0 is not zero
  ! the forces are singular at the beginning. use the startup conditioner to
  ! avoid large accelerations in mean flow at the beginning)
  character(len=strlen),save :: iMeanFlowStartupConditioner
  real(kind=pr) :: tau_meanflow, T_release_meanflow

  ! parameters for passive scalar advection. NOTE: other parameters (such as
  ! diffusivity) are set for each (of up to 9) passive scalar in a derived datatype
  integer, save :: use_passive_scalar
  integer, save :: n_scalars
  logical, save :: compute_scalar
  character(len=strlen),save :: stop_on_fail

  ! solid model main switch
  character(len=strlen),save :: use_solid_model
  !-----------------------------------------------------------------------------
  ! The derived integral quantities for fluid-structure interactions.
  type Integrals
     real(kind=pr) :: time = 0.d0
     real(kind=pr) :: EKin = 0.d0
     real(kind=pr) :: Dissip = 0.d0
     real(kind=pr) :: Divergence = 0.d0
     real(kind=pr) :: Volume = 0.d0
     real(kind=pr) :: APow = 0.d0
     real(kind=pr) :: IPow = 0.d0
     real(kind=pr) :: penalization_power = 0.d0
     real(kind=pr) :: penalization_power_x = 0.d0
     real(kind=pr) :: penalization_power_y = 0.d0
     real(kind=pr) :: penalization_power_z = 0.d0
     real(kind=pr),dimension(1:3) :: Force = 0.d0
     real(kind=pr),dimension(1:3) :: Force_unst = 0.d0
     real(kind=pr),dimension(1:3) :: Torque = 0.d0
     real(kind=pr),dimension(1:3) :: Torque_unst = 0.d0
  end type Integrals
  !-----------------------------------------------------------------------------

  type(Integrals),save :: GlobalIntegrals

  ! Physical parameters for MHD parameters
  real(kind=pr),save :: eta ! magnetic diffusivity
  real(kind=pr),save :: b0, bc ! Boundary condition parameters
  real(kind=pr),save :: ay ! x*x + ay*y*y -r1*r1 == 0 ?in boundary conditions


!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
! Helper routines for general purpose use throughout the code
!*****************************************************************************
!*****************************************************************************
!*****************************************************************************
  interface in_domain
    module procedure in_domain1, in_domain2
  end interface

  interface on_proc
    module procedure on_proc1, on_proc2
  end interface


!!!!!!!!!!
  contains
!!!!!!!!!!


  !-----------------------------------------------------------------------------
  ! convert degree to radiant
  !-----------------------------------------------------------------------------
  real(kind=pr) function deg2rad(deg)
    implicit none
    real(kind=pr), intent(in) :: deg
    deg2rad=deg*pi/180.d0
    return
  end function

  !-----------------------------------------------------------------------------
  ! radiant to degree
  !-----------------------------------------------------------------------------
  real(kind=pr) function rad2deg(deg)
    implicit none
    real(kind=pr), intent(in) :: deg
    rad2deg=deg*180.d0/pi
    return
  end function

  !-----------------------------------------------------------------------------
  ! cross product of two vectors
  !-----------------------------------------------------------------------------
  function cross(a,b)
    implicit none
    real(kind=pr),dimension(1:3),intent(in) :: a,b
    real(kind=pr),dimension(1:3) :: cross
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
  end function

  !-----------------------------------------------------------------------------
  ! 2-norm length of vectors
  !-----------------------------------------------------------------------------
  function norm2(a)
    implicit none
    real(kind=pr),dimension(1:3),intent(in) :: a
    real(kind=pr) :: norm2
    norm2 = sqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )
  end function

  !---------------------------------------------------------------------------
  ! return periodic index, i.e. if we give ix greater than nx, return
  ! smallest image convention. used, e.g., when computing finite difference
  ! operators or interpolations
  !---------------------------------------------------------------------------
  integer function GetIndex(ix,nx)
    implicit none
    integer, intent (in) ::ix,nx
    integer :: tmp
    tmp=ix
    if (tmp<0) tmp = tmp+nx
    if (tmp>nx-1) tmp = tmp-nx
    GetIndex=tmp
    return
  end function GetIndex
  !---------------------------------------------------------------------------
  integer function per(ix,nx)
    implicit none
    integer, intent (in) ::ix,nx
    integer :: tmp
    tmp=ix
    if (tmp<0) tmp = tmp+nx
    if (tmp<0) tmp = tmp+nx
    if (tmp>nx-1) tmp = tmp-nx
    if (tmp>nx-1) tmp = tmp-nx
    if (nx==1) tmp=0
    per=tmp
    return
  end function per

  !---------------------------------------------------------------------------
  ! abort run, with error code and abort-message
  !---------------------------------------------------------------------------
  subroutine abort(code,msg)
    use mpi
    implicit none
    integer :: mpicode
    integer, intent(in) :: code
    character(len=*), intent(in) :: msg

    write(*,*) msg
    call MPI_abort(MPI_COMM_WORLD,code,mpicode)
  end subroutine abort

  !---------------------------------------------------------------------------
  ! wrapper for NaN checking (this may be compiler dependent)
  !---------------------------------------------------------------------------
  logical function is_nan( x )
    implicit none
    real(kind=pr)::x
    is_nan = .false.
    if (.not.(x.eq.x)) is_nan=.true.
  end function

  !---------------------------------------------------------------------------
  ! check wether real coordinates x are in the domain
  !---------------------------------------------------------------------------
  logical function in_domain1( x )
    implicit none
    real(kind=pr),intent(in)::x(1:3)
    in_domain1 = .false.
    if ( ((x(1)>=0.d0).and.(x(1)<xl)).and.&
         ((x(2)>=0.d0).and.(x(2)<yl)).and.&
         ((x(3)>=0.d0).and.(x(3)<zl)) ) in_domain1=.true.
  end function

  !---------------------------------------------------------------------------
  ! check wether integer coordinates x are in the domain
  !---------------------------------------------------------------------------
  logical function in_domain2( x )
    implicit none
    integer,intent(in)::x(1:3)
    in_domain2 = .false.
    if (  ((x(1)>=0).and.(x(1)<nx-1)).and.&
          ((x(2)>=0).and.(x(2)<ny-1)).and.&
          ((x(3)>=0).and.(x(3)<nz-1)) ) in_domain2=.true.
  end function

  !---------------------------------------------------------------------------
  ! check wether real coordinates x are on this mpi-process
  !---------------------------------------------------------------------------
  logical function on_proc1( x )
    implicit none
    real(kind=pr),intent(in)::x(1:3)
    on_proc1 = .false.
    if (  ((x(1)>=ra(1)*dx).and.(x(1)<=rb(1)*dx)).and.&
          ((x(2)>=ra(2)*dy).and.(x(2)<=rb(2)*dy)).and.&
          ((x(3)>=ra(3)*dz).and.(x(3)<=rb(3)*dz)) ) on_proc1=.true.
  end function

  !---------------------------------------------------------------------------
  ! check wether integer coordinates x are on this mpi-process
  !---------------------------------------------------------------------------
  logical function on_proc2( x )
    implicit none
    integer,intent(in)::x(1:3)
    on_proc2 = .false.
    if ( ((x(1)>=ra(1)).and.(x(1)<=rb(1))).and.&
         ((x(2)>=ra(2)).and.(x(2)<=rb(2))).and.&
         ((x(3)>=ra(3)).and.(x(3)<=rb(3))) ) on_proc2=.true.
  end function

  !---------------------------------------------------------------------------
  ! wrapper for random number generator (this may be compiler dependent)
  !---------------------------------------------------------------------------
  real(kind=pr) function rand_nbr()
    implicit none
    call random_number( rand_nbr )
  end function

  !---------------------------------------------------------------------------
  ! soft startup funtion, is zero until time=time_release, then gently goes to
  ! one during the time period time_tau
  !---------------------------------------------------------------------------
  real(kind=pr) function startup_conditioner(time,time_release,time_tau)
    implicit none
    real(kind=pr), intent(in) :: time,time_release,time_tau
    real(kind=pr) :: t

    t = time-time_release

    if (time <= time_release) then
      startup_conditioner = 0.d0
    elseif ( ( time >time_release ).and.(time<(time_release + time_tau)) ) then
      startup_conditioner =  (t**3)/(-0.5d0*time_tau**3) + 3.d0*(t**2)/time_tau**2
    else
      startup_conditioner = 1.d0
    endif
  end function

  !-----------------------------------------------------------------------------
  ! Condition for output conditions.
  ! return true after tfrequ time units or itfrequ time steps or if we're at the
  ! and of the simulation
  !-----------------------------------------------------------------------------
  logical function time_for_output( time, dt, it, tfrequ, ifrequ, tmax, tfirst )
    implicit none
    real(kind=pr), intent(in) :: time, dt, tfrequ, tfirst
    real(kind=pr), intent(in) :: tmax ! final time (if we save at the end of run)
    integer, intent(in) :: it ! time step counter
    integer, intent(in) :: ifrequ ! save every ifrequ time steps

    real(kind=pr) :: tnext1, tnext2

    time_for_output = .false.

    ! we never save before tfirst
    if (time<tfirst) return

    if (intelligent_dt=="yes") then
      ! with intelligent time stepping activated, the time step is adjusted not
      ! to pass by tsave,tintegral,tmax,tslice
      ! this is the next instant we want to save
      tnext1 = dble(ceiling(time/tfrequ))*tfrequ
      tnext2 = dble(floor  (time/tfrequ))*tfrequ
      ! please note that the time actually is very close to the next instant we
      ! want to save. however, it may be slightly less or larger. therefore, we
      ! cannot just check (time-tnext), since tnext may be wrong
      if ((abs(time-tnext1)<=1.0d-6).or.(abs(time-tnext2)<=1.0d-6).or.&
          (modulo(it,ifrequ)==0).or.(abs(time-tmax)<=1.0d-6)) then
        time_for_output = .true.
      endif
    else
      ! without intelligent time stepping, we save output when we're close enough
      if ( (modulo(time,tfrequ)<dt).or.(modulo(it,ifrequ)==0).or.(time==tmax) ) then
        time_for_output = .true.
      endif
    endif
  end function

  !-----------------------------------------------------------------------------
  ! given a point x, check if it lies in the computational domain centered at zero
  ! (note: we assume [-xl/2...+xl/2] size this is useful for insects )
  !-----------------------------------------------------------------------------
  function periodize_coordinate(x_glob, box)
    real(kind=pr),intent(in) :: x_glob(1:3), box(1:3)
    real(kind=pr),dimension(1:3) :: periodize_coordinate

    periodize_coordinate = x_glob

    if (periodic) then
      if (x_glob(1)<-box(1)/2.0) periodize_coordinate(1)=x_glob(1)+box(1)
      if (x_glob(2)<-box(2)/2.0) periodize_coordinate(2)=x_glob(2)+box(2)
      if (x_glob(3)<-box(3)/2.0) periodize_coordinate(3)=x_glob(3)+box(3)

      if (x_glob(1)>box(1)/2.0) periodize_coordinate(1)=x_glob(1)-box(1)
      if (x_glob(2)>box(2)/2.0) periodize_coordinate(2)=x_glob(2)-box(2)
      if (x_glob(3)>box(3)/2.0) periodize_coordinate(3)=x_glob(3)-box(3)
    endif

  end function


!!!!!!!!!!!!!!!
end module vars
!!!!!!!!!!!!!!!
