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
  implicit none

  character(len=1),save:: tab ! Fortran lacks a native tab, so we set one up.
  ! Used in params.f90
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file
  integer,parameter :: strlen=80   ! standard string length

  ! Precision of doubles
  integer,parameter :: pr = 8
  integer,parameter :: i8 = 8

  ! Method variables set in the program file:
  character(len=strlen),save :: method ! mhd  or fsi
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
  ! for simplicity, store what decomposition we use
  character(len=strlen), save :: decomposition

  ! p3dfft domain decomposition parameters and communicators
  integer,save :: mpicommcart,mpicommy,mpicommz,mpitaskid,mpitasks
  integer,dimension(2),save :: mpidims,mpicoords,mpicommslab
  ! only root rank has this true:
  logical, save :: root=.false.

  real(kind=pr),save :: pi ! 3.14....

  real(kind=pr),dimension(:),allocatable,save :: lin ! contains nu and eta

  ! Vabiables timing statistics.  Global to simplify syntax.
  real(kind=pr),save :: time_fft,time_ifft,time_vis,time_mask,time_nlk2
  real(kind=pr),save :: time_vor,time_curl,time_p,time_nlk,time_u,tslices
  real(kind=pr),save :: time_bckp,time_save,time_total,time_fluid,time_nlk_fft
  real(kind=pr),save :: time_sponge,time_insect_head,time_insect_body, time_scalar
  real(kind=pr),save :: time_insect_eye,time_insect_wings, time_insect_vel
  real(kind=pr),save :: time_solid, time_drag, time_surf, time_LAPACK
  real(kind=pr),save :: time_hdf5,time_integrals,time_rhs,time_nlk_scalar,tstart

  ! Variables set via the parameters file
  real(kind=pr),save :: length, alpha_generic

  ! Domain size variables:
  integer,save :: nx,ny,nz
  real(kind=pr),save :: xl,yl,zl,dx,dy,dz,scalex,scaley,scalez

  ! FIXME: please document
  integer,save :: iDealias

  ! Isotropic Turbulence Forcing
  character(len=strlen), save :: forcing_type
  real(kind=pr), save :: kf, eps_forcing


  ! Parameters to set which files are saved and how often:
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask,iSaveMagVorticity
  integer,save :: idobackup = 1
  integer,save :: striding = 1
  integer,save :: iSaveXMF !directly write *.XMF files (1) or not (0)
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
  character(len=strlen),save :: save_only_one_period, field_precision
  real(kind=pr),save :: tsave_period ! then this is period time
  character(len=strlen),save :: iSaveSpectrae

  ! Time-stepping parameters
  real(kind=pr),save :: tmax
  real(kind=pr),save :: dt_fixed
  real(kind=pr),save :: dt_max=0.d0
  real(kind=pr),save :: cfl
  integer,save :: nt
  character(len=strlen),save :: iTimeMethodFluid, intelligent_dt

  ! Physical parameters:
  real(kind=pr),save :: nu

  ! Initial conditions:
  character(len=strlen),save :: inicond, file_ux,file_uy,file_uz
  character(len=strlen),save :: file_bx,file_by,file_bz
  real(kind=pr),save :: omega1, nu_smoothing


  ! Boundary conditions:
  character(len=strlen),save :: iMask
  integer,save :: iMoving,iPenalization
  real(kind=pr),save :: eps
  real(kind=pr),save :: r1,r2,r3 ! Parameters for boundary conditions
  real(kind=pr),save :: pseudoeps, pseudodt, pseudoerrmin, pseudoerrmax

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


  interface abort
    module procedure abort1, abort2
  end interface

  interface in_domain
    module procedure in_domain1, in_domain2
  end interface

  interface on_proc
    module procedure on_proc1, on_proc2
  end interface

  contains

    ! Routines to allocate real and complex arrays using the standard
    ! dimensions.
    subroutine allocreal(u)
      implicit none
      real(kind=pr),dimension(:,:,:),allocatable :: u
      allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
    end subroutine allocreal

    subroutine alloccomplex(u)
      implicit none
      complex(kind=pr),dimension(:,:,:),allocatable :: u
      allocate(u(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
    end subroutine alloccomplex

    subroutine allocrealnd(u)
      implicit none
      real(kind=pr),dimension(:,:,:,:),allocatable :: u
      allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
    end subroutine allocrealnd

    subroutine alloccomplexnd(u)
      implicit none
      complex(kind=pr),dimension(:,:,:,:),allocatable :: u
      allocate(u(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd))
    end subroutine alloccomplexnd
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
    subroutine abort1
      use mpi
      implicit none
      integer :: mpicode

      if (mpirank==0) write(*,*) "Killing run..."
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine abort1
    !---------------------------------------------------------------------------
    subroutine abort2(msg)
      use mpi
      implicit none
      integer :: mpicode
      character(len=*), intent(in) :: msg

      if (mpirank==0) write(*,*) "Killing run..."
      if (mpirank==0) write(*,*) msg
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine abort2
    !---------------------------------------------------------------------------
    ! wrapper for NaN checking (this may be compiler dependent)
    logical function is_nan( x )
      implicit none
      real(kind=pr)::x
      is_nan = .false.
      if (.not.(x.eq.x)) is_nan=.true.
    end function
    !---------------------------------------------------------------------------
    ! check wether real coordinates x are in the domain
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
    real(kind=pr) function rand_nbr()
      implicit none
      call random_number( rand_nbr )
    end function
    !---------------------------------------------------------------------------
    ! soft startup funtion, is zero until time=time_release, then gently goes to
    ! one during the time period time_tau
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
end module vars


! Variables for fsi simulations
module fsi_vars
  use vars
  implicit none

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
  integer,save :: iSaveSolidVelocity
  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle

  ! mean flow control
  real(kind=pr),save :: Uxmean,Uymean,Uzmean, m_fluid
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
     real(kind=pr) :: time
     real(kind=pr) :: EKin
     real(kind=pr) :: Dissip
     real(kind=pr) :: Divergence
     real(kind=pr) :: Volume
     real(kind=pr) :: APow
     real(kind=pr) :: IPow
     real(kind=pr) :: penalization_power
     real(kind=pr) :: penalization_power_x
     real(kind=pr) :: penalization_power_y
     real(kind=pr) :: penalization_power_z
     real(kind=pr),dimension(1:3) :: Force
     real(kind=pr),dimension(1:3) :: Force_unst
     real(kind=pr),dimension(1:3) :: Torque
     real(kind=pr),dimension(1:3) :: Torque_unst
  end type Integrals
  !-----------------------------------------------------------------------------
  ! derived datatype for time
  type timetype
    real(kind=pr) :: time
    real(kind=pr) :: dt1
    real(kind=pr) :: dt0
    integer :: it
    integer :: n0
    integer :: n1
  end type timetype
  !-----------------------------------------------------------------------------

  type(Integrals),save :: GlobalIntegrals

  contains
  !-----------------------------------------------------------------------------

  ! this function simplifies my life with the insects
  real(kind=pr) function deg2rad(deg)
    use vars
    implicit none
    real(kind=pr), intent(in) :: deg
    deg2rad=deg*pi/180.d0
    return
  end function

  real(kind=pr) function rad2deg(deg)
    use vars
    implicit none
    real(kind=pr), intent(in) :: deg
    rad2deg=deg*180.d0/pi
    return
  end function

  function cross(a,b)
    use vars
    implicit none
    real(kind=pr),dimension(1:3),intent(in) :: a,b
    real(kind=pr),dimension(1:3) :: cross
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
  end function

  function norm2(a)
    use vars
    implicit none
    real(kind=pr),dimension(1:3),intent(in) :: a
    real(kind=pr) :: norm2
    norm2 = dsqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )
  end function

end module fsi_vars



! Variables for mhd simulations
module mhd_vars
  use vars
  implicit none

  ! Physical parameters
  real(kind=pr),save :: eta ! magnetic diffusivity
  real(kind=pr),save :: b0, bc ! Boundary condition parameters
  real(kind=pr),save :: ay ! x*x + ay*y*y -r1*r1 == 0 ?in boundary conditions

  ! Determine whether we save various fields
  integer,save :: iSaveMagneticField,iSaveCurrent
end module mhd_vars
