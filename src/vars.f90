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

  ! Method variables set in the program file:
  character(len=strlen),save :: method ! mhd  or fsi
  character(len=strlen),save :: dry_run_without_fluid ! just save mask function
  integer,save :: nf  ! number of linear exponential fields (1 for HYD, 2 for MHD)
  integer,save :: nd  ! number of fields (3 for NS, 6 for MHD)
  integer,save :: neq ! number of fields in u-vector (3 for HYD, 6 for MHD, 4 for 
                      ! passive scalar (ux,uy,uz,theta)
  integer,save :: nrw ! number of real work arrays in work
  integer,save :: ncw ! number of complex work arrays in workc
  integer,save :: ng  ! number of ghostpoints (if used)

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
  
  ! p3dfft only parameters (move to appropraite .f90 file?)
  integer,save :: mpicommcart
  integer,dimension(2),save :: mpidims,mpicoords,mpicommslab
  ! only root rank has this true:
  logical, save :: root=.false.


  real(kind=pr),save :: pi ! 3.14....

  real(kind=pr),dimension(:),allocatable,save :: lin ! contains nu and eta

  ! Vabiables timing statistics.  Global to simplify syntax.
  real(kind=pr),save :: time_fft,time_ifft,time_vis,time_mask,time_nlk2
  real(kind=pr),save :: time_vor,time_curl,time_p,time_nlk,time_u
  real(kind=pr),save :: time_bckp,time_save,time_total,time_fluid,time_nlk_fft
  real(kind=pr),save :: time_sponge,time_insect_head,time_insect_body, time_scalar
  real(kind=pr),save :: time_insect_eye,time_insect_wings, time_insect_vel
  real(kind=pr),save :: time_solid, time_drag, time_surf, time_LAPACK
  real(kind=pr),save :: time_hdf5,time_integrals,time_rhs,time_nlk_scalar
  ! The mask array.  TODO: move out of shave_vars?
  real(kind=pr),dimension (:,:,:),allocatable,save :: mask ! mask function
  integer(kind=2),dimension (:,:,:),allocatable,save :: mask_color ! mask color function
  ! Velocity field inside the solid.  TODO: move out of shave_vars?
  real(kind=pr),allocatable,save :: us(:,:,:,:)  ! Velocity in solid

  ! Variables set via the parameters file
  real(kind=pr),save :: length 

  ! Domain size variables:
  integer,save :: nx,ny,nz
  real(kind=pr),save :: xl,yl,zl,dx,dy,dz,scalex,scaley,scalez

  ! FIXME: please document
  integer,save :: iDealias

  ! Parameters to set which files are saved and how often:
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask
  integer,save :: idobackup
  integer,save :: iSaveXMF !directly write *.XMF files (1) or not (0)
  real(kind=pr),save :: tintegral ! Time between output of integral quantities
  real(kind=pr),save :: tsave ! Time between outpout of entire fields.
  real(kind=pr),save :: tsave_first ! don't save before this time
  ! compute drag force every itdrag time steps and compute unst corrections if
  ! you've told to do so.
  integer,save :: itdrag, unst_corrections
  real(kind=pr),save :: truntime, truntimenext ! Number of hours bet
  real(kind=pr),save :: wtimemax ! Stop after a certain number of hours of wall.
  ! for periodically repeating flows, it may be better to always have only 
  ! one set of files on the disk
  character(len=strlen),save :: save_only_one_period
  real(kind=pr),save :: tsave_period ! then this is period time

  ! Time-stepping parameters
  real(kind=pr),save :: tmax
  real(kind=pr),save :: dt_fixed
  real(kind=pr),save :: dt_max=0.d0
  real(kind=pr),save :: cfl
  integer,save :: nt
  character(len=strlen),save :: iTimeMethodFluid

  ! Physical parameters:
  real(kind=pr),save :: nu

  ! Initial conditions:
  character(len=strlen),save :: inicond, file_ux,file_uy,file_uz
  character(len=strlen),save :: file_bx,file_by,file_bz
  real(kind=pr),save :: omega1 ! FIXME: what is omega1?


  ! Boundary conditions:
  character(len=strlen),save :: iMask
  integer,save :: iMoving,iPenalization
  real(kind=pr),save :: eps
  real(kind=pr),save :: r1,r2,r3 ! Parameters for boundary conditions
  character(len=strlen) :: iSmoothing ! how to smooth the mask
  real(kind=pr),save :: pseudoeps, pseudodt, pseudoerrmin, pseudoerrmax

  
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
      if (tmp>nx-1) tmp = tmp-nx
      if (nx==1) tmp=0
      per=tmp
      return
    end function per
    !---------------------------------------------------------------------------
    subroutine suicide1
      use mpi
      implicit none
      integer :: mpicode
      
      if (mpirank==0) write(*,*) "Killing run..."
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine suicide1
    !---------------------------------------------------------------------------
    subroutine suicide2(msg)
      use mpi
      implicit none
      integer :: mpicode
      character(len=*), intent(in) :: msg
      
      if (mpirank==0) write(*,*) "Killing run..."
      if (mpirank==0) write(*,*) msg
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine suicide2
    !---------------------------------------------------------------------------
    logical function is_nan( x )
      implicit none
      real(kind=pr)::x
      is_nan = .false.
      if (.not.(x.eq.x)) is_nan=.true.
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

  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle
  real(kind=pr),save :: Uxmean,Uymean,Uzmean, m_fluid
  character(len=strlen),save :: iMeanFlow_x,iMeanFlow_y,iMeanFlow_z
  integer,save :: iSaveSolidVelocity
  
  ! parameters for passive scalar advection
  integer, save :: use_passive_scalar
  integer, save :: n_scalars
  real(kind=pr),save :: kappa
  character(len=strlen),save :: inicond_scalar, stop_on_fail
  character(len=strlen),save :: source_term
  real(kind=pr), save :: eps_scalar
  logical, save :: compute_scalar
  real(kind=pr),save :: source_xmin,source_xmax,source_ymin,source_ymax,&
    source_zmin,source_zmax

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
     real(kind=pr),dimension(1:3) :: Force
     real(kind=pr),dimension(1:3) :: Force_unst
     real(kind=pr),dimension(1:3) :: Torque
     real(kind=pr),dimension(1:3) :: Torque_unst
  end type Integrals  
  
  !-----------------------------------------------------------------------------
  ! derived datatype for rigid solid dynamics solver
  type SolidDynType
    ! solid dynamics solver flag (0=off, 1=on)
    integer :: idynamics
    ! vector of unknowns at new time step
    real(kind=pr), dimension(1:4) :: var_new
    ! vector of unknowns at current time step
    real(kind=pr), dimension(1:4) :: var_this
    ! rhs at current time step
    real(kind=pr), dimension(1:4) :: rhs_this
    ! rhs at previous time step
    real(kind=pr), dimension(1:4) :: rhs_old
  end type SolidDynType
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
  type(SolidDynType), save :: SolidDyn
 
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
  
  function cross(a,b)
    use vars
    implicit none
    real(kind=pr),dimension(1:3),intent(in) :: a,b
    real(kind=pr),dimension(1:3) :: cross
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
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


