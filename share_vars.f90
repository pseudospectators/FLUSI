! Variables for pseudospectral simnulations
module vars
  use mpi_header
  implicit none

  character*1,save:: tab ! Fortran lacks a native tab, so we set one up.

  ! Precision of doubles
  integer,parameter :: pr = 8 

  ! Method variables set in the program file:
  character(len=3),save:: method ! mhd  or fsi
  integer,save :: nf ! number of fields (1 for NS, 2 for MHD)
  integer,save :: nd ! number of fields (3 for NS, 6 for MHD)

  ! MPI and p3dfft variables and parameters
  integer,save :: mpisize,mpirank,mpicommcart
  integer,parameter :: mpiinteger=MPI_INTEGER
  integer,parameter :: mpireal=MPI_DOUBLE_PRECISION
  integer,parameter :: mpicomplex=MPI_DOUBLE_COMPLEX
  integer,dimension(2),save :: mpidims,mpicoords,mpicommslab
  integer,dimension (:,:),allocatable,save :: ra_table,rb_table
  ! Local array bounds
  integer,dimension (1:3),save :: ra,rb,rs,ca,cb,cs

  ! Used in params.f90
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file

  real(kind=pr),save :: pi ! 3.14....

  real (kind=pr),dimension(:),allocatable,save :: lin ! contains nu and eta

  ! Vabiables timing statistics.  Global to simplify syntax.
  real (kind=pr),save :: time_fft,time_ifft,time_vis,time_mask,time_fft2
  real (kind=pr),save :: time_vor,time_curl,time_p,time_nlk,time_u, time_ifft2
  real (kind=pr),save :: time_bckp,time_save,time_total,time_fluid,time_nlk_fft
  real (kind=pr),save :: time_sponge
  
  ! The mask array.  TODO: move out of shave_vars?
  real (kind=pr),dimension (:,:,:),allocatable,save :: mask ! mask function
  ! Velocity field inside the solid.  TODO: move out of shave_vars?
  real (kind=pr),allocatable,save :: us(:,:,:,:)  ! Velocity in solid


  ! Variables set via the parameters file
  real(kind=pr),save :: length 
  ! Q: what is length? 
  ! A: a generic lengthscale, for example circle radius or plate spanwise length
  
  ! Domain size variables:
  integer,save :: nx,ny,nz
  real(kind=pr),save :: xl,yl,zl,dx,dy,dz,scalex,scaley,scalez

  ! FIXME: please document
  integer,save :: iDealias

  ! Parameters to set which files are saved and how often:
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask
  integer,save :: iDoBackup
  integer,save :: iSaveXMF !directly write *.XMF files (1) or not (0)
  real(kind=pr),save :: tintegral ! Time between output of integral quantities
  real(kind=pr),save :: tsave ! Time between outpout of entire fields.
  integer,save :: itdrag

  ! Time-stepping parameters
  real(kind=pr),save :: tmax
  real(kind=pr),save :: cfl
  integer,save :: nt
  character(len=80),save :: iTimeMethodFluid

  ! Physical parameters:
  real(kind=pr),save :: nu

  ! Initial conditions:
  character(len=80),save :: inicond, file_ux,file_uy,file_uz
  character(len=80),save :: file_bx,file_by,file_bz
  real(kind=pr),save :: omega1 ! FIXME: what is omega1?

  ! Boundary conditions:
  character(len=80),save :: iMask
  integer,save :: iMoving,iPenalization
  real(kind=pr),save :: dt_fixed
  real(kind=pr),save :: eps
  real(kind=pr),save :: r1,r2,r3 ! Parameters for boundary conditions
  character(len=80) :: iSmoothing ! how to smooth the mask
end module vars


! Variables for fsi simulations
module fsi_vars
  use vars
  implicit none
  
  ! the sponge term is temporarily set in a global variable, but this will
  ! be changed in future versions
  complex(kind=pr), allocatable, save:: sponge(:,:,:,:)
  ! switch for vorticity sponge:
  character(len=80), save :: iVorticitySponge
  real(kind=pr), save :: eps_sponge
  integer, save :: sponge_thickness
  
  ! for periodically repeating flows, it is better to always have only one set
  ! of files on the disk
  character(len=80) :: save_only_one_period
  

  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle
  real(kind=pr),save :: Uxmean,Uymean,Uzmean
  integer,save :: iMeanFlow
  integer,save :: iSaveSolidVelocity

  ! The derived integral quantities for fluid-structure interactions.
  type Integrals
     real(kind=pr) :: time
     real(kind=pr) :: EKin
     real(kind=pr) :: Dissip
     real(kind=pr) :: Divergence
     real(kind=pr) :: Volume
     real(kind=pr),dimension(1:3) :: Force
     real(kind=pr),dimension(1:3) :: Torque
  end type Integrals
  
  ! derived datatype for insect parameters (for readability)
  type InsectParams ! documentaion see insect.f90
    character(len=80) :: WingShape, BodyType, HasHead, HasEye, BodyMotion
    character(len=80) :: FlappingMotion_right, FlappingMotion_left
    ! parameters for body:
    real(kind=pr) :: L_body, b_body, R_head, R_eye
    ! parameters for wing shape:
    real(kind=pr) :: b_top, b_bot, L_chord, L_span, WingThickness
    ! this is a safety distance for smoothing:
    real(kind=pr) :: safety, smooth
    ! vectors desribing the positoions of jerry's key elements
    ! in the body coordinate system
    real(kind=pr), dimension(1:3) :: x_head, x_eye_r, x_eye_l, &
                                     x_pivot_l, x_pivot_r
  end type InsectParams

  type(Integrals),save :: GlobalIntegrals
  type(InsectParams), save :: Insect
  
  contains
  
  ! this function simplifies my life with the insects
  real(kind=pr) function deg2rad(deg)
    use vars
    implicit none
    real(kind=pr), intent(in) :: deg
    deg2rad=deg*pi/180.d0
    return
  end function
end module fsi_vars


! Variables for mhd simulations
module mhd_vars
  use vars
  implicit none

  ! Physical parameters
  real(kind=pr),save :: eta ! magnetic diffusivity
  real(kind=pr),save :: b0, bc ! Boundary condition parameters

  ! Determine whether we save various fields
  integer,save :: iSaveMagneticField,iSaveCurrent
end module mhd_vars


