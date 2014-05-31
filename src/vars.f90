! Variables for pseudospectral simnulations
module vars
  use mpi
  implicit none
  
  interface suicide
    module procedure suicide1,suicide2
  end interface

  character(len=1),save:: tab ! Fortran lacks a native tab, so we set one up.
  ! Used in params.f90
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file
  integer,parameter :: strlen=80   ! standard string length

  ! Precision of doubles
  integer,parameter :: pr = 8 

  ! Method variables set in the program file:
  character(len=strlen),save:: method ! mhd  or fsi
  character(len=strlen), save :: dry_run_without_fluid ! just save mask function
  integer,save :: nf ! number of linear exponential fields (1 for HYD, 2 for MHD)
  integer,save :: nd ! number of fields (3 for NS, 6 for MHD)
  integer,save :: ng ! number of ghostpoints (if used)

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
  

  real(kind=pr),save :: pi ! 3.14....

  real(kind=pr),dimension(:),allocatable,save :: lin ! contains nu and eta

  ! Vabiables timing statistics.  Global to simplify syntax.
  real(kind=pr),save :: time_fft,time_ifft,time_vis,time_mask,time_fft2
  real(kind=pr),save :: time_vor,time_curl,time_p,time_nlk,time_u, time_ifft2
  real(kind=pr),save :: time_bckp,time_save,time_total,time_fluid,time_nlk_fft
  real(kind=pr),save :: time_sponge, time_solid, time_drag, time_surf, time_LAPACK
  
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

  ! Am I a fluid or a solid CPU?
  logical, save :: fluid_cpu, solid_cpu
  ! number of CPU, of CPU for fluid and CPU for solid
  integer, save :: ncpu, ncpu_fluid, ncpu_solid
  ! communicators:
  integer, save :: MPI_COMM_FLUID, MPI_COMM_SOLID
  ! only root rank has this true:
  logical, save :: root=.false.
  
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
    !-----------------------------------------------
        subroutine checknan_real( field, msg)
      implicit none
      integer::ix,iy,iz
      real(kind=pr),intent(in):: field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
      character(len=*), intent(in) :: msg
      integer :: foundnan, foundnans,mpicode
      foundnan = 0
      do ix=ra(1),rb(1)
      do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)  
      if (.not.(field(ix,iy,iz).eq.field(ix,iy,iz))) foundnan = 1
      enddo
      enddo
      enddo  
      
      call MPI_ALLREDUCE (foundnan,foundnans,1,&
                    MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)  
                    
      if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans      
    end subroutine checknan_real
    
    subroutine checknan_cmplx( field, msg )
      implicit none
      integer::ix,iy,iz
      complex(kind=pr),intent(in)::field(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
      character(len=*), intent(in) :: msg
      integer :: foundnan, foundnans,mpicode
      foundnan = 0
      do iz=ca(1),cb(1)
      do iy=ca(2),cb(2)
      do ix=ca(3),cb(3)
      if (.not.(aimag(field(iz,iy,ix)).eq.aimag(field(iz,iy,ix)))) foundnan = 1
      if (.not.(real(field(iz,iy,ix)).eq.real(field(iz,iy,ix)))) foundnan = 1
      enddo
      enddo
      enddo  
      
      call MPI_ALLREDUCE (foundnan,foundnans,1,&
                    MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpicode)  
                    
      if (root.and.foundnans>0) write(*,'("NaN in ",A," sum=",i5)') msg, foundnans
    end subroutine checknan_cmplx
end module vars


! Variables for fsi simulations
module fsi_vars
  use vars
  implicit none
  
  ! the sponge term is temporarily set in a global variable, but this will
  ! be changed in future versions
  complex(kind=pr), allocatable, save:: sponge(:,:,:,:)
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
  
  

  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle
  real(kind=pr),save :: Uxmean,Uymean,Uzmean
  integer,save :: iMeanFlow
  integer,save :: iSaveSolidVelocity

  character(len=strlen),save :: use_solid_model
  !-----------------------------------------------------------------------------
  ! The derived integral quantities for fluid-structure interactions.
  type Integrals
     real(kind=pr) :: time
     real(kind=pr) :: EKin
     real(kind=pr) :: Dissip
     real(kind=pr) :: Divergence
     real(kind=pr) :: Volume
     real(kind=pr),dimension(1:3) :: Force
     real(kind=pr),dimension(1:3) :: Force_unst
     real(kind=pr),dimension(1:3) :: Torque
     real(kind=pr),dimension(1:3) :: Torque_unst
  end type Integrals  
  !-----------------------------------------------------------------------------
  ! derived datatype for insect parameters (for readability)
  type InsectParams ! documentaion see insect.f90
    character(len=strlen) :: WingShape, BodyType, HasHead, HasEye, BodyMotion
    character(len=strlen) :: FlappingMotion_right, FlappingMotion_left
    character(len=strlen) :: KineFromFile ! kinematics from file (Dmitry, 14 Nov 2013)
    ! parameters for body:
    real(kind=pr) :: L_body, b_body, R_head, R_eye
    ! parameters for wing shape:
    real(kind=pr) :: b_top, b_bot, L_chord, L_span, WingThickness
    ! this is a safety distance for smoothing:
    real(kind=pr) :: safety, smooth
    ! vectors desribing the positoions of jerry's key elements
    ! in the body coordinate system
    real(kind=pr), dimension(1:3) :: x_head,x_eye_r,x_eye_l,x_pivot_l,x_pivot_r
    ! parameter for hovering:
    real(kind=pr) :: distance_from_sponge
    ! Wings and body forces
    type(Integrals), dimension(1:2) :: PartIntegrals
    ! Takeoff parameters
    real(kind=pr) :: x_takeoff, z_takeoff, mass_solid, gravity 
    ! Legs model parameters
    integer :: ilegs
    real(kind=pr) :: anglegsend, kzlegsmax, dzlegsmax, t0legs, tlinlegs
  end type InsectParams  
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
  type(InsectParams), save :: Insect
  type(SolidDynType), save :: SolidDyn
 
  contains
  !-----------------------------------------------------------------------------
  
  ! this function simplifies my life with the insects
  function deg2rad(deg)
    use vars
    implicit none
    real(kind=pr), intent(in) :: deg
    real(kind=pr) :: deg2rad
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


