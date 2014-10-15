! Variables for pseudospectral simnulations
module vars
  use mpi
  implicit none
  
  
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file
  integer,parameter :: strlen=80   ! standard string length
  ! Precision of doubles
  integer,parameter :: pr = 8 
  
  !-----------------------------------------------------------------------------
  ! Type declarations
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
    real(kind=pr) :: dt_new
    real(kind=pr) :: dt_old
    integer :: it, it_start 
    integer :: n0
    integer :: n1
  end type timetype
  
  
  
  !-----------------------------------------------------------------------------
  ! Global parameters and variables
  !-----------------------------------------------------------------------------

  ! Method variables set in the program file:
  character(len=strlen),save :: method ! mhd  or fsi
  character(len=strlen),save :: dry_run_without_fluid ! just save mask function
  
  integer,save :: neq ! number of equations
  integer,save :: nrw ! number of real work arrays in work
  integer,save :: ng  ! number of ghostpoints (if used)
  integer,save :: nrhs ! number of registers for right hand side vectors

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

  ! Vabiables timing statistics.  Global to simplify syntax.
  real(kind=pr),save :: time_mask
  real(kind=pr),save :: time_vor, time_p
  real(kind=pr),save :: time_bckp, time_save, time_total, time_fluid
  real(kind=pr),save :: time_insect_body, time_scalar
  real(kind=pr),save :: time_insect_wings, time_insect_vel
  real(kind=pr),save :: time_solid, time_drag, time_surf, time_LAPACK
  real(kind=pr),save :: time_hdf5, time_integrals, time_rhs, time_nlk_scalar

  ! Variables set via the parameters file
  real(kind=pr),save :: length 

  ! Domain size variables:
  integer,save :: nx,ny,nz
  real(kind=pr),save :: xl,yl,zl,dx,dy,dz,scalex,scaley,scalez

  ! Parameters to set which files are saved and how often:
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask
  integer,save :: iSaveSolidVelocity
  integer,save :: idobackup
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
  character(len=strlen),save :: save_only_one_period
  real(kind=pr),save :: tsave_period ! then this is period time

  ! Time-stepping parameters
  real(kind=pr),save :: tmax
  real(kind=pr),save :: dt_fixed
  real(kind=pr),save :: dt_max=0.d0
  real(kind=pr),save :: cfl
  integer,save :: nt
  character(len=strlen),save :: iTimeMethodFluid

  ! viscosity (inverse of Reynolds number:)
  real(kind=pr),save :: nu
  
  ! pseudo speed of sound for the artificial compressibility method
  real(kind=pr), save :: c_0, gamma_p

  ! Initial conditions:
  character(len=strlen),save :: inicond,file_ux,file_uy,file_uz,file_p


  ! Boundary conditions:
  character(len=strlen),save :: iMask
  integer,save :: iMoving,iPenalization
  real(kind=pr),save :: eps
  real(kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle
  ! cavity mask:
  character(len=strlen), save :: iCavity, iChannel
  integer, save :: cavity_size
  ! wall thickness
  real(kind=pr),save :: thick_wall
  ! wall position (solid from pos_wall to pos_wall+thick_wall)
  real(kind=pr),save :: pos_wall 
  

  ! save forces and use unsteady corrections?
  integer, save :: compute_forces  
  
  ! mean flow control
  real(kind=pr),save :: Uxmean,Uymean,Uzmean, m_fluid
  character(len=strlen),save :: iMeanFlow_x,iMeanFlow_y,iMeanFlow_z
  ! mean flow startup conditioner (if "dynamic" and mean flow at t=0 is not zero
  ! the forces are singular at the beginning. use the startup conditioner to 
  ! avoid large accelerations in mean flow at the beginning)
  character(len=strlen),save :: iMeanFlowStartupConditioner
  real(kind=pr) :: tau_meanflow, T_release_meanflow
  
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

  type(Integrals),save :: GlobalIntegrals
  type(SolidDynType), save :: SolidDyn
  
  
  !-----------------------------------------------------------------------------
  ! Small helper functions
  !-----------------------------------------------------------------------------
  contains 

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
      implicit none
      integer :: mpicode
      
      if (mpirank==0) write(*,*) "Killing run..."
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine suicide1
    !---------------------------------------------------------------------------
    subroutine suicide2(msg)
      implicit none
      integer :: mpicode
      character(len=*), intent(in) :: msg
      
      if (mpirank==0) write(*,*) "Killing run..."
      if (mpirank==0) write(*,*) msg
      call MPI_abort(MPI_COMM_WORLD,666,mpicode)
    end subroutine suicide2
    !---------------------------------------------------------------------------
    ! wrapper for NaN checking (this may be compiler dependent)
    logical function is_nan( x )
      implicit none
      real(kind=pr)::x
      is_nan = .false.
      if (.not.(x.eq.x)) is_nan=.true.
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
    !---------------------------------------------------------------------------
    ! this function simplifies my life with the insects
    real(kind=pr) function deg2rad(deg)
      implicit none
      real(kind=pr), intent(in) :: deg
      deg2rad=deg*pi/180.d0
      return
    end function
    !---------------------------------------------------------------------------
    function cross(a,b)
      implicit none
      real(kind=pr),dimension(1:3),intent(in) :: a,b
      real(kind=pr),dimension(1:3) :: cross
      cross(1) = a(2)*b(3)-a(3)*b(2)
      cross(2) = a(3)*b(1)-a(1)*b(3)
      cross(3) = a(1)*b(2)-a(2)*b(1)
    end function
    !---------------------------------------------------------------------------
end module vars