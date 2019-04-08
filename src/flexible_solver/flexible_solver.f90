
module flexible_model
  use vars ! for the precision statement
  ! we need the following line for presribed wings:
  use module_helpers
  use basic_operators
  use interpolation
  use module_ini_files_parser_mpi
  implicit none

  !----------------------------------------------
  ! module global variables
  !----------------------------------------------
  integer,parameter :: nWings = 1
  integer,parameter :: nVeins = 38
  integer,parameter :: nVeins_BC = 5
  integer,parameter :: nMembranes = 1
  integer,parameter :: nMembrane_edges = 10
  ! see "type solid" about nsmax 06 Aug 2014
  integer,parameter :: npmax = 500
  integer,parameter :: nvmax = 2*npmax
  integer,parameter :: nmmax = 3*npmax
  ! TODO: move these into the solid model datastructure

  real(kind=pr),dimension(1:3),save :: grav
  real(kind=pr), parameter :: error_stop = 1.0e-6
  logical :: ActuallyBDF2
  character(len=strlen),save :: TimeMethodFlexibleSolid
  character(len=strlen),save :: use_flexible_wing_model
  character(len=strlen),save :: load_mass_from_file
  character(len=strlen),save :: activate_press_force,activate_noninertial_force,wing_interp

  ! this is a hack to avoid allocating/deallacting these arrays in every time
  ! step. turing fails with memory issues if done otherwise. 7 Aug 2014
  !real(kind=pr),save,dimension(:,:,:,:), allocatable :: surfaces
  !real(kind=pr),save,dimension(:,:,:), allocatable :: p_surface, p_surface_local
  !real(kind=pr),save,dimension(:,:), allocatable :: active_points
  !real(kind=pr),save,dimension(:),allocatable :: heights

  !----------------------------------------------
  ! Wing datatype
  !----------------------------------------------
  type Flexible_wing
    ! These arrays are statically allocated (they lie thus on
    ! the stack), since on turing we had problems with memory fragmentation. It
    ! is still possible to use np<npmax via the params file, but not np>npmax.

    ! coordinates, velocities and phase vector at time t^(n) and t^(n-1):
    real(kind=pr),dimension(1:npmax) :: x,y,z
    real(kind=pr),dimension(1:npmax) :: vx,vy,vz
    real(kind=pr),dimension(1:6*npmax) :: du, u_new, u_old, u_oldold
    real(kind=pr),dimension(-1:0,1:nVeins_BC) :: x_BC,y_BC,z_BC,x0_BC,y0_BC,z0_BC
    integer,dimension(1:nmmax,4) :: tri_elements
    real(kind=pr),dimension(1:nmmax) :: tri_element_areas
    real(kind=pr),dimension(1:nmmax,4) :: tri_element_normals
    integer :: np, ntri

    ! Veins :
      real(kind=pr),dimension(1:nvmax,2,nVeins) :: Veins
      ! Matrix contains identification of vein points connected by bending springs
      ! and their initial and current bending angles on both directions y and z:
      real(kind=pr),dimension(1:nvmax,8,nVeins) :: Veins_bending
      ! Matrix contains identification of vein points connected by extension
      ! springs and their initial and current lengths:
      real(kind=pr),dimension(1:nvmax,5,nVeins) :: Veins_extension
      ! Veins with boundary conditions
      real(kind=pr),dimension(1:nvmax,2,nVeins_BC) :: Veins_BC
      real(kind=pr),dimension(-1:nvmax,8,nVeins_BC) :: Veins_bending_BC
      real(kind=pr),dimension(0:nvmax,5,nVeins_BC) :: Veins_extension_BC
      ! Vein connectors
      real(kind=pr),dimension(1:nvmax,8) :: Vein_connectors

    ! Membrane:
      real(kind=pr),dimension(1:nmmax,2,nMembranes) :: Membranes
      real(kind=pr),dimension(1:nmmax,5,nMembranes) :: Membranes_extension
      real(kind=pr),dimension(1:nmmax,5,nMembrane_edges) :: Membrane_edge
    ! Internal and external forces:
    real(kind=pr),dimension(1:3*npmax) :: Fint, Fext !, Fint_old, Fext_old
    real(kind=pr),dimension(1:nmmax,3) :: press_upside, press_downside
    ! Internal force derivative matrix:
    real(kind=pr),dimension(1:3*npmax,1:3*npmax) :: FJ, FJ_ext
    ! material properties:
    real(kind=pr),dimension(1:nVeins) :: EIy, EIz, kby0, kbz0, ke0_v
    real(kind=pr),dimension(1:nVeins_BC) :: EIy_BC, EIz_BC, kby0_BC, kbz0_BC, ke0_vBC
    real(kind=pr),dimension(1:nVeins) :: rho_v
    real(kind=pr),dimension(1:nVeins_BC) :: rho_vBC
    real(kind=pr),dimension(1:nMembranes) :: rho_m, ke0_m
    real(kind=pr) :: rho_me
    real(kind=pr) :: c0

    real(kind=pr),dimension(1:nvmax,1:nVeins) :: ke_v, kby, kbz
    real(kind=pr),dimension(1:nvmax,1:nVeins) :: m_v
    real(kind=pr),dimension(1:nvmax) :: kby_c, kbz_c
    real(kind=pr),dimension(-1:nvmax,1:nVeins_BC) :: ke_vBC, kby_BC, kbz_BC
    real(kind=pr),dimension(-1:nvmax,1:nVeins_BC) :: m_vBC
    real(kind=pr),dimension(1:nmmax,1:nMembranes) :: ke_m
    real(kind=pr),dimension(1:nmmax,1:nMembranes) :: m_m
    real(kind=pr),dimension(1:nmmax,1:nMembrane_edges) :: ke_me
    real(kind=pr),dimension(1:nmmax) :: m_me
    real(kind=pr),dimension(1:npmax) :: m, c

    ! Position, velocity and acceleration (t: translation, r: rotation) of the
    ! local coordinate system
    real(kind=pr) :: x0, y0, z0
    real(kind=pr) :: WingAngle_x, WingAngle_y, WingAngle_z
    real(kind=pr),dimension(1:3,1:3) :: RotationMat_x,RotationMat_y,RotationMat_z
    integer :: NumCorrection, ControlPoint
    real(kind=pr),dimension(1:3) :: vt0, at0
    real(kind=pr),dimension(1:3) :: vr0, ar0


    ! grid and width in rigid direction:
    real(kind=pr) :: t_wing, wing_smoothing

    ! real(kind=pr),dimension(0:nsmax,1:6) :: wing_oldold

    !real(kind=pr),dimension(1:2) :: Force, Force_unst, Force_press, Inertial_Force
    !real(kind=pr) :: E_kinetic, E_elastic




    ! we need the previous time step for the BDF solver
    real(kind=pr) :: dt_old

    ! these replace the save variables in the unst correction computation:
    !real(kind=pr) :: drag_unst_new, drag_unst_old, lift_unst_new, lift_unst_old
    logical :: StartupStep!, UnsteadyCorrectionsReady
    !character(len=strlen) :: TimeMethodFlexibleSolid
    character(len=strlen) :: Motion
    real(kind=pr),dimension(1:3*npmax) :: RHS_a, RHS_b
    real(kind=pr) :: T_release, tau

  end type flexible_wing


 contains


 !-----------------------------------------------------------------
 include "init_wing.f90"
 include "save_wing.f90"
 include "internal_force.f90"
 include "external_force.f90"
 include "external_force_derivative.f90"
 include "internal_force_derivative.f90"
 include "flexible_wing_motions.f90"
 include "flexible_solver_wrapper.f90"
 include "flexible_solid_time_stepper.f90"
 include "flexible_tri_mask.f90"
 include "prescribed_wing.f90"
 include "supplementary_calc.f90"



!-------------------------------------------------------------------------------
!   energies and stuff for wings
!-------------------------------------------------------------------------------
!subroutine MassSpringEnergies( wing )
!  implicit none
!  type(flexible_wing), intent (inout) :: wing

!  wing%E_kinetic = 0.5d0*ds*sum( wing%mu(0:ns-1) * (wing%vx(0:ns-1)**2 + wing%vy(0:ns-1)**2) )
!  wing%E_pot     = grav*ds *sum( wing%mu(0:ns-1) * (wing%y(0:ns-1)-wing%y0) )
!  wing%E_elastic = 0.5d0*ds*sum( wing%zeta(0:ns-1) * theta_s(0:ns-1)**2 )

!  wing%Inertial_Force(1) = ds* sum( wing%mu(0:ns-1) * wing%ax(0:ns-1) )
!  wing%Inertial_Force(2) = ds* sum( wing%mu(0:ns-1) * wing%ay(0:ns-1) )

!end subroutine MassSpringEnergies


end module flexible_model
