
module flexible_model
  use vars ! for the precision statement
  ! we need the following line for presribed wings:
  use helpers
  use basic_operators
  implicit none

  !----------------------------------------------
  ! module global variables
  !----------------------------------------------
  integer,parameter :: nVeins = 5
  integer,parameter :: nVeins_BC = 2
  integer,parameter :: nMembranes = 1
  integer,save :: np, debug_pressure
  ! see "type solid" about nsmax 06 Aug 2014
  integer,parameter :: npmax = 1000
  integer,parameter :: nvmax = 2*npmax
  integer,parameter :: nmmax = 2*npmax
  ! TODO: move these into the solid model datastructure
  real(kind=pr),save :: EIy
  real(kind=pr),save :: EIz
  real(kind=pr),save :: grav
  real(kind=pr),save :: sigma
  real(kind=pr),save :: t_wing, L_vein, N_smooth
  real(kind=pr),save :: T_release
  real(kind=pr), parameter :: error_stop = 1.0e-6
  character(len=strlen),save :: imposed_motion_leadingedge, TimeMethodSolid
  character(len=strlen),save :: has_cylinder
  character(len=strlen),save :: interp
  character(len=strlen),save :: infinite, plate_shape

  ! this is a hack to avoid allocating/deallacting these arrays in every time
  ! step. turing fails with memory issues if done otherwise. 7 Aug 2014
  real(kind=pr),save,dimension(:,:,:,:), allocatable :: surfaces
  real(kind=pr),save,dimension(:,:,:), allocatable :: p_surface, p_surface_local
  real(kind=pr),save,dimension(:,:), allocatable :: active_points
  real(kind=pr),save,dimension(:),allocatable :: heights

  !----------------------------------------------
  ! Solid datatype
  !----------------------------------------------
  type Wing
    ! These arrays are statically allocated (they lie thus on
    ! the stack), since on turing we had problems with memory fragmentation. It
    ! is still possible to use np<npmax via the params file, but not np>npmax.

    ! coordinates, velocities and phase vector at time t^(n) and t^(n-1):
    real(kind=pr),dimension(1:npmax) :: x,y,z
    real(kind=pr),dimension(1:npmax) :: vx,vy,vz
    real(kind=pr),dimension(1:npmax) :: u_old, u_oldold
    integer,dimension(1:2*npmax,4) :: tri_elements
    real(kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) :: mask
    ! Matrix contains identification of vein points connected by bending springs
    ! and their initial and current bending angles on both directions y and z:
    real(kind=pr),dimension(1:nvmax,8,nVeins) :: Veins_bending
    ! Matrix contains identification of vein points connected by extension
    ! springs and their initial and current lengths:
    real(kind=pr),dimension(1:nvmax,5,nVeins) :: Veins_extension
    ! Veins with boundary conditions
    real(kind=pr),dimension(-1:nvmax,8,nVeins_BC) :: Veins_bending_BC
    real(kind=pr),dimension(-1:nvmax,5,nVeins_BC) :: Veins_extension_BC
    ! Membrane:
    real(kind=pr),dimension(1:nmmax,5,nMembranes) :: Membranes_extension
    real(kind=pr),dimension(1:nmmax,5) :: Membrane_edge
    ! external forces:
    real(kind=pr),dimension(0:nsmax) :: pressure_old, pressure_new
    ! material properties:
    real(kind=pr),dimension(1:nVeins) :: EIy, EIz, kby0, kbz0
    real(kind=pr),dimension(-1:nVeins_BC) :: EIy_BC, EIz_BC, kby0_BC, kbz0_BC
    real(kind=pr),dimension(1:nvmax,1:nVeins) :: ke_v, kb_v
    real(kind=pr),dimension(-1:nvmax,1:nVeins_BC) :: ke_vBC, kb_vBC
    real(kind=pr),dimension(1:nvmax,1:nMembranes) :: kb_m
    real(kind=pr),dimension(1:nvmax) :: kb_me

    ! grid and width in rigid direction:
    real(kind=pr),dimension(1:nVeins+nBCs) :: L_veins
    ! values of theta and theta_dot at times t^(n) and t^(n-1)
    !real(kind=pr),dimension(0:nsmax) :: theta_old, theta_oldold
    !real(kind=pr),dimension(0:nsmax) :: theta_dot_old, theta_dot_oldold

    ! real(kind=pr),dimension(0:nsmax,1:6) :: wing_oldold
    !real(kind=pr),dimension(1:2) :: Force, Force_unst, Force_press, Inertial_Force
    !real(kind=pr) :: E_kinetic, E_elastic
    real(kind=pr) :: x0, y0, z0, Anglewing_y, Anglewing_z, phase
    ! we need the previous time step for the BDF solver
    real(kind=pr) :: dt_old
    ! these replace the save variables in the unst correction computation:
    !real(kind=pr) :: drag_unst_new, drag_unst_old, lift_unst_new, lift_unst_old
    logical :: StartupStep!, UnsteadyCorrectionsReady
  end type wing


 contains


 !-----------------------------------------------------------------
 include "init_wing.f90"
 include "read_wing_mesh_data.f90"
! include "save_wing.f90"
! include "wingForces.f90"
! include "prescribed_wing.f90"
 include "flexible_solver_wrapper.f90"



!-------------------------------------------------------------------------------
!   SOLID SOLVER INITIALIZATION
!-------------------------------------------------------------------------------
subroutine InitializeFlexibleSolidSolver( wings )
  implicit none
  integer :: i
  type(wing), dimension(1:nWings), intent (inout) :: wings

  ! marks all wings to be in the very first time step
  ! the solver then uses CN2 instead of BDF2, because the old old time level
  ! t_n-1 is not available
  do i=1,nWings
    wings(i)%StartupStep = .true.
  enddo
end subroutine InitializeFlexibleSolidSolver


!-------------------------------------------------------------------------------
!   energies and stuff for wings
!-------------------------------------------------------------------------------
!subroutine MassSpringEnergies( wing )
!  implicit none
!  type(wing), intent (inout) :: wing

!  wing%E_kinetic = 0.5d0*ds*sum( wing%mu(0:ns-1) * (wing%vx(0:ns-1)**2 + wing%vy(0:ns-1)**2) )
!  wing%E_pot     = grav*ds *sum( wing%mu(0:ns-1) * (wing%y(0:ns-1)-wing%y0) )
!  wing%E_elastic = 0.5d0*ds*sum( wing%zeta(0:ns-1) * theta_s(0:ns-1)**2 )

!  wing%Inertial_Force(1) = ds* sum( wing%mu(0:ns-1) * wing%ax(0:ns-1) )
!  wing%Inertial_Force(2) = ds* sum( wing%mu(0:ns-1) * wing%ay(0:ns-1) )

!end subroutine MassSpringEnergies


end module flexible_model
