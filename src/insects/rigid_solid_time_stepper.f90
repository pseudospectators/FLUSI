! Rigid solid time stepping routines
! AB2 method with Euler startup
subroutine rigid_solid_time_step(time,dt0,dt1,it)
  use mpi
  use fsi_vars
  implicit none

  real (kind=pr),intent (in) :: time,dt1,dt0
  integer,intent (in) :: it
  real (kind=pr) :: b10,b11

  ! select scheme
  if(it == 0) then
     ! update vector of solid variables
     SolidDyn%var_this(:) = SolidDyn%var_new(:)
     ! compute rhs at this time step
     call rigid_solid_rhs(time,it)
     ! Euler step
     SolidDyn%var_new(:) = SolidDyn%var_this(:) + dt1*SolidDyn%rhs_this(:)
  else
     ! update vectors
     SolidDyn%rhs_old(:) = SolidDyn%rhs_this(:)
     SolidDyn%var_this(:) = SolidDyn%var_new(:)
     ! compute rhs at this time step
     call rigid_solid_rhs(time,it)
     ! adaptive AB2 coefficients 
     b10=dt1/dt0*(0.5*dt1 + dt0)
     b11=-0.5*dt1*dt1/dt0
     ! adaptive AB2 step
     SolidDyn%var_new(:) = SolidDyn%var_this(:) + &
             b10*SolidDyn%rhs_this(:) + b11*SolidDyn%rhs_old(:)
  endif

end subroutine rigid_solid_time_step

! RHS of the ODEs that describe the rigid solid dynamics
! This is a wrapper. Actual implementation is problem-specific
subroutine rigid_solid_rhs(time,it)
  use mpi
  use fsi_vars
  implicit none

  real (kind=pr),intent (in) :: time
  integer,intent (in) :: it

  select case (iMask)
  case ("Insect")
    call dynamics_insect (time,it) ! see insects.f90
  case default
    if (mpirank == 0) then
        write (*,*) &
        "rigidsolidtimestepper.f90::rigid_solid_rhs case not defined"
        call abort()
    endif
  end select
end subroutine rigid_solid_rhs


! Initialize the rigid solid ODE solve
! This is a wrapper. Actual implementation is problem-specific.
! determines if the solid dynamics solver should be activated
! returns idynamics = 0 (do not activate) or 1 (activate)
subroutine rigid_solid_init(idynamics)
  use mpi
  use fsi_vars
  implicit none

  integer, intent(out) :: idynamics

  ! Solid dynamic solver not active by default
  idynamics = 0
  select case (iMask)
  case ("Insect")
    call dynamics_insect_init(idynamics)
  end select
end subroutine rigid_solid_init


