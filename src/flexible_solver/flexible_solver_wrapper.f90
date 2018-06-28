!-------------------------------------------------------------------------------
!   FLEXIBLE SOLID SOLVER WRAPPER
! Input:
!   time: current time level (t^n)
!   dt: time step for new level (t^n+1 - t^n)
!   wings: array of wings, all at the old time level t^n
! Output:
!   wings: array of wings, all at the new time level t^n+1
!-------------------------------------------------------------------------------
!subroutine FlexibleSolidSolverWrapper ( time, dt, wings )
!  use penalization ! mask array etc
!  implicit none
!  real(kind=pr), intent (in) ::  dt, time
!  real(kind=pr) :: t0
!  type(wing), dimension(1:nWings), intent (inout) ::    wings
!  integer :: i
!  t0 = MPI_wtime()

!  do i = 1, nWings
!      if (time>=T_release) then
        !-------------------------------------------
        ! the wings are released (now: active FSI), call IBES solvers
        !-------------------------------------------
!        select case (TimeMethodFlexibleSolid)
!            ! all implicit solvers are in one subroutine
!            call IBES_solver (time, dt, wings(i))
          !case ("prescribed")
            ! this is not a solver, but for passive FSI with prescribed deformation:
          !  call prescribed_beam (time, dt, wings(i))
!          case default
!            call abort(723763,"FlexibleSolidSolver::invalid value of TimeMethodFlexibleSolid"//&
!                 trim(adjustl(TimeMethodSolid)))
!        end select
!      else
        !-------------------------------------------
        ! the wings are not yet released, but their leading edges may move
        ! (passive FSI)
        !-------------------------------------------
!        call integrate_position (time+dt, wings(i))
!      endif

      !-------------------------------------------
      ! compute energies and stuff
      !-------------------------------------------
      !call MassSpringEnergies ( wings(i) )

      !-- check if everything seems okay, if not show beam and abort
      !call show_beam_on_error( wings(i) )
!  enddo
!  time_solid = time_solid + MPI_wtime() - t0
!end subroutine FlexibleSolidSolverWrapper


!-------------------------------------------------------------------------------
!   solid solver main entry point
!-------------------------------------------------------------------------------
subroutine OnlyFlexibleSolidSimulation()
  use vars
  use mpi
  use stl_file_reader
  use ini_files_parser_mpi
  implicit none
  type(wing), dimension(1:nWings) :: wings
  !real (kind=pr) :: time
  !integer :: it!,nsave

  ! in case you forgot to set it
  ! if (itdrag==0) itdrag=10
  ! if (dt_fixed<=1.0d-10) dt_fixed=1.d-3

  write (*,*) "*** information: starting OnlyFlexibleSolidSimulation"
  !time = 0.0
  !it = 0

  !call show_solid_model_information

  !-- initialization
  call init_wings(wings)

  call read_wing_mesh_data(wings)

  

  !--loop over time steps
  !do while ((time<=tmax))
    !-- external loads for testing purposes
    !if (debug_pressure==1) then
    !  wings(1)%pressure_new=0.1*dsin(time)
    !  wings(1)%tau_new=0.1*dsin(time)
    !endif

    !-- time stepping
    !call FlexibleSolidSolverWrapper( time, dt_fixed , wings )

  !  it   = it+1
  !  time = dble(it)*dt_fixed

  !  call SaveWingData( time, wings )
  !enddo
  !call SaveWingData( time, wings )

end subroutine
