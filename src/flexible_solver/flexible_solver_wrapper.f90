!-------------------------------------------------------------------------------
!   FLEXIBLE WING SOLVER WRAPPER
! Input:
!   time: current time level (t^n)
!   dt: time step for new level (t^n+1 - t^n)
!   wings: array of wings, all at the old time level t^n
! Output:
!   wings: array of wings, all at the new time level t^n+1
!-------------------------------------------------------------------------------
subroutine FlexibleSolidSolverWrapper ( time, dt0, dt1, it, wings )

  use penalization ! mask array etc
  use vars

  implicit none
  real(kind=pr), intent (in) ::  dt0, dt1, time
  integer, intent (in) :: it
  real(kind=pr) :: t0
  type(flexible_wing), dimension(1:nWings), intent (inout) ::  wings
  integer :: i
  t0 = MPI_wtime()

  do i = 1, nWings
      !if (time>=T_release) then
        !-------------------------------------------
        ! the wings are released (now: active FSI), call IBES solvers
        !-------------------------------------------
        select case (TimeMethodFlexibleSolid)
        case ("BDF2","EI1")
            ! all implicit solvers are in one subroutine
            call Flexible_solid_time_step(time, dt0, dt1, it, wings(i))
            
        case ("prescribed_wing")
            ! this is not a solver, but for passive FSI with prescribed deformation:
            call prescribed_wing (time, wings(i))

            if (time_for_output(time, dt1, it, 1.d0, 99999999, 6.d0, 0.d0)) then

              if (root) then
              call external_forces_construction(time, it,wings(i))
              call SaveWingData( time, wings )
              endif
            endif

        case default
            call abort(723763,"FlexibleSolidSolver::invalid value of TimeMethodFlexibleSolid"//&
                 trim(adjustl(TimeMethodFlexibleSolid)))
        end select
      !else
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
  enddo

  ! cpu timing
  call toc("Flexible wing (FlexibleSolidSolverWrapper)", MPI_wtime() - t0)

end subroutine FlexibleSolidSolverWrapper


!-------------------------------------------------------------------------------
!   solid solver main entry point
!-------------------------------------------------------------------------------
subroutine OnlyFlexibleSolidSimulation()
  use vars
  use mpi
  implicit none
  type(flexible_wing), dimension(1:nWings) :: wings
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
  !call init_wings(wings)

  !call read_mesh_data(wings)



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
