!-------------------------------------------------------------------------------
!   FLEXIBLE WING SOLVER WRAPPER
! Input:
!   time: current time level (t^n)
!   dt: time step for new level (t^n+1 - t^n)
!   wings: array of wings, all at the old time level t^n
! Output:
!   wings: array of wings, all at the new time level t^n+1
!-------------------------------------------------------------------------------
subroutine FlexibleSolidSolverWrapper ( time, dt0, dt1, it, wings, Insect )

  use penalization ! mask array etc
  use vars
  use module_insects

  implicit none
  real(kind=pr), intent (in) ::  dt0, dt1, time
  integer, intent (in) :: it
  real(kind=pr) :: t0
  type(diptera), intent(inout) :: Insect
  type(flexible_wing), dimension(1:nWings), intent (inout) ::  wings
  integer :: i, mpicode
  logical :: run_in_parallel
  t0 = MPI_wtime()

  run_in_parallel = .true.

  do i = 1, nWings
    if (run_in_parallel) then
      if (mpirank==i-1) then
        !if (time>=T_release) then
          !-------------------------------------------
          ! the wings are released (now: active FSI), call IBES solvers
          !-------------------------------------------
          select case (TimeMethodFlexibleSolid)
          case ("BDF2","EI1")
              ! all implicit solvers are in one subroutine
              call Flexible_solid_time_step(time, dt0, dt1, it, wings(i), Insect)

          case ("prescribed_wing")
              ! this is not a solver, but for passive FSI with prescribed deformation:
              call prescribed_wing (time, wings(i), Insect)

              if (time_for_output(time, dt1, it, 1.d0, 99999999, 6.d0, 0.d0)) then

                if (root) then
                call external_forces_construction(time, it,wings(i),Insect)
                call SaveWingData( time, it, dt1, wings )
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
        endif
  else
      select case (TimeMethodFlexibleSolid)
      case ("BDF2","EI1")
          ! all implicit solvers are in one subroutine
          call Flexible_solid_time_step(time, dt0, dt1, it, wings(i), Insect)

      case ("prescribed_wing")
          ! this is not a solver, but for passive FSI with prescribed deformation:
          call prescribed_wing (time, wings(i), Insect)

          if (time_for_output(time, dt1, it, 1.d0, 99999999, 6.d0, 0.d0)) then

            if (root) then
            call external_forces_construction(time, it,wings(i),Insect)
            call SaveWingData( time, it, dt1, wings )
            endif
          endif

      case default
          call abort(723763,"FlexibleSolidSolver::invalid value of TimeMethodFlexibleSolid"//&
               trim(adjustl(TimeMethodFlexibleSolid)))
      end select
  endif
  enddo



  ! cpu timing
  call toc("Flexible wing (FlexibleSolidSolverWrapper)", MPI_wtime() - t0)

  run_in_parallel = .true.

  if (run_in_parallel) then
  t0 = MPI_wtime()


  do i =1, nwings
    !call MPI_Bcast( wings(i), nWings, Flexible_wing, i, MPI_COMM_WORLD, mpicode )
    !call MPI_Bcast( wings(i), nWings, Flexible_wing, i, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%y, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%z, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vx, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vy, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vz, npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%u_old, 6*npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%u_oldold, 6*npmax, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x0, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x_pivot_b, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x_pivot_g, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%phi, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%alpha, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%theta, 1, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vt0, 3, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%at0, 3, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vr0, 3, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%ar0, 3, MPI_DOUBLE_PRECISION, i-1, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%StartupStep, 1, MPI_LOGICAL, i-1, MPI_COMM_WORLD, mpicode )

  enddo

  call toc("Flexible wing (Transfer_solid_data_to_all_CPUs)", MPI_wtime() - t0)
  endif


end subroutine FlexibleSolidSolverWrapper


!-------------------------------------------------------------------------------
!   solid solver main entry point
!-------------------------------------------------------------------------------
subroutine OnlyFlexibleSolidSimulation()
  use vars
  use module_helpers
  use module_ini_files_parser_mpi
  use p3dfft_wrapper
  use penalization, only : mask_color
  use module_insects
  use mpi
  implicit none
  character(len=strlen)  :: infile
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the wings we're using (object oriented)
  type(flexible_wing), dimension(1:nWings) :: Wings
  real(kind=pr) :: t1,t2
  real(kind=pr) :: time, dt0, dt1
  integer :: it!,nsave
  real(kind=pr), allocatable :: work(:,:,:)

  ! in case you forgot to set it
  ! if (itdrag==0) itdrag=10
  if (dt_fixed<=1.0d-10) dt_fixed=1.d-3
  dt0 = dt_fixed
  dt1 = dt_fixed

  write (*,*) "*** information: starting OnlyFlexibleSolidSimulation"

  !-----------------------------------------------------------------------------
  ! Read input parameters and mesh data
  !-----------------------------------------------------------------------------
  if (root) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(2,infile)




  ! prints the profiling on screen
  call summarize_profiling( MPI_COMM_WORLD )

end subroutine
