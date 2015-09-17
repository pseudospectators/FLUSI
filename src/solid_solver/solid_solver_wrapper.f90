!-------------------------------------------------------------------------------
!   SOLID SOLVER WRAPPER
! Input:
!   time: current time level (t^n)
!   dt: time step for new level (t^n+1 - t^n)
!   Beams: array of beams, all at the old time level t^n
! Output:
!   Beams: array of beams, all at the new time level t^n+1
!-------------------------------------------------------------------------------
subroutine SolidSolverWrapper ( time, dt, beams )
  use penalization ! mask array etc
  implicit none
  real(kind=pr), intent (in) ::  dt, time
  real(kind=pr) :: t0
  type(solid), dimension(1:nBeams), intent (inout) ::    beams
  integer :: i
  t0 = MPI_wtime()

  do i = 1, nBeams
      !------------------------------------------
      ! check if input values are okay
      !------------------------------------------
      if (Vector_isNAN(beams(i)%pressure_new).or.&
          Vector_isNAN(beams(i)%pressure_old).or.&
          Vector_isNAN(beams(i)%tau_new).or.&
          Vector_isNAN(beams(i)%tau_old) ) then
        if (root) write(*,*) "SolidSolver: input values contain NaNs"
        ! at this occasion, check if mask and us contain NaNs
        if (root) write(*,*) "Checking mask and us for NaNs"
        call checknan( mask, "mask")
        call checknan( us(:,:,:,1), "usx")
        call checknan( us(:,:,:,2), "usy")
        call checknan( us(:,:,:,3), "usz")
        ! time to go..
        call abort("The input forces for the solid solver contain NaN..abort")
      endif

      if (time>=T_release) then
        !-------------------------------------------
        ! the beams are released (now: active FSI), call IBES solvers
        !-------------------------------------------
        select case (TimeMethodSolid)
          case ("CN2","BDF2","EI1")
            ! all implicit solvers are in one subroutine
            call IBES_solver (time, dt, beams(i))
          case ("RK4")
            ! note this solver may be very unstable:
            ! call RK4_wrapper (time, dt, beams(i))
            call abort("RK4 for beam equation is not implemented")
          case ("EE1")
            ! note this solver may be very unstable:
            ! call EE1_wrapper (time, dt, beams(i))
            call abort("EE1 for beam equation is not implemented")
          case ("prescribed")
            ! this is not a solver, but for passive FSI with prescribed deformation:
            call prescribed_beam (time, dt, beams(i))
          case default
            call abort("SolidSolverWrapper::invalid value of TimeMethodSolid"//&
                 trim(adjustl(TimeMethodSolid)))
        end select
      else
        !-------------------------------------------
        ! the beams are not yet released, but their leading edges may move
        ! (passive FSI)
        !-------------------------------------------
        call integrate_position (time+dt, beams(i))
      endif

      !-------------------------------------------
      ! compute energies and stuff
      !-------------------------------------------
      call SolidEnergies ( beams(i) )

      !-- check if everything seems okay, if not show beam and abort
      call show_beam_on_error( beams(i) )
  enddo
  time_solid = time_solid + MPI_wtime() - t0
end subroutine SolidSolverWrapper
