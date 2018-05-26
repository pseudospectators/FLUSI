subroutine init_beams ( beams )
  !---------------------------------------------------
  ! initializes an array of beams. the initial state is always
  ! straight lines, possible oriented with different angles, at rest.
  !---------------------------------------------------
  implicit none
  integer :: n, i, a
  type(solid), dimension (1:nBeams), intent (out) :: beams
  real (kind=pr) :: alpha, alpha_t, alpha_tt
  ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real (kind=pr), dimension(1:6) :: LeadingEdge
  character(len=strlen) :: fname
  character(len=1)  :: beamstr
  character(len=16) :: frmt

  if (mpirank==0) then
    write (*,'("Initializing ",i1," beams")')  nBeams
    write (*,'("Beam width (2*t_beam) is ",es12.4)') 2.0d0*t_beam
    write (*,'("Beam width (2*t_beam) covers ",(es12.4)," points")') &
    2.0d0*t_beam/max(dy,dz)
  endif

  call lapack_unit_test()

  !-------------------------------------------
  ! allocate beam storage for each beam
  !-------------------------------------------
  do i = 1, nBeams
    !------------------------------------------------
    !-- overwrite beam files if not retaking a backup
    !------------------------------------------------
    if ((index(inicond,'backup::') == 0).and.(root)) then
      write (beamstr,'(i1)') i

      open  (14,file='IBES_iter.t',status='replace')
      write (14,'(21(A15,1x))') "%          time","iterations"
      close (14)

      open  (14,file='beam_data'//beamstr//'.t',status='replace')
      write (14,'(21(A15,1x))') "%          time","press(ns-1)","press(3ns/4)","press(ns/2)",&
      "x(ns-1)","y(ns-1)","vx(ns-1)","vy(ns-1)","theta(ns-1)","theta_dot(ns-1)",&
      "Force(1)","Force(2)","Force_unst(1)","Force_unst(2)","Force_press(1)","Force_press(2)",&
      "Inert_Force(1)","Inert_Force(2)","E_kinetic","E_elastic","E_pot"
      close (14)

      open  (14,file='mouvement'//beamstr//'.t',status='replace')
      write (14,'(10(A15,1x))') "%          time","alpha","alpha_t","alpha_tt","leadingedge_x","leadingedge_y",&
      "leadingedge_vx","leadingedge_vy","leadingedge_ax","leadingedge_ay"
      close (14)

      open (15,file='beam_x'//beamstr//'.t',status='replace')
      open (16,file='beam_y'//beamstr//'.t',status='replace')
      open (17,file='beam_vx'//beamstr//'.t',status='replace')
      open (18,file='beam_vy'//beamstr//'.t',status='replace')
      open (19,file='beam_theta'//beamstr//'.t',status='replace')
      open (20,file='beam_p'//beamstr//'.t',status='replace')

      do n = 15,20
        write(n,'(A15,1x)',advance='no') "%          time"
        do a = 0,ns-1
          write(n,'("     value(",i3,")",1x)',advance='no') a
        enddo
        close(n)
      enddo

      ! call init_empty_file('beam_x'//beamstr//'.t')
      ! call init_empty_file('beam_y'//beamstr//'.t')
      ! call init_empty_file('beam_vx'//beamstr//'.t')
      ! call init_empty_file('beam_vy'//beamstr//'.t')
      ! call init_empty_file('beam_theta'//beamstr//'.t')
      ! call init_empty_file('beam_p'//beamstr//'.t')
    endif

    !---------------------------------------------
    ! define adjustable parameters for each beam
    ! this is position and motion protocoll
    !--------------------------------------------
    beams(i)%x0 = 0.d0 ! always zero, translation and rotation is handled elsewhere
    beams(i)%y0 = 0.d0
    beams(i)%AngleBeam = AngleBeam
    beams(i)%phase = 0.d0

    !--------------------------------------
    !-- initialize beam
    !--------------------------------------
    ! fetch leading edge position
    call mouvement(0.d0, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )
    ! initialize as zero
    beams(i)%x = 0.d0
    beams(i)%y = 0.d0
    beams(i)%vx = 0.d0
    beams(i)%vy = 0.d0
    beams(i)%ax = 0.d0
    beams(i)%ay = 0.d0
    beams(i)%theta = 0.d0
    beams(i)%theta_dot = 0.d0
    beams(i)%theta_old = 0.d0
    beams(i)%theta_dot_old = 0.d0
    beams(i)%theta_oldold = 0.d0
    beams(i)%theta_dot_oldold = 0.d0
    beams(i)%pressure_old = 0.d0
    beams(i)%pressure_new = 0.d0
    beams(i)%tau_old = 0.d0
    beams(i)%tau_new = 0.d0
    beams(i)%Inertial_Force = 0.d0
    ! this is used for unst correction computation:
    beams(i)%drag_unst_new = 0.d0
    beams(i)%drag_unst_old = 0.d0
    beams(i)%lift_unst_new = 0.d0
    beams(i)%lift_unst_old = 0.d0
    beams(i)%Force = 0.d0
    beams(i)%Force_unst = 0.d0
    beams(i)%Force_press = 0.d0
    beams(i)%Inertial_Force = 0.d0
    beams(i)%E_kinetic = 0.d0
    beams(i)%E_elastic = 0.d0
    beams(i)%E_pot = 0.d0
    beams(i)%UnsteadyCorrectionsReady = .false.
    beams(i)%StartupStep = .true.
    beams(i)%dt_old = 0.d0

    ! remember that what follows is described in relative coordinates
    do n = 0, ns-1
      beams(i)%x(n) = LeadingEdge(1) + dble(n)*ds
      beams(i)%y(n) = LeadingEdge(2)
      beams(i)%vx(n)= 0.d0
      beams(i)%vy(n)= 0.d0
      beams(i)%ax(n)= 0.d0
      beams(i)%ay(n)= 0.d0
      ! the grid:
      beams(i)%s(n) = dble(n)*ds
      ! width in rigid (usually span) direction:
      beams(i)%L_rigid(n) = z_top(beams(i)%s(n)) + z_bottom(beams(i)%s(n))
    enddo

    !---------------------------------------------------------------------------
    ! set up material and derivatives
    !---------------------------------------------------------------------------
    beams(i)%zeta(0:ns-1) = eta0 * beams(i)%L_rigid(0:ns-1)
    beams(i)%mu(0:ns-1)   = mue0 * beams(i)%L_rigid(0:ns-1)
    call Differentiate1D (beams(i)%zeta, beams(i)%zeta_s, ns, ds, 1)
    call Differentiate1D (beams(i)%zeta, beams(i)%zeta_ss, ns, ds, 2)
    call Differentiate1D (beams(i)%zeta, beams(i)%zeta_sss, ns, ds, 3)
    call Differentiate1D (beams(i)%mu, beams(i)%mu_s, ns, ds, 1)
    beams(i)%mu_star(0:ns-1) = beams(i)%mu_s(0:ns-1) / beams(i)%mu(0:ns-1)


    if (mpirank ==0) then
      write(*,'(80("-"))')
      write(*,'("setting up solid material with mu=",es12.4," and eta=",es12.4)') mue0, eta0
      write(frmt,'("(",i3.3,"(es12.4,1x))")') ns
      write(*,*) "Length in rigid direction:"
      write(*,frmt) beams(i)%L_rigid(0:ns-1)
      write(*,*) "density coefficient mu:"
      write(*,frmt) beams(i)%mu(0:ns-1)
      write(*,*) "stiffness coefficient eta:"
      write(*,frmt) beams(i)%zeta(0:ns-1)
      write(*,'(80("-"))')
    endif


    if (TimeMethodSolid=="prescribed") then
      if(mpirank==0) write(*,*) "prescribed deformation: initializing..."
      call prescribed_beam ( 0.d0, 0.d0, beams(i) )
    endif

     ! to take static angle into account
    call integrate_position (0.d0, beams(i))
  enddo


  !-------------------------------------------
  ! If we resume a backup, read from file (all ranks do that)
  !-------------------------------------------
  if ( index(inicond,'backup::') /= 0 ) then
    fname = inicond(index(inicond,'::')+2:index(inicond,'.'))//'fsi_bckp'
    call read_solid_backup( beams, trim(adjustl(fname)) )
  endif
end subroutine init_beams
