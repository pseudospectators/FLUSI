subroutine init_wings ( wings )
  !---------------------------------------------------
  ! initializes an array of wings. the initial state is always
  ! straight lines, possible oriented with different angles, at rest.
  !---------------------------------------------------
  implicit none
  integer :: n, i, a
  type(wing), dimension (1:nWings), intent (out) :: wings
  real (kind=pr) :: alpha
  ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  !real (kind=pr), dimension(1:6) :: LeadingEdge
  character(len=strlen) :: fname
  character(len=1)  :: wingstr
  character(len=16) :: frmt

  !if (mpirank==0) then
  !  write (*,'("Initializing ",i1," wings")')  nWings
  !  write (*,'("wing width (2*t_wing) is ",es12.4)') 2.0d0*t_wing
  !  write (*,'("wing width (2*t_wing) covers ",(es12.4)," points")') &
  !  2.0d0*t_wing/max(dy,dz)
  !endif

  !call lapack_unit_test()

  !-------------------------------------------
  ! allocate wing storage for each wing
  !-------------------------------------------
  do i = 1, nWings
    !TO DO Add reading from backup file procedure

    !---------------------------------------------
    ! define adjustable parameters for each wing
    ! this is position and motion protocoll
    !--------------------------------------------
    wings(i)%x0 = 0.d0 ! always zero, translation and rotation is handled elsewhere
    wings(i)%y0 = 0.d0
    wings(i)%z0 = 0.d0
    !rotation angles only used to determine the starting position of the wings
    wings(i)%Anglewing_y = 0.d0
    wings(i)%Anglewing_z = 0.d0

    !--------------------------------------
    !-- initialize wing
    !--------------------------------------
    ! fetch leading edge position
    !call mouvement(0.d0, alpha, alpha_t, alpha_tt, LeadingEdge, wings(i) )
    ! initialize as zero
    wings(i)%x = 0.d0
    wings(i)%y = 0.d0
    wings(i)%z = 0.d0
    wings(i)%vx = 0.d0
    wings(i)%vy = 0.d0
    wings(i)%vz = 0.d0
    wings(i)%u_old = 0.d0
    wings(i)%u_oldold = 0.d0
    wings(i)%tri_elements = 0
    wings(i)%Veins_bending = 0.d0
    wings(i)%Veins_extension = 0.d0
    wings(i)%Veins_bending_BC = 0.d0
    wings(i)%Veins_extension_BC = 0.d0
    wings(i)%Membranes_extension = 0.d0
    wings(i)%Membrane_edge = 0.d0
    !wings(i)%pressure_old = 0.d0
    !wings(i)%pressure_new = 0.d0
    ! this is used for unst correction computation:
    !wings(i)%drag_unst_new = 0.d0
    !wings(i)%drag_unst_old = 0.d0
    !wings(i)%lift_unst_new = 0.d0
    !wings(i)%lift_unst_old = 0.d0
    !wings(i)%Force = 0.d0
    !wings(i)%Force_unst = 0.d0
    !wings(i)%Force_press = 0.d0
    !wings(i)%Inertial_Force = 0.d0
    !wings(i)%E_kinetic = 0.d0
    !wings(i)%E_elastic = 0.d0
    !wings(i)%E_tot = 0.d0
    !wings(i)%UnsteadyCorrectionsReady = .false.
    wings(i)%StartupStep = .true.
    wings(i)%dt_old = 0.d0

    !---------------------------------------------------------------------------
    ! set up material and derivatives
    !---------------------------------------------------------------------------
    wings(i)%EIy = 0.d0
    wings(i)%EIz = 0.d0
    wings(i)%kby0 = 0.d0
    wings(i)%kbz0 = 0.d0
    wings(i)%EIy_BC = 0.d0
    wings(i)%EIz_BC = 0.d0
    wings(i)%kby0_BC = 0.d0
    wings(i)%kbz0_BC = 0.d0
    wings(i)%ke_v = 0.d0
    wings(i)%kb_v = 0.d0
    wings(i)%ke_vBC = 0.d0
    wings(i)%kb_vBC = 0.d0
    wings(i)%kb_m = 0.d0
    wings(i)%kb_me = 0.d0
    !wings(i)%zeta(0:ns-1) = eta0 * wings(i)%L_rigid(0:ns-1)
    !wings(i)%mu(0:ns-1)   = mue0 * wings(i)%L_rigid(0:ns-1)

    ! Reading mesh data from ASCII files
    call read_mesh_data(wings)

    !if (mpirank ==0) then
    !  write(*,'(80("-"))')
    !  write(*,'("setting up solid material with mu=",es12.4," and eta=",es12.4)') mue0, eta0
    !  write(frmt,'("(",i3.3,"(es12.4,1x))")') ns
    !  write(*,*) "Length in rigid direction:"
    !  write(*,frmt) wings(i)%L_rigid(0:ns-1)
    !  write(*,*) "density coefficient mu:"
    !  write(*,frmt) wings(i)%mu(0:ns-1)
    !  write(*,*) "stiffness coefficient eta:"
    !  write(*,frmt) wings(i)%zeta(0:ns-1)
    !  write(*,'(80("-"))')
    !endif


    !if (TimeMethodSolid=="prescribed") then
    !  if(mpirank==0) write(*,*) "prescribed deformation: initializing..."
    !  call prescribed_wing ( 0.d0, 0.d0, wings(i) )
    !endif

     ! to take static angle into account
    !call integrate_position (0.d0, wings(i))
  enddo

  !-------------------------------------------
  ! If we resume a backup, read from file (all ranks do that)
  !-------------------------------------------
  !if ( index(inicond,'backup::') /= 0 ) then
  !  fname = inicond(index(inicond,'::')+2:index(inicond,'.'))//'fsi_bckp'
  !  call read_solid_backup( wings, trim(adjustl(fname)) )
  !endif

  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Flexible wings initialization is complete."
    write(*,'(80("<"))')
  endif

end subroutine init_wings

subroutine read_mesh_data(wings)

  use vars
  implicit none
  type (wing), dimension(1:nWings), intent (inout) :: wings
  character(len=strlen) :: data_file
  character(len=1)  :: wingstr
  integer :: i, j
  real(kind=pr), allocatable :: tmp(:,:)

  do i=1, nWings

    !-- for naming files..
    write (wingstr,'(i1)') i

    ! Read initial coordinates x,y,z of all points in 2nd,3rd,4th columms
    ! respectively in the points_coor.t data file
    data_file = 'points_coor'//wingstr//'.t'
    call  read_mesh_data_2D_array(data_file, tmp)
    do j=1, maxval(tmp(:,1))
      wings(i)%x(j) = tmp(j,2)
      wings(i)%y(j) = tmp(j,3)
      wings(i)%z(j) = tmp(j,4)
    end do

    ! Read indices of three vertices (correnponding to 3rd, 4th and 5tn columms)
    ! of all triangle elements of the mesh
    data_file = 'mesh_triangle_elements'//wingstr//'.t'
    call  read_mesh_data_2D_array(data_file, tmp)
    do j=1, size(tmp,DIM=1)
      wings(i)%tri_elements(j,1) = j
      wings(i)%tri_elements(j,2) = int(tmp(j,3))
      wings(i)%tri_elements(j,3) = int(tmp(j,4))
      wings(i)%tri_elements(j,4) = int(tmp(j,5))
    end do

  end do

end subroutine read_mesh_data

subroutine read_mesh_data_2D_array(data_file, data_2D_array)

character(len=strlen), intent(in) :: data_file
real(kind=pr),allocatable,intent(inout) :: data_2D_array(:,:)
integer :: num_lines, num_cols, n_header=0

call count_lines_in_ascii_file_mpi(data_file, num_lines, n_header)
call count_cols_in_ascii_file_mpi(data_file, num_cols, n_header)
allocate(data_2D_array(1:num_lines, 1:num_cols) )
call read_array_from_ascii_file_mpi(data_file, data_2D_array, n_header)

end subroutine read_mesh_data_2D_array
