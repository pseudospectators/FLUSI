subroutine init_wings ( fname, wings, dx_reference)
  !---------------------------------------------------
  ! initializes an array of wings. the initial state is always
  ! straight lines, possible oriented with different angles, at rest.
  !---------------------------------------------------
  implicit none
  integer :: n, i, a, j, ind, itri
  character(len=strlen), intent(in) :: fname
  type(flexible_wing), dimension (1:nWings), intent (inout) :: Wings
  real(kind=pr), intent(in) :: dx_reference
  character(len=strlen) :: filename
  real(kind=pr) :: alpha
  real(kind=pr) :: delta(1:3)
  real(kind=pr), dimension(1:3) :: u
  real(kind=pr), allocatable :: normal(:,:)

  type(inifile) :: PARAMS
  ! LeadingEdge: x, y, vx, vy, ax, ay (Array)
  !real (kind=pr), dimension(1:6) :: LeadingEdge
  character(len=1)  :: wingstr
  character(len=16) :: frmt


  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Initializing flexible wing module!"
    write(*,*) "*.ini file is: "//trim(adjustl(fname))
    write(*,'(80("<"))')
  endif


  !-------------------------------------------
  ! allocate wing storage for each wing
  !-------------------------------------------
    !TODO Add reading from backup file procedure

    do i = 1, nWings
    !---------------------------------------------
    ! define adjustable parameters for each wing
    ! this is position and motion protocoll
    !--------------------------------------------

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
    wings(i)%tri_element_areas = 0.d0
    wings(i)%tri_element_normals = 0.d0
    wings(i)%Veins_bending = 0.d0
    wings(i)%Veins_extension = 0.d0
    wings(i)%Veins_bending_BC = 0.d0
    wings(i)%Veins_extension_BC = 0.d0
    wings(i)%Membranes_extension = 0.d0
    wings(i)%Membrane_edge = 0.d0
    wings(i)%m=0.d0
    wings(i)%c=0.d0
    wings(i)%StartupStep = .true.
    wings(i)%dt_old = 0.d0
    wings(i)%press_upside = 0.d0
    wings(i)%press_downside = 0.d0
    wings(i)%x0 = 0.d0
    wings(i)%y0 = 0.d0
    wings(i)%z0 = 0.d0
    wings(i)%WingAngle_x = 0.d0
    wings(i)%WingAngle_y = 0.d0
    wings(i)%WingAngle_z = 0.d0
    wings(i)%RotationMat_x = 0.d0
    wings(i)%RotationMat_y = 0.d0
    wings(i)%RotationMat_z = 0.d0
    wings(i)%vt0 = 0.d0
    wings(i)%at0 = 0.d0
    wings(i)%vr0 = 0.d0
    wings(i)%ar0 = 0.d0

    ! Reading mesh data from ASCII files and save them into the state vector of
    ! the mass-spring model u_old
    call read_wing_mesh_data(wings(i), i)


    !-----------------------------------------------------------------------------
    ! read in parameters form ini file
    !-----------------------------------------------------------------------------
    ! read in the complete ini file, from which we initialize the flexible wings
    call read_ini_file_mpi(PARAMS, fname, verbose=.true.)

    call read_param_mpi(PARAMS,"Flexible_wing","x0",wings(i)%x0, 0.d0)
    call read_param_mpi(PARAMS,"Flexible_wing","y0",wings(i)%y0, 0.d0)
    call read_param_mpi(PARAMS,"Flexible_wing","z0",wings(i)%z0, 0.d0)
    !call read_param_mpi(PARAMS,"Flexible_wing","v0",wings(i)%v0, (/0.d0, 0.d0, 0.d0/))
    call read_param_mpi(PARAMS,"Flexible_wing","t_wing",wings(i)%t_wing, 2.0d0*dx_reference)
    call read_param_mpi(PARAMS,"Flexible_wing","wing_smoothing",wings(i)%wing_smoothing, 1.0d0*dx_reference)

    call read_param_mpi(PARAMS,"Flexible_wing","EIy",wings(i)%EIy)
    call read_param_mpi(PARAMS,"Flexible_wing","EIz",wings(i)%EIz)
    call read_param_mpi(PARAMS,"Flexible_wing","EIy_with_BC",wings(i)%EIy_BC)
    call read_param_mpi(PARAMS,"Flexible_wing","EIz_with_BC",wings(i)%EIz_BC)

    call read_param_mpi(PARAMS,"Flexible_wing","ke_veins",wings(i)%ke0_v)
    call read_param_mpi(PARAMS,"Flexible_wing","ke_veins_with_BC",wings(i)%ke0_vBC)
    call read_param_mpi(PARAMS,"Flexible_wing","ke_membranes",wings(i)%ke0_m)

    if (load_mass_from_file == 'no') then
    call read_param_mpi(PARAMS,"Flexible_wing","density_veins",wings(i)%rho_v)
    call read_param_mpi(PARAMS,"Flexible_wing","density_veins_with_BC",wings(i)%rho_vBC)
    call read_param_mpi(PARAMS,"Flexible_wing","density_membranes",wings(i)%rho_m)
    endif

    call read_param_mpi(PARAMS,"Flexible_wing","damping",wings(i)%c0, 0.d0)

    call read_param_mpi(PARAMS,"Flexible_wing","Wing_angle_x",wings(i)%WingAngle_x, 0.d0)
    call read_param_mpi(PARAMS,"Flexible_wing","Wing_angle_y",wings(i)%WingAngle_y, 0.d0)
    call read_param_mpi(PARAMS,"Flexible_wing","Wing_angle_z",wings(i)%WingAngle_z, 0.d0)

    call read_param_mpi(PARAMS,"Flexible_wing","Motion",wings(i)%Motion,"stationary")

    call read_param_mpi(PARAMS,"Flexible_wing","Gravity",grav, (/0.d0, 0.d0, -9.8d0/))
    call read_param_mpi(PARAMS,"Flexible_wing","TimeMethodFlexibleSolid",TimeMethodFlexibleSolid,"EI1")

    call read_param_mpi(PARAMS,"Flexible_wing","T_release",wings(i)%T_release,0.0d0)
    call read_param_mpi(PARAMS,"Flexible_wing","tau",wings(i)%tau,0.0d0)
    call read_param_mpi(PARAMS,"Flexible_wing","ControlPoint",wings(i)%ControlPoint,0)


    ! clean ini file
    call clean_ini_file_mpi(PARAMS)

    !--------------------------------------------------------------------------
    ! Move the wing to the desired position X0
    !--------------------------------------------------------------------------
    wings(i)%x(1:wings(i)%np) = wings(i)%u_old(1:wings(i)%np)
    wings(i)%y(1:wings(i)%np) = wings(i)%u_old(wings(i)%np+1:2*wings(i)%np)
    wings(i)%z(1:wings(i)%np) = wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np)

    call rotate_wing(wings(i))

    call translate_wing(wings(i))

   if (root) then
           write(*,*) maxval(wings(i)%x), minval(wings(i)%x)
	         write(*,*) maxval(wings(i)%y), minval(wings(i)%y)
	         write(*,*) maxval(wings(i)%z), minval(wings(i)%z)
    endif

    call determine_boundary_points_from_origin(wings(i))

    !--------------------------------------------------------------------------
    ! Determine initial geometrical properties of the wings: initial lengths,
    !  angles of springs and orientation of the wings
    !--------------------------------------------------------------------------
    allocate(normal(1:wings(i)%ntri,1:3))
    do itri=1,wings(i)%ntri
        ! Calculate the normal vector of one triangle
        normal(itri,1:3) = cross((/wings(i)%x(wings(i)%tri_elements(itri,2)) - &
                                   wings(i)%x(wings(i)%tri_elements(itri,3)),  &
                                   wings(i)%y(wings(i)%tri_elements(itri,2)) - &
                                   wings(i)%y(wings(i)%tri_elements(itri,3)),  &
                                   wings(i)%z(wings(i)%tri_elements(itri,2)) - &
                                   wings(i)%z(wings(i)%tri_elements(itri,3))/),&
                                 (/wings(i)%x(wings(i)%tri_elements(itri,3)) - &
                                   wings(i)%x(wings(i)%tri_elements(itri,4)),  &
                                   wings(i)%y(wings(i)%tri_elements(itri,3)) - &
                                   wings(i)%y(wings(i)%tri_elements(itri,4)),  &
                                   wings(i)%z(wings(i)%tri_elements(itri,3)) - &
                                   wings(i)%z(wings(i)%tri_elements(itri,4))/))

        ! dimentionalized to get a unit vector
        wings(i)%tri_element_normals(itri,1) = normal(itri,1)/norm2(normal(itri,1:3))
        wings(i)%tri_element_normals(itri,2) = normal(itri,2)/norm2(normal(itri,1:3))
        wings(i)%tri_element_normals(itri,3) = normal(itri,3)/norm2(normal(itri,1:3))

        !Calculate area of triangle elements
        wings(i)%tri_element_areas(itri) = 0.5*norm2(normal(itri,1:3))

        ! Check the orientation of the normal vectors comparing with Oz axis. This is
        ! done only at the first time step of the simulation.
        if (dot_product(wings(i)%tri_element_normals(itri,1:3),(/0.0d0,1.0d0,0.0d0/))<-1.0d-10) then
            wings(i)%tri_element_normals(itri,4) = -1
        elseif (dot_product(wings(i)%tri_element_normals(itri,1:3),(/0.0d0,1.0d0,0.0d0/))>1.0d-10) then
            wings(i)%tri_element_normals(itri,4) =  1
        else
            call abort(1412, "Wing normal vector is perpendicular with the Oz axis. &
                              The wing should be placed on the Oxy plane for the best performance of the solver.")
        endif
    enddo
    deallocate(normal)



    do j=1,nMembranes
        call length_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                                wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                                wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np),    &
                                wings(i)%membranes_extension(:,:,j))

        wings(i)%membranes_extension(:,4,j) = wings(i)%membranes_extension(:,5,j)
    enddo

    do j=1,nMembrane_edges
        call length_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                            wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                            wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np), &
                            wings(i)%membrane_edge(:,:,j))

        wings(i)%membrane_edge(:,4,j) = wings(i)%membrane_edge(:,5,j)
    enddo

    do j=1,nVeins

        call length_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                                wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                                wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np), &
                                wings(i)%veins_extension(:,:,j))


        call angle_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                               wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                               wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np), &
                               wings(i)%veins_bending(:,:,j))

        wings(i)%veins_extension(:,4,j) = wings(i)%veins_extension(:,5,j)
        wings(i)%veins_bending(:,5,j) = wings(i)%veins_bending(:,7,j)
        wings(i)%veins_bending(:,6,j) = wings(i)%veins_bending(:,8,j)

    enddo

    do j=1,nVeins_BC
        call length_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                                wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                                wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np),   &
                                wings(i)%veins_extension_BC(1:,:,j))
        call angle_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                               wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                               wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np), &
                               wings(i)%veins_bending_BC(1:,:,j))

        wings(i)%veins_extension_BC(1:,4,j) = wings(i)%veins_extension_BC(1:,5,j)
        wings(i)%veins_bending_BC(1:,5,j) = wings(i)%veins_bending_BC(1:,7,j)
        wings(i)%veins_bending_BC(1:,6,j) = wings(i)%veins_bending_BC(1:,8,j)
    end do

    call angle_calculation_wrapper(wings(i)%u_old(1:wings(i)%np), &
                           wings(i)%u_old(wings(i)%np+1:2*wings(i)%np), &
                           wings(i)%u_old(2*wings(i)%np+1:3*wings(i)%np), &
                           wings(i)%vein_connectors)

        wings(i)%vein_connectors(1:,5) = wings(i)%vein_connectors(1:,7)
        wings(i)%vein_connectors(1:,6) = wings(i)%vein_connectors(1:,8)

    !--------------------------------------------------------------------------
    ! Set up material properties
    !--------------------------------------------------------------------------

    do j=1,nMembranes
    wings(i)%ke_m(:,j) = wings(i)%ke0_m(j)
      do ind=1,nint(maxval(wings(i)%membranes(:,1,j)))
          wings(i)%c(nint(wings(i)%membranes(ind,2,j))) = 2.d-3
      enddo



      if (load_mass_from_file == 'no') then
      do ind=1,nint(maxval(wings(i)%membranes(:,1,j)))
          wings(i)%m(nint(wings(i)%membranes(ind,2,j))) = wings(i)%rho_m(j)
      enddo
      endif
    enddo

    do j=1,nMembrane_edges
      wings(i)%ke_me(:,j) = wings(i)%ke0_m(1)
    enddo

    do j=1,nVeins
      call convert_flexural_rigidity_into_spring_stiffness(wings(i)%EIy(j), wings(i)%EIz(j),  &
                                                          wings(i)%kby0(j), wings(i)%kbz0(j), &
                                                          wings(i)%veins_extension(:,:,j))

      wings(i)%kby(:,j) = wings(i)%kby0(j)
      wings(i)%kbz(:,j) = wings(i)%kbz0(j)
      wings(i)%ke_v(:,j) = wings(i)%ke0_v(j)
      if (load_mass_from_file == 'no') then
      do ind=1,nint(maxval(wings(i)%veins(:,1,j)))
          wings(i)%m(nint(wings(i)%veins(ind,2,j))) = wings(i)%rho_v(j)
      enddo
      endif
    enddo

    do j=1,nVeins_BC
      call convert_flexural_rigidity_into_spring_stiffness(wings(i)%EIy_BC(j), wings(i)%EIz_BC(j),  &
                                                          wings(i)%kby0_BC(j), wings(i)%kbz0_BC(j), &
                                                          wings(i)%veins_extension_BC(1:,:,j))

      wings(i)%kby_BC(:,j) = wings(i)%kby0_BC(j)
      wings(i)%kbz_BC(:,j) = wings(i)%kbz0_BC(j)
      wings(i)%ke_vBC(:,j) = wings(i)%ke0_vBC(j)
      wings(i)%kby_BC(-1,j) = 50.d0*wings(i)%kby_BC(1,j)
      wings(i)%kby_BC(0,j) = 50.d0*wings(i)%kby_BC(1,j)
      wings(i)%kbz_BC(-1,j) = 50.d0*wings(i)%kbz_BC(1,j)
      wings(i)%kbz_BC(0,j) = 50.d0*wings(i)%kbz_BC(1,j)
      if (load_mass_from_file == 'no') then
      do ind=1,nint(maxval(wings(i)%veins_BC(:,1,j)))
          wings(i)%m(nint(wings(i)%veins_BC(ind,2,j))) = wings(i)%rho_vBC(j)
      enddo
      endif
    enddo

    !TODO: set stiffness according to which veins they connect
    wings(i)%kby_c(1:nvmax) = wings(i)%kby_BC(1:nvmax,1)
    wings(i)%kbz_c(1:nvmax) = wings(i)%kbz_BC(1:nvmax,1)

    if (mpirank ==0) then
      write(*,'(80("-"))')
      write(*,'("Setting up material properties for the wing number ",i2.2," with")') i
      write(frmt,'("(",i3.3,"(es12.4,1x))")') wings(i)%np
      write(*,*) "Mass points:"
      write(*,frmt) wings(i)%m(1:wings(i)%np)
      write(*,*) "Damping coeficients:"
      write(*,frmt) wings(i)%c(1:wings(i)%np)
      do j=1,nVeins_BC
        write(frmt,'("(",i3.3,"(es12.4,1x))")') nint(maxval(wings(i)%veins_bending_BC(:,1,j)))+2
        write(*,'("bending stiffness of y-direction bending springs of the vein with BC number ",i2.2,":")',advance='yes') j
        write(*,frmt) wings(i)%kby_BC(-1:nint(maxval(wings(i)%veins_bending_BC(:,1,j))),j)
        write(*,'("bending stiffness of z-direction bending springs of the vein with BC number ",i2.2,":")') j
        write(*,frmt) wings(i)%kbz_BC(-1:nint(maxval(wings(i)%veins_bending_BC(:,1,j))),j)
        write(frmt,'("(",i3.3,"(es12.4,1x))")') nint(maxval(wings(i)%veins_extension_BC(:,1,j)))+1
        write(*,'("extension stiffness of extension springs of the vein with BC number ",i2.2,":")') j
        write(*,frmt) wings(i)%ke_vBC(0:nint(maxval(wings(i)%veins_extension_BC(:,1,j))),j)
      enddo
      do j=1,nVeins
        write(frmt,'("(",i3.3,"(es12.4,1x))")') nint(maxval(wings(i)%veins_bending(:,1,j)))
        write(*,'("bending stiffness of y-direction bending springs of the vein number ",i2.2,":")',advance='yes') j
        write(*,frmt) wings(i)%kby(1:nint(maxval(wings(i)%veins_bending(:,1,j))),j)
        write(*,'("bending stiffness of z-direction bending springs of the vein number ",i2.2,":")') j
        write(*,frmt) wings(i)%kbz(1:nint(maxval(wings(i)%veins_bending(:,1,j))),j)
        write(frmt,'("(",i3.3,"(es12.4,1x))")') nint(maxval(wings(i)%veins_extension(:,1,j)))
        write(*,'("extension stiffness of extension springs of the vein number ",i2.2,":")') j
        write(*,frmt) wings(i)%ke_v(1:nint(maxval(wings(i)%veins_extension(:,1,j))),j)
      enddo
      !do j=1,nMembranes
        !write(*,*) nint(maxval(wings(i)%membranes_extension(:,1,j)))
      !  write(frmt,'("(",i6.6,"(es12.4,1x))")') nint(maxval(wings(i)%membranes_extension(:,1,j)))
      !  write(*,'("extension stiffness of extension springs of the membrane number ",i2.2,":")',advance='yes') j
      !  write(*,frmt) wings(i)%ke_m(1:nint(maxval(wings(i)%membranes_extension(:,1,j))),j)
      !enddo
        write(frmt,'("(",i3.3,"(es12.4,1x))")') nint(maxval(wings(i)%vein_connectors(:,1)))
        write(*,*) "bending stiffness of y-direction bending springs of the vein connectors:"
        write(*,frmt) wings(i)%kby_c(1:nint(maxval(wings(i)%vein_connectors(:,1))))
        write(*,*) "bending stiffness of z-direction bending springs of the vein connectors:"
        write(*,frmt) wings(i)%kbz_c(1:nint(maxval(wings(i)%vein_connectors(:,1))))

    endif

  enddo

  !-------------------------------------------
  ! If we resume a backup, read from file (all ranks do that)
  !-------------------------------------------
  if ( index(inicond,'backup::') /= 0 ) then
    filename = inicond(index(inicond,'::')+2:index(inicond,'.'))//'fsi_bckp'
    call read_flexible_wing_backup( wings, trim(adjustl(filename)) )
  endif

  if (root) then
    write(*,'(80("<"))')
    write(*,*) "Flexible wings initialization is complete."
    write(*,'(80("<"))')
  endif

end subroutine init_wings

subroutine read_wing_mesh_data(wing, i)

  use vars
  implicit none
  integer, intent(in) :: i !ordinal number of the current wing
  type(flexible_wing), intent (inout) :: wing !for the ith wing
  character(len=strlen) :: data_file
  character(len=1)  :: wingstr
  integer :: j, np
  real(kind=pr), allocatable :: tmp2D(:,:)
  real(kind=pr), allocatable :: tmp1D(:,:)

    !-- for naming files..
    write (wingstr,'(i1)') i

    ! Read initial coordinates x,y,z of all points in 2nd,3rd,4th columms
    ! respectively in the points_coor.t data file
        data_file = 'points_coor'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        ! Saving number of points
        wing%np = nint(maxval(tmp2D(:,1)))
        np = wing%np

        do j=1, nint(maxval(tmp2D(:,1)))
          wing%u_old(j) = tmp2D(j,2)!*2*pi/xl*20
          wing%u_old(j + np) = tmp2D(j,3)!*2*pi/yl*20
          wing%u_old(j + 2*np) = tmp2D(j,4)!*2*pi/zl
        end do

        deallocate(tmp2D)
    ! Read indices of three vertices (correnponding to 3rd, 4th and 5tn columms)
    ! of all triangle elements of the mesh
        data_file = 'mesh_triangle_elements'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        do j=1, size(tmp2D,DIM=1)
          wing%tri_elements(j,1) = j
          wing%tri_elements(j,2) = int(tmp2D(j,3))
          wing%tri_elements(j,3) = int(tmp2D(j,4))
          wing%tri_elements(j,4) = int(tmp2D(j,5))
        end do

        deallocate(tmp2D)
        ! Saving number of triangle elements
        wing%ntri = maxval(wing%tri_elements(:,1))

    ! Read identification numbers of all points belonging to veins
        data_file = 'veins'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%veins(1:int((size(tmp2D,DIM=1))*(1.0/nVeins)),1:2,1:nVeins) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins)),2,nVeins/))

        deallocate(tmp2D)

    ! Read bending springs information of veins without boundary conditions
        data_file = 'veins_bending'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%veins_bending(1:int((size(tmp2D,DIM=1))*(1.0/nVeins)),1:8,1:nVeins) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins)),8,nVeins/))

        deallocate(tmp2D)

     ! Read bending springs information of veins without boundary conditions
         data_file = 'vein_connectors'//wingstr//'.dat'
         call  read_mesh_data_2D_array(data_file, tmp2D)

         wing%vein_connectors(1:int((size(tmp2D,DIM=1))),1:8) = tmp2D

         deallocate(tmp2D)

     ! Read identification numbers of all points belonging to veins with BCs
         data_file = 'veins_BC'//wingstr//'.dat'
         call  read_mesh_data_2D_array(data_file, tmp2D)

         wing%veins_BC(1:int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),1:2,1:nVeins_BC) = &
         reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),2,nVeins_BC/))

         deallocate(tmp2D)

     ! Read bending springs information of veins with boundary conditions
        data_file = 'veins_bending_BC'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%veins_bending_BC(1:int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),1:8,1:nVeins_BC) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),8,nVeins_BC/))

        deallocate(tmp2D)

     ! Read extension springs information of veins without boundary conditions
        data_file = 'veins_extension'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%veins_extension(1:int((size(tmp2D,DIM=1))*(1.0/nVeins)),1:5,1:nVeins) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins)),5,nVeins/))

        deallocate(tmp2D)

     ! Read extension springs information of veins with boundary conditions
        data_file = 'veins_extension_BC'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%veins_extension_BC(1:int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),1:5,1:nVeins_BC) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nVeins_BC)),5,nVeins_BC/))

        deallocate(tmp2D)

      ! TODO change back to general case when we have nMembranes membranes with 2D array
      ! Read identification numbers of all points belonging to membranes
         data_file = 'membranes'//wingstr//'.dat'
         call  read_mesh_data_1D_array(data_file, tmp1D)

         do j=1,int((size(tmp1D)))
         wing%membranes(j,2,1) = tmp1D(j,1)
         wing%membranes(j,1,1) = j
        enddo

         deallocate(tmp1D)

      ! Read extension springs information of membranes
         data_file = 'membranes_extension'//wingstr//'.dat'
         call  read_mesh_data_2D_array(data_file, tmp2D)

         wing%membranes_extension(1:int((size(tmp2D,DIM=1))*(1.0/nMembranes)),1:5,1:nMembranes) = &
         reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nMembranes)),5,nMembranes/))

         deallocate(tmp2D)

      ! Read extension springs information of the edge of the wing
        data_file = 'membrane_edge'//wingstr//'.dat'
        call  read_mesh_data_2D_array(data_file, tmp2D)

        wing%membrane_edge(1:int((size(tmp2D,DIM=1))*(1.0/nMembrane_edges)),1:5,1:nMembrane_edges) = &
        reshape(tmp2D,(/int((size(tmp2D,DIM=1))*(1.0/nMembrane_edges)),5,nMembrane_edges/))

        deallocate(tmp2D)

      ! Read mass
        if (load_mass_from_file == 'yes') then
          data_file = 'mass'//wingstr//'.dat'
          call  read_mesh_data_1D_array(data_file, tmp1D)

          do j=1,int((size(tmp1D)))
          wing%m(j) = tmp1D(j,1)
          !if (root) write(*,*) tmp1D(j,1)
          enddo

          deallocate(tmp1D)
      endif


end subroutine read_wing_mesh_data

subroutine read_mesh_data_1D_array(data_file, data_1D_array)

character(len=strlen), intent(in) :: data_file
real(kind=pr),allocatable,intent(inout) :: data_1D_array(:,:)
integer :: num_lines, n_header=0

call count_lines_in_ascii_file_mpi(data_file, num_lines, n_header)
allocate(data_1D_array(1:num_lines,1) )
call read_array_from_ascii_file_mpi(data_file, data_1D_array, n_header)

end subroutine read_mesh_data_1D_array

subroutine read_mesh_data_2D_array(data_file, data_2D_array)

character(len=strlen), intent(in) :: data_file
real(kind=pr),allocatable,intent(inout) :: data_2D_array(:,:)
integer :: num_lines, num_cols, n_header=0

call count_lines_in_ascii_file_mpi(data_file, num_lines, n_header)
call count_cols_in_ascii_file_mpi(data_file, num_cols, n_header)
allocate(data_2D_array(1:num_lines, 1:num_cols) )
call read_array_from_ascii_file_mpi(data_file, data_2D_array, n_header)

end subroutine read_mesh_data_2D_array

subroutine determine_boundary_points_from_origin(wing)

  implicit none
  type(flexible_wing), intent (inout) :: wing
  integer :: i, np
  real(kind=pr), dimension(1:3) :: delta

  np = wing%np

  ! Calculate the second boundary point for the Leading edge vein from the first
  ! point which is read from param file since the first point of the LE vein is
  ! where we define the root of the wing at (0, 0, 0)
      wing%x_BC(-1,1) = 0.d0
      wing%y_BC(-1,1) = 0.d0
      wing%z_BC(-1,1) = 0.d0

      wing%x0_BC(-1,1) = wing%x_BC(-1,1)
      wing%y0_BC(-1,1) = wing%y_BC(-1,1)
      wing%z0_BC(-1,1) = wing%z_BC(-1,1)

      wing%x_BC(0,1) = wing%u_old(nint(wing%veins_bending_BC(1,2,1)))/2
      wing%y_BC(0,1) = wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + np)/2
      wing%z_BC(0,1) = wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + 2*np)/2

      wing%x0_BC(0,1) = wing%x_BC(0,1)
      wing%y0_BC(0,1) = wing%y_BC(0,1)
      wing%z0_BC(0,1) = wing%z_BC(0,1)

      wing%veins_extension_BC(0,4,1) = sqrt(((wing%u_old(nint(wing%veins_bending_BC(1,2,1))))/2)**2 + &
                                            ((wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + np))/2)**2 + &
                                            ((wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + 2*np))/2)**2)

      ! Calculate initial angles
      call angle_calculation(wing%x_BC(0,1),wing%u_old(nint(wing%veins_bending_BC(1,2,1))), &
                             wing%u_old(nint(wing%veins_bending_BC(1,3,1))), wing%y_BC(0,1),&
                             wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + np), &
                             wing%u_old(nint(wing%veins_bending_BC(1,3,1)) + np), &
                             wing%z_BC(0,1),wing%u_old(nint(wing%veins_bending_BC(1,2,1)) + 2*np), &
                             wing%u_old(nint(wing%veins_bending_BC(1,3,1)) + 2*np), &
                             wing%veins_bending_BC(0,5,1),wing%veins_bending_BC(0,6,1))


  ! Calculate boundary for other veins
      do i=2,nVeins_BC

        delta(1) = wing%u_old(nint(wing%veins_bending_BC(1,3,i))) - &
                       wing%u_old(nint(wing%veins_bending_BC(1,2,i)))
        delta(2) = wing%u_old(nint(wing%veins_bending_BC(1,3,i)) + np) - &
                       wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + np)
        delta(3) = wing%u_old(nint(wing%veins_bending_BC(1,3,i)) + 2*np) - &
                       wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + 2*np)

        wing%x_BC(-1,i) = wing%u_old(nint(wing%veins_bending_BC(1,2,i))) - 2*delta(1)
        wing%y_BC(-1,i) = wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + np) - 2*delta(2)
        wing%z_BC(-1,i) = wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + 2*np) - 2*delta(3)

        wing%x0_BC(-1,i) = wing%x_BC(-1,i)
        wing%y0_BC(-1,i) = wing%y_BC(-1,i)
        wing%z0_BC(-1,i) = wing%z_BC(-1,i)

        wing%x_BC(0,i) = (wing%x_BC(-1,i) + wing%u_old(nint(wing%veins_bending_BC(1,2,i))))/2
        wing%y_BC(0,i) = (wing%y_BC(-1,i) + wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + np))/2
        wing%z_BC(0,i) = (wing%z_BC(-1,i) + wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + 2*np))/2

        wing%x0_BC(0,i) = wing%x_BC(0,i)
        wing%y0_BC(0,i) = wing%y_BC(0,i)
        wing%z0_BC(0,i) = wing%z_BC(0,i)

        ! Calculate initial lengths of springs connecting veins with the BC
        wing%veins_extension_BC(0,4,i) = sqrt((delta(1))**2 + (delta(2))**2 + (delta(3))**2)

        ! Calculate initial angles
        call angle_calculation(wing%x_BC(0,i),wing%u_old(nint(wing%veins_bending_BC(1,2,i))), &
                               wing%u_old(nint(wing%veins_bending_BC(1,3,i))), wing%y_BC(0,i),&
                               wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + np),&
                               wing%u_old(nint(wing%veins_bending_BC(1,3,i)) + np), &
                               wing%z_BC(0,i),wing%u_old(nint(wing%veins_bending_BC(1,2,i)) + 2*np), &
                               wing%u_old(nint(wing%veins_bending_BC(1,3,i)) + 2*np), &
                               wing%veins_bending_BC(0,5,i),wing%veins_bending_BC(0,6,i))


      enddo

end subroutine
