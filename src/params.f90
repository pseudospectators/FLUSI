! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile,Insect,verbose)
  use vars
  use ini_files_parser_mpi
  use insect_module
  ! The file we read the PARAMS from
  character(len=strlen),intent(in) :: paramsfile
  ! print a copy of the parameters read or not?
  logical, intent(in) :: verbose
  ! the insect we initialize here
  type(diptera), intent(inout) :: Insect
  ! this array contains the entire ascii-params file
  type(inifile) :: PARAMS

  ! Read the paramsfile to the derived dataytpe for ini-files
  call read_ini_file_mpi(PARAMS,paramsfile,verbose)

  ! Get parameter values from PARAMS
  call get_params_common(PARAMS)

  ! read all the specific parameters from the file
  select case(method)
  case("fsi")
    ! Get fsi-specific parameter values from PARAMS
    call get_params_fsi(PARAMS,Insect)
  case("mhd")
    ! Get mhd-specific parameter values from PARAMS
    call get_params_mhd(PARAMS)
  case default
    if(mpirank==0) call abort("Error! Unkonwn method in get_params; stopping.")
  end select


  if (mpirank==0 .and. verbose) then
    write (*,*) "*************************************************"
    write (*,'(A,i3)') " *** DONE READING PARAMETERS"
    write (*,*) "*************************************************"
  endif

  ! clean ini file
  call clean_ini_file_mpi(PARAMS)

end subroutine get_params




! Read individual parameter values from the PARAMS string for the vars
! module.
subroutine get_params_common(PARAMS)
  use ini_files_parser_mpi
  use vars
  implicit none
  character(len=strlen) :: dummystr

  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS
  character(len=strlen) :: dummy

  ! Resolution section
  call read_param_mpi(PARAMS,"Resolution","nx",nx, 4)
  call read_param_mpi(PARAMS,"Resolution","ny",ny, 4)
  call read_param_mpi(PARAMS,"Resolution","nz",nz, 4)

  ! Geometry section
  call read_param_mpi(PARAMS,"Geometry","xl",xl, 6.283185307179586d0)
  call read_param_mpi(PARAMS,"Geometry","yl",yl, 6.283185307179586d0)
  call read_param_mpi(PARAMS,"Geometry","zl",zl, 6.283185307179586d0)

  ! lattice spacing is global (since we allow to specify reals in multiples of
  ! grid points, we nedd that value now.)
  dx=xl/dble(nx)
  dy=yl/dble(ny)
  dz=zl/dble(nz)

  if (nx==1) then
    if (root) write(*,*) "2D run: setting x coordinate accordingly (OVERWRITE!!!)"
    dx = 1.d0!max(dz,dy)
    xl = 1.d0!dx
    if (root) write(*,'("xl=",es12.4," dx=",es12.4)') xl,dx
  endif

  ! Geometry section
  call read_param_mpi(PARAMS,"Geometry","Size",length, 0.d0)
  call read_param_mpi(PARAMS,"Geometry","alpha",alpha_generic, 0.d0)
  call read_param_mpi(PARAMS,"Geometry","r1",r1,1.d0)
  call read_param_mpi(PARAMS,"Geometry","r2",r2,1.0681415d0)
  call read_param_mpi(PARAMS,"Geometry","r3",r3,1.206371d0)

  ! Time section
  call read_param_mpi(PARAMS,"Time","nt",nt, 9999999)
  call read_param_mpi(PARAMS,"Time","iTimeMethodFluid",iTimeMethodFluid,"AB2")
  call read_param_mpi(PARAMS,"Time","intelligent_dt",intelligent_dt,"no")
  call read_param_mpi(PARAMS,"Time","Tmax",Tmax,1.d9)
  call read_param_mpi(PARAMS,"Time","CFL",cfl,0.1d0)
  call read_param_mpi(PARAMS,"Time","dt_max",dt_max,0.d0)
  call read_param_mpi(PARAMS,"Time","dt_fixed",dt_fixed,0.d0)

  ! Reynolds number section:
  call read_param_mpi(PARAMS,"ReynoldsNumber","nu",nu,1.d-2)

  ! Initial conditions section
  call read_param_mpi(PARAMS,"InitialCondition","inicond",inicond, "none")
  call read_param_mpi(PARAMS,"InitialCondition","omega1",omega1,0.d0)
  call read_param_mpi(PARAMS,"InitialCondition","nu_smoothing",nu_smoothing,1.0d-13)
  ! if reading from file, which files?
  call read_param_mpi(PARAMS,"InitialCondition","file_ux",file_ux, "none")
  call read_param_mpi(PARAMS,"InitialCondition","file_uy",file_uy, "none")
  call read_param_mpi(PARAMS,"InitialCondition","file_uz",file_uz, "none")
  ! if running in MHD mode, we also need the B-field initialized
  call read_param_mpi(PARAMS,"InitialCondition","file_bx",file_bx, "none")
  call read_param_mpi(PARAMS,"InitialCondition","file_by",file_by, "none")
  call read_param_mpi(PARAMS,"InitialCondition","file_bz",file_bz, "none")

  ! Dealasing section
  call read_param_mpi(PARAMS,"Dealiasing","iDealias",iDealias, 1)

  ! Penalization section
  call read_param_mpi(PARAMS,"Penalization","iPenalization",iPenalization, 0)
  call read_param_mpi(PARAMS,"Penalization","iMoving",iMoving, 0)
  call read_param_mpi(PARAMS,"Penalization","iMask",iMask, "none")
  call read_param_mpi(PARAMS,"Penalization","eps",eps, 1.d-2)
  call read_param_mpi(PARAMS,"Penalization","pseudoeps",pseudoeps, 1.d-2)
  call read_param_mpi(PARAMS,"Penalization","pseudodt",pseudodt, 1.d-2)
  call read_param_mpi(PARAMS,"Penalization","pseuderrmin",pseudoerrmin,3d-4)
  call read_param_mpi(PARAMS,"Penalization","pseuderrmax",pseudoerrmax,5d-4)
  call read_param_mpi(PARAMS,"Penalization","periodic",dummystr,"no")
  if (dummystr=="yes") then
    periodic=.true.
  endif

  ! turbulent inlet
  call read_param_mpi(PARAMS,"TurbulentInlet","use_turbulent_inlet",use_turbulent_inlet, "no")
  call read_param_mpi(PARAMS,"TurbulentInlet","rescale",rescale, 1.d0)
  call read_param_mpi(PARAMS,"TurbulentInlet","inlet_thickness",inlet_thickness,48)

  ! averaging in time
  call read_param_mpi(PARAMS,"Averaging","time_avg",time_avg, "no")
  call read_param_mpi(PARAMS,"Averaging","vel_avg",vel_avg, "yes")
  call read_param_mpi(PARAMS,"Averaging","ekin_avg",ekin_avg, "no")
  call read_param_mpi(PARAMS,"Averaging","enstrophy_avg",enstrophy_avg, "no")
  call read_param_mpi(PARAMS,"Averaging","save_one_only",save_one_only, "yes")
  call read_param_mpi(PARAMS,"Averaging","tstart_avg",tstart_avg,0.d0)
  if (vel_avg/="yes" .and. ekin_avg/="yes") then
    if (mpirank==0) write(*,*) "You specified to avg, but you want neither u_avg nor ekin_avg"
    if (mpirank==0) write(*,*) "Thus no averaging is done"
    time_avg="no"
  endif
  ! Saving section
  call read_param_mpi(PARAMS,"Saving","iDoBackup",iDoBackup, 1)
  call read_param_mpi(PARAMS,"Saving","iSaveVelocity",iSaveVelocity, 0)
  call read_param_mpi(PARAMS,"Saving","iSavePress",iSavePress, 0)
  call read_param_mpi(PARAMS,"Saving","iSaveVorticity",iSaveVorticity, 0)
  call read_param_mpi(PARAMS,"Saving","iSaveMagVorticity",iSaveMagVorticity, 0)
  call read_param_mpi(PARAMS,"Saving","iSaveMask",iSaveMask, 0)
  call read_param_mpi(PARAMS,"Saving","tsave",tsave, 9.d9)
  call read_param_mpi(PARAMS,"Saving","truntime",truntime, 1.d0)
  call read_param_mpi(PARAMS,"Saving","wtimemax",wtimemax, 8760.d0) ! 1 year
  call read_param_mpi(PARAMS,"Saving","tintegral",tintegral,0.01d0)
  call read_param_mpi(PARAMS,"Saving","tsave_first",tsave_first,0.0d0)
  call read_param_mpi(PARAMS,"Saving","tsave_period",tsave_period,1.0d0)
  call read_param_mpi(PARAMS,"Saving","save_only_one_period",&
  save_only_one_period,"no")
  call read_param_mpi(PARAMS,"Saving","naming",naming,"time")
  call read_param_mpi(PARAMS,"Saving","field_precision",&
  field_precision,"single")
  call read_param_mpi(PARAMS,"Saving","itdrag",itdrag,99999)
  call read_param_mpi(PARAMS,"Saving","itbeam",itbeam,99999)
  call read_param_mpi(PARAMS,"Saving","iSaveSpectrae",iSaveSpectrae,"no")
  call read_param_mpi(PARAMS,"Saving","striding",striding,1)

  ! Forcing section (for isotropic turbulence)
  call read_param_mpi(PARAMS,"Forcing","forcing_type",forcing_type, "none")
  call read_param_mpi(PARAMS,"Forcing","kf",kf, 0.d0)
  call read_param_mpi(PARAMS,"Forcing","eps_forcing",eps_forcing, 0.d0)

  !-- dry run, just the mask function
  call read_param_mpi(PARAMS,"DryRun","dry_run_without_fluid",dummy,"no")
  if (dummy=="yes") then
    call abort('dry_run_without_fluid is depecreated; run flusi with ./flusi --dry-run PARAMS.ini instead')
  endif


  ! Set other parameters (all procs)
  pi=4.d0 *datan(1.d0)
  ! scaling for FFTs
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
end subroutine get_params_common


! Read individual parameter values from the PARAMS string for fsi.
subroutine get_params_fsi(PARAMS,Insect)
  use ini_files_parser_mpi
  use vars
  use insect_module
  implicit none

  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS
  ! the insect to initialize
  type(diptera), intent(inout) :: Insect
  real(kind=pr), dimension(1:3) :: defaultvec
  character(len=strlen) :: old_meanflow
  ! ---------------------------------------------------
  ! penalization / cavity
  ! ---------------------------------------------------
  call read_param_mpi(PARAMS,"Penalization","iCavity",iCavity,"no")
  call read_param_mpi(PARAMS,"Penalization","cavity_size",cavity_size,0)
  call read_param_mpi(PARAMS,"Penalization","compute_forces",compute_forces,1)
  call read_param_mpi(PARAMS,"Penalization","unst_corrections",unst_corrections,0)
  call read_param_mpi(PARAMS,"Penalization","iChannel",iChannel,"no")
  if (iChannel=="0") iChannel="no" ! for downward compatibility with older ini files
  if (iChannel=="1") iChannel="xy" ! for downward compatibility with older ini files
  call read_param_mpi(PARAMS,"Penalization","thick_wall",thick_wall,0.2d0)
  call read_param_mpi(PARAMS,"Penalization","pos_wall",pos_wall,0.3d0)
  call read_param_mpi(PARAMS,"Penalization","us_fixed",us_fixed,(/0.d0,0.d0,0.d0/))
  ! ---------------------------------------------------
  ! sponge
  ! ---------------------------------------------------
  call read_param_mpi(PARAMS,"Sponge","iVorticitySponge",iVorticitySponge,"no")
  call read_param_mpi(PARAMS,"Sponge","iSpongeType",iSpongeType,"top_cover")
  call read_param_mpi(PARAMS,"Sponge","eps_sponge",eps_sponge, 0.d0)
  call read_param_mpi(PARAMS,"Sponge","sponge_thickness",sponge_thickness,0)


  ! ---------------------------------------------------
  ! Geometry section
  ! ---------------------------------------------------
  call read_param_mpi(PARAMS,"Geometry","x0",x0, xl/2.d0)
  call read_param_mpi(PARAMS,"Geometry","y0",y0, yl/2.d0)
  call read_param_mpi(PARAMS,"Geometry","z0",z0, zl/2.d0)

  ! ---------------------------------------------------
  ! Saving section
  ! ---------------------------------------------------
  call read_param_mpi(PARAMS,"Saving","iSaveSolidVelocity",iSaveSolidVelocity,0)

  ! ---------------------------------------------------
  ! MeanFlow section
  ! ---------------------------------------------------
  call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_x",iMeanFlow_x,"free")
  call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_y",iMeanFlow_y,"free")
  call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow_z",iMeanFlow_z,"free")
  ! for compatibility with old params files:
  call read_param_mpi(PARAMS,"MeanFlow","iMeanFlow",old_meanflow,"unused")
  call read_param_mpi(PARAMS,"MeanFlow","m_fluid",m_fluid, 999.d9)
  call read_param_mpi(PARAMS,"MeanFlow","ux",uxmean, 1.d0)
  call read_param_mpi(PARAMS,"MeanFlow","uy",uymean, 1.d0)
  call read_param_mpi(PARAMS,"MeanFlow","uz",uzmean, 1.d0)
  call read_param_mpi(PARAMS,"MeanFlow","iMeanFlowStartupConditioner",&
  iMeanFlowStartupConditioner,"no")
  call read_param_mpi(PARAMS,"MeanFlow","tau_meanflow",tau_meanflow, 0.d0)
  call read_param_mpi(PARAMS,"MeanFlow","T_release_meanflow",T_release_meanflow,0.d0)

  ! for compatibility with old files:
  if (old_meanflow=="1") then
    iMeanFlow_x="fixed"
    iMeanFlow_y="fixed"
    iMeanFlow_z="fixed"
  endif

  ! ---------------------------------------------------
  ! Insects section
  ! ---------------------------------------------------
  if (iMask=="Insect") then
    call get_params_insect( PARAMS,Insect )
  endif

  ! ---------------------------------------------------
  ! solid model (TODO: SAME LEVEL OF OBJECT ORIENTATION AS INSECT)
  ! ---------------------------------------------------
  call get_params_solid( PARAMS )
  call get_params_scalars( PARAMS )

  ! ----------------------------------------------------
  ! slice extraction
  ! ----------------------------------------------------
  ! the first one is for compatibility; it is overwritten if the 2nd is found
  call read_param_mpi(PARAMS,"SaveSlices","use_slicing",use_slicing,"xxx")
  if (use_slicing=="xxx") then
    call read_param_mpi(PARAMS,"SaveSlices","save_slices",use_slicing,"no")
  endif
  call read_param_mpi(PARAMS,"SaveSlices","slice1",slices_to_save(1),-2)
  call read_param_mpi(PARAMS,"SaveSlices","slice2",slices_to_save(2),-2)
  call read_param_mpi(PARAMS,"SaveSlices","slice3",slices_to_save(3),-2)
  call read_param_mpi(PARAMS,"SaveSlices","slice4",slices_to_save(4),-2)
  call read_param_mpi(PARAMS,"SaveSlices","itslice",itslice,9999900)
  call read_param_mpi(PARAMS,"SaveSlices","ncache_slices",ncache_slices,nx)
  call read_param_mpi(PARAMS,"SaveSlices","tslice",tslice,99999.9d0)
  call read_param_mpi(PARAMS,"SaveSlices","tslice_first",tslice_first,0.0d0)
  ! ---------------------------------------------------
  ! DONE..
  ! ---------------------------------------------------

  lin(1)=nu
end subroutine get_params_fsi


!-------------------------------------------------------------------------------
! This routine reads in the parameters that describe the inscet from the
! parameter.ini file. it is outsourced from params.f90
!-------------------------------------------------------------------------------
subroutine get_params_insect( PARAMS,Insect )
  use vars
  use insect_module
  use ini_files_parser_mpi
  implicit none

  type(diptera),intent(inout) :: Insect
  integer :: j, tmp
  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS
  real(kind=pr),dimension(1:3)::defaultvec
  character(len=strlen) :: DoF_string, dummystr

  call read_param_mpi(PARAMS,"Insects","WingShape",Insect%WingShape,"none")
  call read_param_mpi(PARAMS,"Insects","b_top",Insect%b_top, 0.d0)
  call read_param_mpi(PARAMS,"Insects","b_bot",Insect%b_bot, 0.d0)
  call read_param_mpi(PARAMS,"Insects","L_chord",Insect%L_chord, 0.d0)
  call read_param_mpi(PARAMS,"Insects","L_span",Insect%L_span, 0.d0)
  call read_param_mpi(PARAMS,"Insects","FlappingMotion_right",Insect%FlappingMotion_right,"none")
  call read_param_mpi(PARAMS,"Insects","FlappingMotion_left",Insect%FlappingMotion_left,"none")
  ! this file is used in old syntax form for both wings:
  call read_param_mpi(PARAMS,"Insects","infile",Insect%infile,"none.in")

  !-----------------------------------------------------------------------------
  ! compatibility for wingbeat kinematics read from file
  !-----------------------------------------------------------------------------

  if ( index(Insect%FlappingMotion_right,"from_file::") /= 0 ) then
    ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
    Insect%kine_wing_r%infile = Insect%FlappingMotion_right( 12:strlen  )
    Insect%FlappingMotion_right = "from_file"

  elseif ( index(Insect%FlappingMotion_right,"kinematics_loader::") /= 0 ) then
    ! new syntax, uses the kinematics loader for non-periodic kinematics
    Insect%kine_wing_r%infile = Insect%FlappingMotion_right( 20:strlen )
    Insect%FlappingMotion_right = "kinematics_loader"

  elseif ( Insect%FlappingMotion_right == "from_file" ) then
    ! old syntax, implies symmetric periodic motion, read from *.ini file
    Insect%kine_wing_r%infile = Insect%infile

  elseif ( Insect%FlappingMotion_right == "kinematics_loader" ) then
    ! old syntax, implies symmetric non-periodic motion, read from *.dat file
    Insect%kine_wing_r%infile = Insect%infile
  endif


  if ( index(Insect%FlappingMotion_left,"from_file::") /= 0 ) then
    ! new syntax, uses fourier/hermite periodic kinematics read from *.ini file
    Insect%kine_wing_l%infile = Insect%FlappingMotion_left( 12:strlen  )
    Insect%FlappingMotion_left = "from_file"

  elseif ( index(Insect%FlappingMotion_left,"kinematics_loader::") /= 0 ) then
    ! new syntax, uses the kinematics loader for non-periodic kinematics
    Insect%kine_wing_l%infile = Insect%FlappingMotion_left( 20:strlen )
    Insect%FlappingMotion_left = "kinematics_loader"

  elseif ( Insect%FlappingMotion_left == "from_file" ) then
    ! old syntax, implies symmetric periodic motion, read from *.ini file
    Insect%kine_wing_l%infile = Insect%infile

  elseif ( Insect%FlappingMotion_left == "kinematics_loader" ) then
    ! old syntax, implies symmetric non-periodic motion, read from *.dat file
    Insect%kine_wing_l%infile = Insect%infile
  endif

  if (root) then
    write(*,*) "Left wing: "//trim(adjustl(Insect%FlappingMotion_left))
    write(*,*) "Left wing: "//trim(adjustl(Insect%kine_wing_l%infile))
    write(*,*) "Right wing: "//trim(adjustl(Insect%FlappingMotion_right))
    write(*,*) "Right wing: "//trim(adjustl(Insect%kine_wing_r%infile))
  endif

  ! these flags trigger reading the kinematics from file when the Flapping
  ! motion is first called
  Insect%kine_wing_l%initialized = .false.
  Insect%kine_wing_r%initialized = .false.

  call read_param_mpi(PARAMS,"Insects","BodyType",Insect%BodyType,"ellipsoid")
  call read_param_mpi(PARAMS,"Insects","HasDetails",Insect%HasDetails,"all")
  call read_param_mpi(PARAMS,"Insects","BodyMotion",Insect%BodyMotion,"yes")
  call read_param_mpi(PARAMS,"Insects","LeftWing",Insect%LeftWing,"yes")
  call read_param_mpi(PARAMS,"Insects","RightWing",Insect%RightWing,"yes")
  call read_param_mpi(PARAMS,"Insects","b_body",Insect%b_body, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","L_body",Insect%L_body, 1.d0)
  call read_param_mpi(PARAMS,"Insects","R_head",Insect%R_head, 0.1d0)
  call read_param_mpi(PARAMS,"Insects","R_eye",Insect%R_eye, 0.d1)
  call read_param_mpi(PARAMS,"Insects","mass",Insect%mass, 1.d0)
  call read_param_mpi(PARAMS,"Insects","gravity",Insect%gravity, 1.d0)
  call read_param_mpi(PARAMS,"Insects","WingThickness",Insect%WingThickness, 4.0d0*dx)
  call read_param_mpi(PARAMS,"Insects","J_body_yawpitchroll",defaultvec, (/0.d0,0.d0,0.d0/))
  Insect%Jroll_body  = defaultvec(3)
  Insect%Jyaw_body   = defaultvec(1)
  Insect%Jpitch_body = defaultvec(2)
  call read_param_mpi(PARAMS,"Insects","x0",Insect%x0, (/0.5d0*xl,0.5d0*yl,0.5d0*zl/))
  call read_param_mpi(PARAMS,"Insects","v0",Insect%v0, (/0.d0, 0.d0, 0.d0/))
  call read_param_mpi(PARAMS,"Insects","yawpitchroll_0",Insect%yawpitchroll_0,&
  (/0.d0, 0.d0, 0.d0/))
  ! convert yawpitchroll to radiants
  Insect%yawpitchroll_0 = Insect%yawpitchroll_0 * (pi/180.d0)
  call read_param_mpi(PARAMS,"Insects","eta0",Insect%eta0, 0.0d0)
  Insect%eta0 = Insect%eta0*(pi/180.d0)


  ! degrees of freedom for free flight solver. The string from ini file contains
  ! 6 characters 1 or 0 that turn on/off x,y,z,yaw,pitch,roll degrees of freedom
  ! by multiplying the respective RHS by zero, keeping the value thus constant
  call read_param_mpi(PARAMS,"Insects","DoF",DoF_string, "111111")
  do j=1,6
    read (DoF_string(j:j), '(i1)') tmp
    Insect%DoF_on_off(j) = dble(tmp)
  enddo
  if (root) write(*,'(6(f4.2,1x))') Insect%DoF_on_off


  ! wing inertia tensor (we currently assume two identical wings)
  ! this allows computing inertial power
  call read_param_mpi(PARAMS,"Insects","Jxx",Insect%Jxx,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jyy",Insect%Jyy,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jzz",Insect%Jzz,0.d0)
  call read_param_mpi(PARAMS,"Insects","Jxy",Insect%Jxy,0.d0)

  call read_param_mpi(PARAMS,"Insects","startup_conditioner",Insect%startup_conditioner,"no")

  ! position vector of the head
  call read_param_mpi(PARAMS,"Insects","x_head",&
  Insect%x_head, (/0.5d0*Insect%L_body,0.d0,0.d0 /) )

  ! eyes
  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.0d0,+1.0d0,1.0d0/)
  call read_param_mpi(PARAMS,"Insects","x_eye_r",Insect%x_eye_r, defaultvec)

  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.0d0,-1.0d0,1.0d0/)
  call read_param_mpi(PARAMS,"Insects","x_eye_l",Insect%x_eye_l, defaultvec)

  ! wing hinges (root points)
  defaultvec=(/0.d0, +Insect%b_body, 0.d0 /)
  call read_param_mpi(PARAMS,"Insects","x_pivot_l",Insect%x_pivot_l, defaultvec)

  defaultvec=(/0.d0, -Insect%b_body, 0.d0 /)
  call read_param_mpi(PARAMS,"Insects","x_pivot_r",Insect%x_pivot_r, defaultvec)

  Insect%smooth = 2.0d0*dz
end subroutine get_params_insect




!-------------------------------------------------------------------------------
! Read individual parameter values that are specific to the solid model only
!-------------------------------------------------------------------------------
subroutine get_params_solid(PARAMS)
  use ini_files_parser_mpi
  use solid_model
  implicit none

  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS

  !-- solid model is deactivated by default
  call read_param_mpi(PARAMS,"SolidModel","use_solid_model",use_solid_model,"no")

  !-- if using the solid model, look for other parameters
  if (use_solid_model=="yes") then
    !-- beam resolution
    call read_param_mpi(PARAMS,"SolidModel","ns",ns, 32)
    !-- interpolation method
    call read_param_mpi(PARAMS,"SolidModel","interp",interp,"delta")
    !-- density / stiffness / gravity
    call read_param_mpi(PARAMS,"SolidModel","mue",mue0,1.0d0)
    call read_param_mpi(PARAMS,"SolidModel","eta",eta0,1.0d0)
    call read_param_mpi(PARAMS,"SolidModel","f",frequ,1.0d0)
    call read_param_mpi(PARAMS,"SolidModel","angle",AngleBeam,1.0d0)
    call read_param_mpi(PARAMS,"SolidModel","gravity",grav,0.0d0)
    !-- beam thickness
    call read_param_mpi(PARAMS,"SolidModel","t_beam",t_beam,0.05d0)
    !-- damping
    call read_param_mpi(PARAMS,"SolidModel","sigma",sigma,0.0d0)
    !-- timing
    call read_param_mpi(PARAMS,"SolidModel","T_release",T_release,0.0d0)
    call read_param_mpi(PARAMS,"SolidModel","tau",tau,0.0d0)
    call read_param_mpi(PARAMS,"SolidModel","N_smooth",N_smooth,3.0d0)
    call read_param_mpi(PARAMS,"SolidModel","L_span",L_span,1.0d0)
    call read_param_mpi(PARAMS,"SolidModel","has_cylinder",has_cylinder,"no")
    call read_param_mpi(PARAMS,"SolidModel","R_cylinder",R_cylinder,0.0d0)
    !-- time marching method for the solid
    call read_param_mpi(PARAMS,"SolidModel","TimeMethodSolid",TimeMethodSolid,"BDF2")

    select case (TimeMethodSolid)
    case ("RK4","CN2","BDF2","EI1","EE1","prescribed")
      if (mpirank==0) write(*,*) "Solid solver is ", TimeMethodSolid
    case default
      if (mpirank==0) write(*,*) "Solid solver is UNDEFINED, using BDF2"
      TimeMethodSolid="BDF2"
    end select

    call read_param_mpi(PARAMS,"SolidModel","imposed_motion_leadingedge",&
    imposed_motion_leadingedge,"fixed_middle")
    call read_param_mpi(PARAMS,"SolidModel","infinite",infinite,"no")
    call read_param_mpi(PARAMS,"SolidModel","plate_shape",plate_shape,"rectangular")
    call read_param_mpi(PARAMS,"SolidModel","debug_pressure",debug_pressure,0)

    !-- grid spacing
    ds = 1.d0/dble(ns-1)
  endif

end subroutine get_params_solid


!-------------------------------------------------------------------------------
! Read individual parameter values that are specific to passive scalars
!-------------------------------------------------------------------------------
subroutine get_params_scalars(PARAMS)
  use ini_files_parser_mpi
  use solid_model
  use passive_scalar_module
  implicit none

  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS
  character(len=7) :: name
  integer :: j

  call read_param_mpi(PARAMS,"PassiveScalar","use_passive_scalar",use_passive_scalar,0)
  call read_param_mpi(PARAMS,"PassiveScalar","n_scalars",n_scalars,0)
  call read_param_mpi(PARAMS,"PassiveScalar","stop_on_fail",stop_on_fail,"no")

  !-- if using passive scalars, read their individual parameters
  if (use_passive_scalar==1) then
    allocate(scalar_props(1:n_scalars))
    ! loop over scalars
    do j=1,n_scalars
      ! we now read sections "Scalar1" "Scalar2" and so on
      write (name,'("Scalar",i1)') j
      if (mpirank==0) write(*,*) "reading scalar "//name

      call read_param_mpi(PARAMS,name,"kappa",scalar_props(j)%kappa,0.1d0)
      ! defalut is same penalization as fluid
      call read_param_mpi(PARAMS,name,"eps",scalar_props(j)%eps, eps)
      call read_param_mpi(PARAMS,name,"inicond",scalar_props(j)%inicond,"right_left_discontinuous")
      call read_param_mpi(PARAMS,name,"sourceterm",scalar_props(j)%sourceterm,"none")
      call read_param_mpi(PARAMS,name,"sourceterm_x0",scalar_props(j)%source_x0,(/0.d0,0.d0,0.d0,0.d0/))
    enddo
  endif

end subroutine get_params_scalars



! Read individual parameter values from the PARAMS string for fmhd
subroutine get_params_mhd(PARAMS)
  use ini_files_parser_mpi
  use vars
  implicit none

  ! Contains the ascii-params file
  type(inifile), intent(inout) :: PARAMS

  ! MHD section
  call read_param_mpi(PARAMS,"MHD","eta",eta,4.5d-2)

  ! MHDGeometry section
  call read_param_mpi(PARAMS,"MHDGeometry","b0",b0,4.5d0)
  call read_param_mpi(PARAMS,"MHDGeometry","bc",bc,3.88888888888d0)
  call read_param_mpi(PARAMS,"MHDGeometry","ay",ay,1.0d0)

  ! Saving section
  call read_param_mpi(PARAMS,"Saving","iSaveMagneticField",&
  iSaveMagneticField, 0)
  call read_param_mpi(PARAMS,"Saving","iSaveCurrent",iSaveCurrent, 0)

  lin(1)=nu
  lin(2)=eta
end subroutine get_params_mhd
