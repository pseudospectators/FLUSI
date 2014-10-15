! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile,Insect) 
  use vars
  use insect_module
  ! The file we read the PARAMS from
  character(len=strlen),intent(in) :: paramsfile
  ! the insect we initialize here
  type(diptera), intent(inout) :: Insect
  integer :: i  
  ! this array contains the entire ascii-params file
  character(len=strlen), dimension(1:nlines) :: PARAMS
  real(kind=pr), dimension(1:3) :: defaultvec
  character(len=strlen) :: old_meanflow

  
  ! check if the specified file exists
  call check_file_exists( paramsfile )
  
  ! Read the paramsfile and put the length i and the text in PARAMS
  call read_params_file(PARAMS,i,paramsfile,.true.)


  ! Resolution section
  call param_int(PARAMS,i,"Resolution","nx",nx, 4)
  call param_int(PARAMS,i,"Resolution","ny",ny, 4)
  call param_int(PARAMS,i,"Resolution","nz",nz, 4)
    
  ! Geometry section
  call param_dbl(PARAMS,i,"Geometry","xl",xl, 1.d0)
  call param_dbl(PARAMS,i,"Geometry","yl",yl, 1.d0)
  call param_dbl(PARAMS,i,"Geometry","zl",zl, 1.d0)
      
  ! lattice spacing is global (since we allow to specify reals in multiples of
  ! grid points, we nedd that value now.)
  dx=xl/dble(nx)
  dy=yl/dble(ny)
  dz=zl/dble(nz) 
  pi=4.d0 *datan(1.d0)
  ! scaling for FFTs
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  
  if (nx==1) then
    if (root) write(*,*) "2D run: setting x coordinate accordingly (OVERWRITE!!!)"    
    dx = 1.d0
    xl = 1.d0
    if (root) write(*,'("xl=",es12.4," dx=",es12.4)') xl,dx
  endif
    
  ! Geometry section
  call param_dbl(PARAMS,i,"Geometry","Size",length, 0.d0)

  !-----------------------------------------------------------------------------
  ! Time section
  !-----------------------------------------------------------------------------
  call param_int(PARAMS,i,"Time","nt",nt, 9999999)
  call param_str(PARAMS,i,"Time","iTimeMethodFluid",iTimeMethodFluid,"AB2")  
  call param_dbl(PARAMS,i,"Time","Tmax",Tmax,1.d9)
  call param_dbl(PARAMS,i,"Time","CFL",cfl,0.1d0)
  call param_dbl(PARAMS,i,"Time","dt_max",dt_max,0.d0)
  call param_dbl(PARAMS,i,"Time","dt_fixed",dt_fixed,0.d0)

  !-----------------------------------------------------------------------------
  ! Reynolds number section:
  !-----------------------------------------------------------------------------
  call param_dbl(PARAMS,i,"ReynoldsNumber","nu",nu,1.d-2)  

  !-----------------------------------------------------------------------------
  ! Initial conditions section
  !-----------------------------------------------------------------------------
  call param_str(PARAMS,i,"InitialCondition","inicond",inicond, "none")
  ! if inicond reading from file, which files?
  if (inicond=="infile") then 
    call param_str(PARAMS,i,"InitialCondition","file_ux",file_ux, "none")
    call param_str(PARAMS,i,"InitialCondition","file_uy",file_uy, "none")
    call param_str(PARAMS,i,"InitialCondition","file_uz",file_uz, "none")
    call param_str(PARAMS,i,"InitialCondition","file_p",file_p, "none")
  endif

  !-----------------------------------------------------------------------------
  ! Saving section
  !-----------------------------------------------------------------------------
  call param_int(PARAMS,i,"Saving","iDoBackup",iDoBackup, 1)
  call param_int(PARAMS,i,"Saving","iSaveVelocity",iSaveVelocity, 0) 
  call param_int(PARAMS,i,"Saving","iSavePress",iSavePress, 0)
  call param_int(PARAMS,i,"Saving","iSaveVorticity",iSaveVorticity, 0)  
  call param_int(PARAMS,i,"Saving","iSaveMask",iSaveMask, 0)
  call param_int(PARAMS,i,"Saving","iSaveXMF",iSaveXMF, 0) ! default is no
  call param_dbl(PARAMS,i,"Saving","tsave",tsave, 9.d9)
  call param_dbl(PARAMS,i,"Saving","truntime",truntime, 1.d0)
  call param_dbl(PARAMS,i,"Saving","wtimemax",wtimemax, 8760.d0) ! 1 year
  call param_dbl(PARAMS,i,"Saving","tintegral",tintegral,0.01d0)
  call param_dbl(PARAMS,i,"Saving","tsave_first",tsave_first,0.0d0)
  call param_dbl(PARAMS,i,"Saving","tsave_period",tsave_period,1.0d0)
  call param_str(PARAMS,i,"Saving","save_only_one_period",&
       save_only_one_period,"no")  
  call param_int(PARAMS,i,"Saving","itdrag",itdrag,99999)
  call param_int(PARAMS,i,"Saving","itbeam",itbeam,99999)
  
  !-- dry run, just the mask function
  call param_str(PARAMS,i,"DryRun","dry_run_without_fluid",dry_run_without_fluid,"no")
  if (dry_run_without_fluid=="yes") then
    write(*,*) "Attention! This is a dry run without fluid"
    write(*,*) "Deactivating all useless save-switches..."
    idobackup=0
    iSavePress=0
    iSaveVelocity=0
    iSaveVorticity=0
  endif
  
  !-----------------------------------------------------------------------------
  ! penalization / cavity
  !-----------------------------------------------------------------------------
  call param_int(PARAMS,i,"Penalization","iPenalization",iPenalization, 0)
  call param_int(PARAMS,i,"Penalization","iMoving",iMoving, 0)
  call param_str(PARAMS,i,"Penalization","iMask",iMask, "none")
  call param_dbl(PARAMS,i,"Penalization","eps",eps, 1.d-2)
  call param_str(PARAMS,i,"Penalization","iCavity",iCavity,"no") 
  call param_int(PARAMS,i,"Penalization","cavity_size",cavity_size,0)
  call param_int(PARAMS,i,"Penalization","compute_forces",compute_forces,1)   
  call param_int(PARAMS,i,"Penalization","unst_corrections",unst_corrections,0)        
  call param_str(PARAMS,i,"Penalization","iChannel",iChannel,"no") 
  if (iChannel=="0") iChannel="no" ! for downward compatibility with older ini files
  if (iChannel=="1") iChannel="xy" ! for downward compatibility with older ini files
  call param_dbl(PARAMS,i,"Penalization","thick_wall",thick_wall,0.2d0)
  call param_dbl(PARAMS,i,"Penalization","pos_wall",pos_wall,0.3d0)
  
  !-----------------------------------------------------------------------------
  ! Geometry section
  !-----------------------------------------------------------------------------
  call param_dbl(PARAMS,i,"Geometry","x0",x0, 0.d0)
  call param_dbl(PARAMS,i,"Geometry","y0",y0, 0.d0)
  call param_dbl(PARAMS,i,"Geometry","z0",z0, 0.d0)
  
  !-----------------------------------------------------------------------------
  ! Saving section
  !-----------------------------------------------------------------------------
  call param_int(PARAMS,i,"Saving","iSaveSolidVelocity",iSaveSolidVelocity,0)

  !-----------------------------------------------------------------------------
  ! MeanFlow section
  !-----------------------------------------------------------------------------
  call param_str(PARAMS,i,"MeanFlow","iMeanFlow_x",iMeanFlow_x,"free")
  call param_str(PARAMS,i,"MeanFlow","iMeanFlow_y",iMeanFlow_y,"free")
  call param_str(PARAMS,i,"MeanFlow","iMeanFlow_z",iMeanFlow_z,"free")
  ! for compatibility with old params files:
  call param_str(PARAMS,i,"MeanFlow","iMeanFlow",old_meanflow,"unused")
  call param_dbl(PARAMS,i,"MeanFlow","m_fluid",m_fluid, 999.d9) 
  call param_dbl(PARAMS,i,"MeanFlow","ux",uxmean, 1.d0) 
  call param_dbl(PARAMS,i,"MeanFlow","uy",uymean, 1.d0) 
  call param_dbl(PARAMS,i,"MeanFlow","uz",uzmean, 1.d0) 
  call param_str(PARAMS,i,"MeanFlow","iMeanFlowStartupConditioner",&
       iMeanFlowStartupConditioner,"no")
  call param_dbl(PARAMS,i,"MeanFlow","tau_meanflow",tau_meanflow, 0.d0) 
  call param_dbl(PARAMS,i,"MeanFlow","T_release_meanflow",T_release_meanflow,0.d0) 

  ! for compatibility with old files:
  if (old_meanflow=="1") then 
    iMeanFlow_x="fixed"
    iMeanFlow_y="fixed"
    iMeanFlow_z="fixed"
  endif
    
  !-----------------------------------------------------------------------------
  ! Insects section
  !-----------------------------------------------------------------------------
  if (iMask=="Insect") then
    call read_insect_parameters( PARAMS,i,Insect ) 
  endif

  !-----------------------------------------------------------------------------
  ! Incompressibility
  !-----------------------------------------------------------------------------
  call param_dbl(PARAMS,i,"Incompressibility","c_0",&
       c_0, 0.d0)
  call param_dbl(PARAMS,i,"Incompressibility","gamma_p",&
       gamma_p, 1.d0)     
  call param_str(PARAMS,i,"Incompressibility","method",&
       method,"centered_2nd")
       
  !-----------------------------------------------------------------------------
  ! solid model (TODO: SAME LEVEL OF OBJECT ORIENTATION AS INSECT)
  !-----------------------------------------------------------------------------
  call get_params_solid( PARAMS, i )

  !-----------------------------------------------------------------------------
  ! passive scalar section
  !-----------------------------------------------------------------------------
  call param_int(PARAMS,i,"PassiveScalar","use_passive_scalar",&
       use_passive_scalar,0)
  call param_int(PARAMS,i,"PassiveScalar","n_scalars",&
       n_scalars,0)
  call param_dbl(PARAMS,i,"PassiveScalar","kappa",&
       kappa, 0.1d0)   
  call param_dbl(PARAMS,i,"PassiveScalar","eps_scalar",&
       eps_scalar, eps)       
  call param_str(PARAMS,i,"PassiveScalar","inicond_scalar",&
       inicond_scalar,"right_left_discontinuous") 
  call param_str(PARAMS,i,"PassiveScalar","stop_on_fail",&
       stop_on_fail,"no")
  call param_str(PARAMS,i,"PassiveScalar","source_term",&
       source_term,"no")
  call param_dbl(PARAMS,i,"PassiveScalar","source_xmin",source_xmin,0.d0)
  call param_dbl(PARAMS,i,"PassiveScalar","source_xmax",source_xmax,0.d0)
  call param_dbl(PARAMS,i,"PassiveScalar","source_ymin",source_ymin,0.d0)
  call param_dbl(PARAMS,i,"PassiveScalar","source_ymax",source_ymax,0.d0)
  call param_dbl(PARAMS,i,"PassiveScalar","source_zmin",source_zmin,0.d0)
  call param_dbl(PARAMS,i,"PassiveScalar","source_zmax",source_zmax,0.d0)

  !-----------------------------------------------------------------------------
  ! DONE..
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
     write (*,*) "*************************************************"
     write (*,'(A,i3)') " *** DONE READING PARAMETERS"
     write (*,*) "*************************************************"
  endif
end subroutine 


!-------------------------------------------------------------------------------
! This routine reads in the parameters that describe the inscet from the 
! parameter.ini file. it is outsourced from params.f90
!-------------------------------------------------------------------------------
subroutine read_insect_parameters( PARAMS,i,Insect ) 
  use vars
  use insect_module
  implicit none
  
  type(diptera),intent(inout) :: Insect
  integer,intent(in) :: i
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  real(kind=pr),dimension(1:3)::defaultvec
  
  call param_str(PARAMS,i,"Insects","WingShape",Insect%WingShape,"none")
  call param_dbl(PARAMS,i,"Insects","b_top",Insect%b_top, 0.d0) 
  call param_dbl(PARAMS,i,"Insects","b_bot",Insect%b_bot, 0.d0) 
  call param_dbl(PARAMS,i,"Insects","L_chord",Insect%L_chord, 0.d0) 
  call param_dbl(PARAMS,i,"Insects","L_span",Insect%L_span, 0.d0) 
  call param_str(PARAMS,i,"Insects","FlappingMotion_right",Insect%FlappingMotion_right,"none")
  call param_str(PARAMS,i,"Insects","FlappingMotion_left",Insect%FlappingMotion_left,"none")
  call param_str(PARAMS,i,"Insects","BodyType",Insect%BodyType,"ellipsoid")  
  call param_str(PARAMS,i,"Insects","BodyMotion",Insect%BodyMotion,"yes")
  call param_str(PARAMS,i,"Insects","LeftWing",Insect%LeftWing,"yes")
  call param_str(PARAMS,i,"Insects","RightWing",Insect%RightWing,"yes") 
  call param_dbl(PARAMS,i,"Insects","b_body",Insect%b_body, 0.1d0) 
  call param_dbl(PARAMS,i,"Insects","L_body",Insect%L_body, 1.d0)
  call param_dbl(PARAMS,i,"Insects","R_head",Insect%R_head, 0.1d0) 
  call param_dbl(PARAMS,i,"Insects","R_eye",Insect%R_eye, 0.d1) 
  call param_dbl(PARAMS,i,"Insects","distance_from_sponge",Insect%distance_from_sponge, 1.d0) 
  call param_dbl(PARAMS,i,"Insects","WingThickness",Insect%WingThickness, 4.0*dx) 
  ! wing inertia tensor (we currentlz assume two identical wings)
  ! this allows computing inertial power
  call param_dbl(PARAMS,i,"Insects","Jxx",Insect%Jxx,0.d0)
  call param_dbl(PARAMS,i,"Insects","Jyy",Insect%Jyy,0.d0)
  call param_dbl(PARAMS,i,"Insects","Jzz",Insect%Jzz,0.d0)
  call param_dbl(PARAMS,i,"Insects","Jxy",Insect%Jxy,0.d0)
  call param_str(PARAMS,i,"Insects","infile",Insect%infile,"none.in")
  
  ! position vector of the head
  call param_vct(PARAMS,i,"Insects","x_head",&
       Insect%x_head, (/0.5*Insect%L_body,0.d0,0.d0 /) ) 
  
  ! eyes
  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8*(/1.,+1.,1./)
  call param_vct(PARAMS,i,"Insects","x_eye_r",Insect%x_eye_r, defaultvec) 
       
  defaultvec = Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8*(/1.,-1.,1./)
  call param_vct(PARAMS,i,"Insects","x_eye_l",Insect%x_eye_l, defaultvec) 
       
  ! wing hinges (root points)    
  defaultvec=(/0.d0, +Insect%b_body, 0.d0 /)    
  call param_vct(PARAMS,i,"Insects","x_pivot_l",Insect%x_pivot_l, defaultvec)  
       
  defaultvec=(/0.d0, -Insect%b_body, 0.d0 /)
  call param_vct(PARAMS,i,"Insects","x_pivot_r",Insect%x_pivot_r, defaultvec)        
     
  Insect%smooth = 2.0*dz

  ! flag: read kinematics from file (Dmitry, 14 Nov 2013)
  call param_str(PARAMS,i,"Insects","KineFromFile",Insect%KineFromFile,"no")     
  
 
  ! Takeoff 
  call param_dbl(PARAMS,i,"Insects","x_takeoff",Insect%x_takeoff, 2.0d0)
  call param_dbl(PARAMS,i,"Insects","z_takeoff",Insect%z_takeoff, 0.86d0)
  call param_dbl(PARAMS,i,"Insects","mass_solid",&
       Insect%mass_solid, 54.414118839786745d0)
  call param_dbl(PARAMS,i,"Insects","gravity",&
       Insect%gravity, -0.055129281110537755d0)

  ! Legs model parameters
  call param_int(PARAMS,i,"Insects","ilegs",Insect%ilegs, 1)
  call param_dbl(PARAMS,i,"Insects","anglegsend",&
       Insect%anglegsend, 0.7853981633974483d0)
  call param_dbl(PARAMS,i,"Insects","kzlegsmax",&
       Insect%kzlegsmax,64.24974647375242d0)
  call param_dbl(PARAMS,i,"Insects","dzlegsmax",&
       Insect%dzlegsmax,0.2719665271966527d0)
  call param_dbl(PARAMS,i,"Insects","t0legs",&
       Insect%t0legs,0.13643141797265643d0)
  call param_dbl(PARAMS,i,"Insects","tlinlegs",&
       Insect%tlinlegs,0.3547216867289067d0)
  
end subroutine read_insect_parameters


!-------------------------------------------------------------------------------
! Read individual parameter values that are specific to the solid model only
!-------------------------------------------------------------------------------
subroutine get_params_solid(PARAMS,i)
  use mpi
  use solid_model
  implicit none

  integer,intent(in) :: i
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS

  !-- solid model is deactivated by default
  call param_str(PARAMS,i,"SolidModel","use_solid_model",use_solid_model,"no")
  
  !-- if using the solid model, look for other parameters
  if (use_solid_model=="yes") then
    !-- beam resolution
    call param_int(PARAMS,i,"SolidModel","ns",ns, 32)
    !-- interpolation method
    call param_str(PARAMS,i,"SolidModel","interp",interp,"delta")
    !-- density / stiffness / gravity
    call param_dbl(PARAMS,i,"SolidModel","mue",mue,1.0d0)
    call param_dbl(PARAMS,i,"SolidModel","eta",eta,1.0d0)
    call param_dbl(PARAMS,i,"SolidModel","f",frequ,1.0d0)
    call param_dbl(PARAMS,i,"SolidModel","angle",AngleBeam,1.0d0)
    call param_dbl(PARAMS,i,"SolidModel","gravity",grav,0.0d0)
    !-- beam thickness
    call param_dbl(PARAMS,i,"SolidModel","t_beam",t_beam,0.05d0)
    !-- damping
    call param_dbl(PARAMS,i,"SolidModel","sigma",sigma,0.0d0)
    !-- timing
    call param_dbl(PARAMS,i,"SolidModel","T_release",T_release,0.0d0)
    call param_dbl(PARAMS,i,"SolidModel","tau",tau,0.0d0)
    call param_dbl(PARAMS,i,"SolidModel","N_smooth",N_smooth,3.0d0)
    call param_dbl(PARAMS,i,"SolidModel","L_span",L_span,1.0d0)
    call param_str(PARAMS,i,"SolidModel","has_cylinder",has_cylinder,"no")
    call param_dbl(PARAMS,i,"SolidModel","R_cylinder",R_cylinder,0.0d0)
    !-- time marching method for the solid
    call param_str(PARAMS,i,"SolidModel","TimeMethodSolid",TimeMethodSolid,"BDF2")
    
    select case (TimeMethodSolid)
      case ("RK4","CN2","BDF2","EI1","EE1")
        if (mpirank==0) write(*,*) "Solid solver is ", TimeMethodSolid
      case default
        if (mpirank==0) write(*,*) "Solid solver is UNDEFINED, using BDF2"
        TimeMethodSolid="BDF2"
    end select
    
    call param_str(PARAMS,i,"SolidModel","imposed_motion_leadingedge",&
         imposed_motion_leadingedge,"fixed_middle")
    call param_str(PARAMS,i,"SolidModel","infinite",infinite,"no")
    call param_str(PARAMS,i,"SolidModel","plate_shape",plate_shape,"rectangular")
    !-- grid spacing
    ds = 1.d0/dble(ns-1)
  endif

end subroutine get_params_solid




!-------------------------------------------------------------------------------
! Fetches a REAL VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_real: this is the parameter you were looking for
!-------------------------------------------------------------------------------
subroutine param_dbl (PARAMS, actual_lines, section, keyword, params_real, &
     defaultvalue)
  use vars
  use mpi
  implicit none
  character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
  character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
  character(len=strlen) :: value    ! returns the value
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  real (kind=pr) :: params_real, defaultvalue 
  integer actual_lines
  integer mpicode

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        if ( index(value,'*dx') == 0 ) then
          !-- value to be read is in absolute form (i.e. we just read the value)
          read (value, *) params_real    
          write (value,'(g10.3)') params_real
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        else
          !-- the value is given in gridpoints (e.g. thickness=5*dx)
          read (value(1:index(value,'*dx')-1),*) params_real
          params_real = params_real*max(dx,dy,dz)
          write (value,'(g10.3,"(=",g10.3,"*dx)")') params_real, params_real/max(dx,dy,dz)
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif 
     else
        write (value,'(g10.3)') defaultvalue
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_real = defaultvalue
     endif
  endif
  
  ! And then broadcast
  call MPI_BCAST( params_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode ) 
end subroutine param_dbl




!-------------------------------------------------------------------------------
! Fetches a STRING VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_string: this is the parameter you were looking for
!-------------------------------------------------------------------------------
subroutine param_str (PARAMS, actual_lines, section, keyword, &
     params_string, defaultvalue)
  use vars
  use mpi
  implicit none
  
  character(len=*), intent(in) :: section ! what section do you look for? for example [Resolution]
  character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
  character(len=strlen)  value    ! returns the value
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  character(len=strlen), intent (inout) :: params_string
  character(len=*), intent (in) :: defaultvalue 
  integer actual_lines
  integer mpicode
  
  !------------------
  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  !------------------  
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        params_string = value
        ! its a bit dirty but it avoids filling the screen with
        ! "nothing" anytime we check the runtime control file
        if (keyword.ne."runtime_control") then 
           write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
     else
        value = defaultvalue
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_string = defaultvalue
     endif
  endif

  ! And then broadcast
  call MPI_BCAST( params_string, strlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpicode) 
end subroutine param_str



!-------------------------------------------------------------------------------
! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find a vector, we return this and warn
! Output:
!       params_vector: this is the parameter you were looking for
!-------------------------------------------------------------------------------
subroutine param_vct (PARAMS, actual_lines, section, keyword, params_vector, &
     defaultvalue)
  use vars
  use mpi
  implicit none
  character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
  character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
  character(len=strlen)  value    ! returns the value
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  real (kind=pr) :: params_vector(1:3), defaultvalue(1:3)
  integer, intent(in) :: actual_lines
  integer :: mpicode

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        ! read the three values from the vector string
        read (value, *) params_vector
        write (value,'(3(g10.3,1x))') params_vector
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
     else
        write (value,'(3(g10.3,1x))') defaultvalue
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_vector = defaultvalue
     endif
  endif
  
  ! And then broadcast
  call MPI_BCAST( params_vector, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode ) 
end subroutine param_vct



!-------------------------------------------------------------------------------
! Fetches a INTEGER VALUED parameter from the PARAMS.ini file.
! Displays what it does on stdout (so you can see whats going on)
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
!       defaultvalue: if the we can't find the parameter, we return this and warn
! Output:
!       params_int: this is the parameter you were looking for
!-------------------------------------------------------------------------------
subroutine param_int(PARAMS, actual_lines, section, keyword, params_int,&
     defaultvalue)
  use mpi
  use vars
  implicit none

  character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
  character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
  character(len=strlen) ::  value    ! returns the value
  ! Contains the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  integer :: params_int, actual_lines, defaultvalue
  integer :: mpicode

  !------------------
  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  !------------------
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        read (value, *) params_int
        write (value,'(i7)') params_int
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
     else
        write (value,'(i7)') defaultvalue
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_int = defaultvalue
     endif
  endif

  ! And then broadcast
  call MPI_BCAST( params_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )  
end subroutine param_int





!-------------------------------------------------------------------------------
! Extracts a value from the PARAMS.ini file, which is in "section" and
! which is named "keyword"
! Input:
!       PARAMS: the complete *.ini file as array of characters*256
!               with maximum nlines=2048 lines (this value is set in share_vars)
!       actual_lines: number of lines in the params file
!                     because we don't need to loop over empty lines
!       section: the section we're looking for
!       keyword: the keyword we're looking for
! Output:
!       value: is a string contaiing everything between '=' and ';'
!              to be processed further, depending on the expected type
!              of variable (e.g. you read an integer from this string)
!-------------------------------------------------------------------------------
subroutine GetValue (PARAMS, actual_lines, section, keyword, value)
  use vars
  use mpi
  implicit none

  character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
  character(len=*), intent(in) :: keyword   ! what keyword do you look for? for example nx=128
  character(len=*), intent(inout) :: value   ! returns the value
  character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  integer, intent(in) ::  actual_lines   ! how many lines did you actually read?  
  integer :: i
  integer :: index1,index2
  logical :: foundsection

  foundsection = .false.
  value = ''

  !-- loop over the lines of PARAMS.ini file
  do i=1, actual_lines
     !-- ignore commented lines completely
     if ((PARAMS(i)(1:1).ne.'#').and.(PARAMS(i)(1:1).ne.';').and.&
         (PARAMS(i)(1:1).ne.'!')) then

        !-- does this line contain the "[section]" statement?
        if (index(PARAMS(i),'['//section//']')==1) then
          ! yes, it does
          foundsection = .true.
        elseif (PARAMS(i)(1:1) == '[') then
          ! we're already at the next section mark, so we leave the section we
          ! were looking for again
          foundsection = .false.          
        endif
                
        !-- we're inside the section we want
        if (foundsection) then
          if (index(PARAMS(i),keyword//'=')==1) then
            if (index(PARAMS(i),';')/=0) then
              ! found "keyword=" in this line, as well as the delimiter ";" 
              index1 = index(PARAMS(i),'=')+1
              index2 = index(PARAMS(i),';')-1
              value = PARAMS(i)(index1:index2)
              exit ! leave do loop
            else              
              ! found "keyword=" in this line, but delimiter ";" is missing.
              ! proceed, but this may fail (e.g. value=999commentary or 
              ! value=999\newline fails)
              write (*,'(A)') "Though found the keyword, I'm unable to find &
                   &value for variable --> "//trim(keyword)//" <-- &
                   &missing delimiter (;)"              
              index1 = index(PARAMS(i),'=')+1
              index2 = index(PARAMS(i),' ')-1
              value = PARAMS(i)(index1:index2)
              write (*,'(A)') "As you forgot to set ; at the end of your value,&
                    & I try to extract a value nonetheless. Proceed, with&
                    & fingers crossed (there is no guarantee this will work)"
              write (*,'(A)') "Extracted --->"//trim(value)//"<-----"       
              exit
            endif
          endif
        endif
        
     endif
  enddo
end subroutine getvalue

! Read the file paramsfile, count the lines (output in i) and put the
! text in PARAMS.
subroutine read_params_file(PARAMS,i,paramsfile, verbose)
  use mpi
  use vars
  implicit none
  
  integer,intent(out) :: i
  integer :: io_error
  ! This array will contain the ascii-params file
  character(len=strlen), dimension(1:nlines), intent(inout) :: PARAMS
  character(len=strlen) :: paramsfile ! this is the file we read the PARAMS from
  logical, intent(in) :: verbose

  ! Read in the params file (root only)
  io_error=0
  if (mpirank==0) then
     if (verbose) then
     write (*,*) "*************************************************"
     write (*,'(A,i3)') " *** info: reading params from "//trim(paramsfile)//" rank=", mpirank
     write (*,*) "*************************************************"
     endif
     i = 1
     open(unit=14,file=trim(adjustl(paramsfile)),action='read',status='old')    
     do while ((io_error==0).and.(i<=nlines))
        read (14,'(A)',iostat=io_error) PARAMS(i)  
        i = i+1
     enddo
     close (14)
     i = i-1 ! counted one too far
  endif

end subroutine read_params_file