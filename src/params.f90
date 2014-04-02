! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile) 
  use vars

  character (len=80) :: paramsfile  ! The file we read the PARAMS from
  integer :: i  
  character PARAMS(nlines)*256 ! this array will contain the ascii-params file
  logical :: exist1
  
  ! check if the specified file exists
  inquire ( file=paramsfile, exist=exist1 )
  
  if (((paramsfile=="").or.(exist1.eqv..false.)).and.(mpirank==0)) then
    write(*,*) "Please specify the params file!"
    write(*,*) "eg: mhd PARAMS or flusi PARAMS"
    stop
  endif
  
  ! Read the paramsfile and put the length i and the text in PARAMS
  call read_params_file(PARAMS,i,paramsfile,.true.)

  ! Get parameter values from PARAMS
  call get_params_common(PARAMS,i)

  select case(method(1:3))
     case("fsi") 
        ! Get fsi-specific parameter values from PARAMS
        call get_params_fsi(PARAMS,i)
     case("mhd") 
        ! Get mhd-specific parameter values from PARAMS
        call get_params_mhd(PARAMS,i)
  case default
     if(mpirank == 0) then
        write(*,*) "Error! Unkonwn method in get_params; stopping."
        stop
     end if
  end select
end subroutine get_params


! Read the file paramsfile, count the lines (output in i) and put the
! text in PARAMS.
subroutine read_params_file(PARAMS,i,paramsfile, verbose)
  use mpi
  use vars
  implicit none
  
  integer,intent(out) :: i
  integer :: io_error
  ! This array will contain the ascii-params file
  character,intent(inout) :: PARAMS(nlines)*256
  character(len=80) :: paramsfile ! this is the file we read the PARAMS from
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


! Read individual parameter values from the PARAMS string for the vars
! module.
subroutine get_params_common(PARAMS,i)
  use mpi
  use vars
  implicit none

  integer,intent(in) :: i
  character,intent(in) :: PARAMS(nlines)*256 ! Contains the ascii-params file

  ! Resolution section
  call GetValue_Int(PARAMS,i,"Resolution","nx",nx, 4)
  call GetValue_Int(PARAMS,i,"Resolution","ny",ny, 4)
  call GetValue_Int(PARAMS,i,"Resolution","nz",nz, 4)

  ! Time section
  call GetValue_Int(PARAMS,i,"Time","nt",nt, 9999999)
  ! This is the default value. ifort complains when just putting it in the call
  iTimeMethodFluid = "AB2" 
  call GetValue_String(PARAMS,i,"Time","iTimeMethodFluid", iTimeMethodFluid,&
       iTimeMethodFluid)  
  call GetValue_Real(PARAMS,i,"Time","Tmax",Tmax,1.d9)
  call GetValue_Real(PARAMS,i,"Time","CFL",cfl,0.1d0)
  call GetValue_Real(PARAMS,i,"Time","dt_fixed",dt_fixed,0.d0)

  ! Reynolds number section:
  call GetValue_Real(PARAMS,i,"ReynoldsNumber","nu",nu,1.d-2)  

  ! Initial conditions section
  inicond = "none"
  call GetValue_String(PARAMS,i,"InitialCondition","inicond",inicond, inicond)
  call GetValue_Real(PARAMS,i,"InitialCondition","omega1",omega1,0.d0) 
  ! if reading from file, which files?
  if (inicond=="infile") then 
    file_ux="none"
    call GetValue_String(PARAMS,i,"InitialCondition","file_ux",file_ux, file_ux)
    file_uy="none"
    call GetValue_String(PARAMS,i,"InitialCondition","file_uy",file_uy, file_uy)
    file_uz="none"
    call GetValue_String(PARAMS,i,"InitialCondition","file_uz",file_uz, file_uz)
    ! if running in MHD mode, we also need the B-field initialized
    if (method=="mhd") then
      file_bx="none"
      call GetValue_String(PARAMS,i,"InitialCondition","file_bx",file_bx, file_bx)
      file_by="none"
      call GetValue_String(PARAMS,i,"InitialCondition","file_by",file_by, file_by)
      file_bz="none"
      call GetValue_String(PARAMS,i,"InitialCondition","file_bz",file_bz, file_bz)
    endif
  endif

  ! Dealasing section
  call GetValue_Int(PARAMS,i,"Dealiasing","iDealias",iDealias, 1)

  ! Penalization section
  call GetValue_Int(PARAMS,i,"Penalization","iPenalization",iPenalization, 0)
  call GetValue_Int(PARAMS,i,"Penalization","iMoving",iMoving, 0)
  imask="none"
  call GetValue_String(PARAMS,i,"Penalization","iMask",iMask, iMask)
  iSmoothing="erf" ! std choice
  call GetValue_String(PARAMS,i,"Penalization","iSmoothing",iSmoothing,iSmoothing)
  call GetValue_Real(PARAMS,i,"Penalization","eps",eps, 1.d-2)
  call GetValue_Real(PARAMS,i,"Penalization","pseudoeps",pseudoeps, 1.d-2)
  call GetValue_Real(PARAMS,i,"Penalization","pseudodt",pseudodt, 1.d-2)
  call GetValue_Real(PARAMS,i,"Penalization","pseuderrmin",pseudoerrmin,3d-4)
  call GetValue_Real(PARAMS,i,"Penalization","pseuderrmax",pseudoerrmax,5d-4)

  ! Geometry section
  call GetValue_Real(PARAMS,i,"Geometry","xl",xl, 1.d0)
  call GetValue_Real(PARAMS,i,"Geometry","yl",yl, 1.d0)
  call GetValue_Real(PARAMS,i,"Geometry","zl",zl, 1.d0)
  call GetValue_Real(PARAMS,i,"Geometry","Size",length, 0.d0)
  call GetValue_Real(PARAMS,i,"Geometry","r1",r1,1.d0)
  call GetValue_Real(PARAMS,i,"Geometry","r2",r2,1.0681415d0)
  call GetValue_Real(PARAMS,i,"Geometry","r3",r3,1.206371d0)


  ! Saving section
  call GetValue_Int(PARAMS,i,"Saving","iDoBackup",iDoBackup, 1)
  call GetValue_Int(PARAMS,i,"Saving","iSaveVelocity",iSaveVelocity, 0) 
  call GetValue_Int(PARAMS,i,"Saving","iSavePress",iSavePress, 0)
  call GetValue_Int(PARAMS,i,"Saving","iSaveVorticity",iSaveVorticity, 1)  
  call GetValue_Int(PARAMS,i,"Saving","iSaveMask",iSaveMask, 0)
  call GetValue_Int(PARAMS,i,"Saving","iSaveXMF",iSaveXMF, 1) ! default is yes
  call GetValue_Real(PARAMS,i,"Saving","tsave",tsave, 9.d9)
  call GetValue_Real(PARAMS,i,"Saving","tintegral",tintegral,0.01d0)
  call GetValue_Int(PARAMS,i,"Saving","itdrag",itdrag,99999)
  
  !-- dry run, just the mask function
  dry_run_without_fluid="no" ! std choice
  call GetValue_String(PARAMS,i,"DryRun","dry_run_without_fluid",&
       dry_run_without_fluid,dry_run_without_fluid)
  if (dry_run_without_fluid=="yes") then
    write(*,*) "Attention! This is a dry run without fluid"
    write(*,*) "Deactivating all useless save-switches..."
    iDoBackup=0
    iSavePress=0
    iSaveVelocity=0
    iSaveVorticity=0
  endif
  
  ! Set other parameters (all procs)
  pi=4.d0 *datan(1.d0)
  ! scaling for FFTs
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl  
  ! lattice spacing is global
  dx=xl/dble(nx)
  dy=yl/dble(ny)
  dz=zl/dble(nz) 
end subroutine get_params_common


! Read individual parameter values from the PARAMS string for fsi.
subroutine get_params_fsi(PARAMS,i)
  use mpi
  use fsi_vars
  implicit none

  integer,intent(in) :: i
  character,intent(in) :: PARAMS(nlines)*256 ! Contains the ascii-params file
  
  ! ---------------------------------------------------
  ! penalization / cavity
  ! ---------------------------------------------------
  iCavity = "no"
  call GetValue_String(PARAMS,i,"Penalization","iCavity",iCavity,iCavity) 
  call GetValue_Int(PARAMS,i,"Penalization","cavity_size",cavity_size,0)
  call GetValue_Int(PARAMS,i,"Penalization","compute_forces",compute_forces,1)   
  call GetValue_Int(PARAMS,i,"Penalization","unst_corrections",unst_corrections,0)        
  iChannel = "no"
  call GetValue_String(PARAMS,i,"Penalization","iChannel",iChannel,iChannel) 
  if (iChannel=="0") iChannel="no" ! for downward compatibility with older ini files
  if (iChannel=="1") iChannel="xy" ! for downward compatibility with older ini files
  call GetValue_Real(PARAMS,i,"Penalization","thick_wall",thick_wall,0.2d0)
  call GetValue_Real(PARAMS,i,"Penalization","pos_wall",pos_wall,0.3d0)
  
  ! ---------------------------------------------------
  ! sponge
  ! ---------------------------------------------------
  iVorticitySponge = "no"
  call GetValue_String(PARAMS,i,"Sponge","iVorticitySponge",&
       iVorticitySponge,iVorticitySponge)  
  iSpongeType = "top_cover"
  call GetValue_String(PARAMS,i,"Sponge","iSpongeType",&
       iSpongeType,iSpongeType)         
  call GetValue_Real(PARAMS,i,"Sponge","eps_sponge",eps_sponge, 0.d0)   
  call GetValue_Int(PARAMS,i,"Sponge","sponge_thickness",sponge_thickness,0)
  
  
  ! ---------------------------------------------------
  ! Geometry section
  ! ---------------------------------------------------
  call GetValue_Real(PARAMS,i,"Geometry","x0",x0, 0.d0)
  call GetValue_Real(PARAMS,i,"Geometry","y0",y0, 0.d0)
  call GetValue_Real(PARAMS,i,"Geometry","z0",z0, 0.d0)
  
  ! ---------------------------------------------------
  ! Saving section
  ! ---------------------------------------------------
  call GetValue_Int(PARAMS,i,"Saving","iSaveSolidVelocity",iSaveSolidVelocity,0)
  save_only_one_period = "no"
  call GetValue_String(PARAMS,i,"Saving","save_only_one_period",&
       save_only_one_period,save_only_one_period)  

  ! ---------------------------------------------------
  ! MeanFlow section
  ! ---------------------------------------------------
  call GetValue_Int(PARAMS,i,"MeanFlow","iMeanFlow",iMeanFlow, 3)  
  call GetValue_Real(PARAMS,i,"MeanFlow","ux",uxmean, 1.d0) 
  call GetValue_Real(PARAMS,i,"MeanFlow","uy",uymean, 1.d0) 
  call GetValue_Real(PARAMS,i,"MeanFlow","uz",uzmean, 1.d0) 

  ! ---------------------------------------------------
  ! Insects section
  ! ---------------------------------------------------
  Insect%WingShape="none"
  call GetValue_String(PARAMS,i,"Insects","WingShape",&
       Insect%WingShape,Insect%WingShape)
  call GetValue_Real(PARAMS,i,"Insects","b_top",Insect%b_top, 0.d0) 
  call GetValue_Real(PARAMS,i,"Insects","b_bot",Insect%b_bot, 0.d0) 
  call GetValue_Real(PARAMS,i,"Insects","L_chord",Insect%L_chord, 0.d0) 
  call GetValue_Real(PARAMS,i,"Insects","L_span",Insect%L_span, 0.d0) 
  
  ! for string parameters, set the default in the first place
  Insect%FlappingMotion_right="none"
  call GetValue_String(PARAMS,i,"Insects","FlappingMotion_right",&
       Insect%FlappingMotion_right,Insect%FlappingMotion_right)
       
  Insect%FlappingMotion_left="none"
  call GetValue_String(PARAMS,i,"Insects","FlappingMotion_left",&
       Insect%FlappingMotion_left,Insect%FlappingMotion_left)
       
  Insect%BodyType="ellipsoid"
  call GetValue_String(PARAMS,i,"Insects","BodyType",&
       Insect%BodyType,Insect%BodyType)  
       
  Insect%HasEye="yes"
  call GetValue_String(PARAMS,i,"Insects","HasEye",&
       Insect%HasEye,Insect%HasEye)
       
  Insect%HasHead="yes"
  call GetValue_String(PARAMS,i,"Insects","HasHead",&
       Insect%HasHead,Insect%HasHead)    
       
  Insect%BodyMotion="yes"
  call GetValue_String(PARAMS,i,"Insects","BodyMotion",&
       Insect%BodyMotion,Insect%BodyMotion)        
       
  call GetValue_Real(PARAMS,i,"Insects","b_body",Insect%b_body, 0.1d0) 
  call GetValue_Real(PARAMS,i,"Insects","L_body",Insect%L_body, 1.d0)
  call GetValue_Real(PARAMS,i,"Insects","R_head",Insect%R_head, 0.1d0) 
  call GetValue_Real(PARAMS,i,"Insects","R_eye",Insect%R_eye, 0.d1) 
  call GetValue_Real(PARAMS,i,"Insects","distance_from_sponge",Insect%distance_from_sponge, 1.d0) 
  call GetValue_Real(PARAMS,i,"Insects","WingThickness",Insect%WingThickness, 4.0*dx) 
  
  ! position vector of the head
  Insect%x_head=(/0.5*Insect%L_body,0.d0,0.d0 /)
  call GetValue_Vector(PARAMS,i,"Insects","x_head",&
       Insect%x_head, Insect%x_head) 
  
  ! eyes
  Insect%x_eye_r=Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8*(/1,+1,1/)
  call GetValue_Vector(PARAMS,i,"Insects","x_eye_r",&
       Insect%x_eye_r, Insect%x_eye_r) 
       
  Insect%x_eye_l=Insect%x_head+sin(45.d0*pi/180.d0)*Insect%R_head*0.8*(/1,-1,1/)
  call GetValue_Vector(PARAMS,i,"Insects","x_eye_l",&
       Insect%x_eye_l, Insect%x_eye_l) 
       
  ! wing hinges (root points)    
  Insect%x_pivot_l=(/0.d0, +Insect%b_body, 0.d0 /)    
  call GetValue_Vector(PARAMS,i,"Insects","x_pivot_l",&
       Insect%x_pivot_l, Insect%x_pivot_l)  
       
  Insect%x_pivot_r=(/0.d0, -Insect%b_body, 0.d0 /)
  call GetValue_Vector(PARAMS,i,"Insects","x_pivot_r",&
       Insect%x_pivot_r, Insect%x_pivot_r)        
              
     
  Insect%smooth = 2.0*dz

  ! flag: read kinematics from file
  Insect%KineFromFile="no"
  call GetValue_String(PARAMS,i,"Insects","KineFromFile",&
       Insect%KineFromFile,Insect%KineFromFile)    
 
  ! Takeoff 
  call GetValue_Real(PARAMS,i,"Insects","x_takeoff",Insect%x_takeoff, 2.0d0)
  call GetValue_Real(PARAMS,i,"Insects","z_takeoff",Insect%z_takeoff, 0.86d0)
  call GetValue_Real(PARAMS,i,"Insects","mass_solid",&
       Insect%mass_solid, 54.414118839786745d0)
  call GetValue_Real(PARAMS,i,"Insects","gravity",&
       Insect%gravity, -0.055129281110537755d0)

  ! Legs model parameters
  call GetValue_Int(PARAMS,i,"Insects","ilegs",Insect%ilegs, 1)
  call GetValue_Real(PARAMS,i,"Insects","anglegsend",&
       Insect%anglegsend, 0.7853981633974483d0)
  call GetValue_Real(PARAMS,i,"Insects","kzlegsmax",&
       Insect%kzlegsmax,64.24974647375242d0)
  call GetValue_Real(PARAMS,i,"Insects","dzlegsmax",&
       Insect%dzlegsmax,0.2719665271966527d0)
  call GetValue_Real(PARAMS,i,"Insects","t0legs",&
       Insect%t0legs,0.13643141797265643d0)
  call GetValue_Real(PARAMS,i,"Insects","tlinlegs",&
       Insect%tlinlegs,0.3547216867289067d0)

  ! ---------------------------------------------------
  ! DONE..
  ! ---------------------------------------------------
       
  lin(1)=nu

  if (mpirank==0) then
     write (*,*) "*************************************************"
     write (*,'(A,i3)') " *** DONE READING PARAMETERS"
     write (*,*) "*************************************************"
  endif
end subroutine get_params_fsi


! Read individual parameter values from the PARAMS string for fmhd
subroutine get_params_mhd(PARAMS,i)
  use mpi
  use mhd_vars
  implicit none

  integer,intent(in) :: i
  character,intent(in) :: PARAMS(nlines)*256 ! Contains the ascii-params file

  ! MHD section
  call GetValue_Real(PARAMS,i,"MHD","eta",eta,4.5d-2)

  ! MHDGeometry section
  call GetValue_Real(PARAMS,i,"MHDGeometry","b0",b0,4.5d0)
  call GetValue_Real(PARAMS,i,"MHDGeometry","bc",bc,3.88888888888d0)
  call GetValue_Real(PARAMS,i,"MHDGeometry","ay",ay,1.0d0)

  ! Saving section
  call GetValue_Int(PARAMS,i,"Saving","iSaveMagneticField",&
       iSaveMagneticField, 0)
  call GetValue_Int(PARAMS,i,"Saving","iSaveCurrent",iSaveCurrent, 0)

  lin(1)=nu
  lin(2)=eta
end subroutine get_params_mhd




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
subroutine GetValue_real (PARAMS, actual_lines, section, keyword, params_real, &
     defaultvalue)
  use vars
  use mpi
  implicit none
  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  real (kind=pr) :: params_real, defaultvalue 
  integer actual_lines
  integer mpicode

  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        read (value, *) params_real    
        write (value,'(g10.3)') params_real
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
     else
        write (value,'(g10.3)') defaultvalue
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_real = defaultvalue
     endif
  endif
  
  ! And then broadcast
  call MPI_BCAST( params_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode ) 
end subroutine GetValue_real




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
subroutine GetValue_string (PARAMS, actual_lines, section, keyword, &
     params_string, defaultvalue)
  use vars
  use mpi
  implicit none
  
  character section*(*) ! what section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  character (len=80), intent (inout) :: params_string
  character (len=80), intent (inout) :: defaultvalue 
  integer actual_lines
  integer mpicode

  !------------------
  ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
  !------------------  
  if (mpirank==0) then
     call GetValue(PARAMS, actual_lines, section, keyword, value)
     if (value .ne. '') then
        params_string = value
        ! its a bit dirty but it avoids filling the screen with "nothing" anytime we check
        ! the runtime control file
        if (keyword.ne."runtime_control") then 
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
     else
        write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
             " (THIS IS THE DEFAULT VALUE!)"
        params_string = defaultvalue
     endif
  endif

  ! And then broadcast
  call MPI_BCAST( params_string, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpicode) 
end subroutine GetValue_string



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
subroutine GetValue_vector (PARAMS, actual_lines, section, keyword, params_vector, &
     defaultvalue)
  use vars
  use mpi
  implicit none
  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  real (kind=pr) :: params_vector(1:3), defaultvalue(1:3)
  integer actual_lines
  integer mpicode

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
end subroutine GetValue_vector



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
subroutine GetValue_Int(PARAMS, actual_lines, section, keyword, params_int,&
     defaultvalue)
  use mpi
  use vars
  implicit none

  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character (len=80)  value    ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  integer params_int, actual_lines, defaultvalue
  integer mpicode

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
end subroutine GetValue_Int





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
subroutine GetValue (PARAMS, actual_lines, section, keyword, value)
  use vars
  use mpi
  implicit none

  character section*(*) ! What section do you look for? for example [Resolution]
  character keyword*(*)   ! what keyword do you look for? for example nx=128
  character value*(*)   ! returns the value
  character PARAMS(nlines)*256  ! this is the complete PARAMS.ini file
  integer actual_lines   ! how many lines did you actually read?  
  integer :: maxline = 256  ! how many characters per line?
  integer i, j,k
  logical foundsection

  foundsection = .false.
  value = ''

  !------------------------------------------------------------------
  do i=1, actual_lines     ! loop over the lines of PARAMS.ini file
     if ((PARAMS(i)(1:1).ne.'#').and.&
          (PARAMS(i)(1:1).ne.';').and.&
          (PARAMS(i)(1:1).ne.'!')) then   ! ignore commented lines compleetly


        if (PARAMS(i)(1:1) == '[') then   ! the first char would have to be '['
           do j = 2, maxline    ! then we look fot the corrresponding ']'
              if (PARAMS(i)(j:j) == ']') then  ! we found it
                 if (section == PARAMS(i)(2:j-1)) then ! is this the section we"re looking for?
                    foundsection = .true.   ! yes, it is
                    exit
                 endif
              endif
           enddo
        else      
           if (foundsection .eqv. .true.) then  ! yes we found the section, now we're looking for the keyword
              do j=1, maxline    ! scan the line
                 if (PARAMS(i)(j:j) == '=') then  ! found the '='
                    if (keyword == PARAMS(i)(1:j-1)) then ! is this the keyword you're looking for?
                       do k = j+1, maxline   ! everything behind the '=' and before ';' is the value
                          if (PARAMS(i)(k:k) == ';') then  ! found the delimiter
                             value = PARAMS(i)(j+1:k-1)  ! value is between '=', and ';'   
                             exit
                          endif
                       enddo
                       if ((value == '').and.(mpirank==0)) then
                          write (*,'(A)') "??? Though found the keyword, I'm unable to find value for variable --> "& 
                               //trim(keyword)//" <-- maybe missing delimiter (;)?"    
                       endif
                    endif
                 endif
              enddo
           endif
        endif

     endif
  enddo
end subroutine getvalue
