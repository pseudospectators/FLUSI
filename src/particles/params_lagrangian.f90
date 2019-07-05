subroutine get_params_lagrangian(PARAMS)
! Read individual parameter values from the PARAMS string for particles
  use mpi
  use vars
  use module_ini_files_parser_mpi
  use interpolation
  !use lagrangian_vars
  implicit none

  integer:: n_part_loc
  ! Contains the ascii-params file
  !character(len=strlen), dimension(1:nlines), intent(in) :: PARAMS
  type(inifile) :: PARAMS

  ! Lagrangian section
  call read_param_mpi(PARAMS,"Lagrangian","use_lagrangian",use_lagrangian,0)

   
  ! reclc is the number of bytes for I/O access direct
  call read_param_mpi(PARAMS,"Lagrangian","reclc",reclc,16)

  ! Total number of particles
  call read_param_mpi(PARAMS,"Lagrangian","npart",npart,1)
   
  ! Number of particles allowed by proc
  call read_param_mpi(PARAMS,"Lagrangian","nper_part",nper_part,100)
  n_part_loc=int(npart/mpisize)+(npart-int(npart/mpisize)*mpisize)
  npart_proc=int(dble(nper_part)/100.0*dble(n_part_loc))+n_part_loc

  if (mpisize==1) npart_proc=npart
  if (mpirank==0) print*,'read Lagrangian::npart_proc=',npart_proc

  ! interpolation method
  call read_param_mpi(PARAMS,"Lagrangian","lagr_interp",lagr_interp,"linear")

  ! Choice of Lagrangian quatities to compute
  call read_param_mpi(PARAMS,"Lagrangian","ilagr_acc",ilagr_acc,0)
  call read_param_mpi(PARAMS,"Lagrangian","ilagr_vort",ilagr_vort,0)

  ! itraj_test=1 -> save the first particle in ascci file 'traj_test'
  call read_param_mpi(PARAMS,"Lagrangian","itraj_test",itraj_test,0)

  ! Number of data to interpolate (interpolation module)
  if (nx==1) then
     n_interp=2*(1+ilagr_acc) + ilagr_vort
  else
     n_interp=3*(1+ilagr_acc+ilagr_vort)
  endif
  if (mpirank==0) print*,'read Lagrangian:: Number of data to interpolate n_interp=',n_interp

  ! Choose initial positions of the particles
  call read_param_mpi(PARAMS,"Lagrangian","inicond_lagrangian",inicond_lagrangian,"random")

  ! Saving 
  call read_param_mpi(PARAMS,"Lagrangian","tsave_part",tsave_part,9.d9)

end subroutine get_params_lagrangian
