!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! can be done in parallel. the flag --second order can be used for filtering
subroutine convert_pressure(help)
  use vars
  use p3dfft_wrapper
  use helpers
  use basic_operators
  use penalization
  use insect_module
  use solid_model
  use turbulent_inlet_module
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, fname_p, fname_ini
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk,nlk
  complex(kind=pr),dimension(:,:,:,:),allocatable :: workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, workr, vort
  real(kind=pr),dimension(:,:,:),allocatable :: press
  real(kind=pr) :: time, divu_max
  type(diptera) :: Insect
  type(solid), dimension(1) :: beams

  method="fsi"
  allocate(lin(1)) ! Set up the linear term

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pressure ux uy uz p"
    write(*,*) "./flusi -p --pressure PARAMS.ini ux uy uz p"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Compute pressure from velocity."
    write(*,*) " "
    write(*,*) " If the first argument is a parameter file, the code will read it and use it to"
    write(*,*) " construct the mask function, if possible. In some cases, where FSI eqns are solved"
    write(*,*) " it might not be possible to do so currently. Note if you plan on supplying the mask"
    write(*,*) " you need the solid body velocity field as well."
    write(*,*) " "
    write(*,*) " NOTE: the field has to be periodic, as we solve a poisson eqn here. It cannot work"
    write(*,*) "       on cropped fields."
    write(*,*) " "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3, fname_ini)
  if (index(fname_ini,".ini") /= 0) then
    ! first argument is an ini file.
    call get_command_argument(4, fname_ux)
    call get_command_argument(5, fname_uy)
    call get_command_argument(6, fname_uz)
    call get_command_argument(7, fname_p)
  else
    ! first argument is NOT an ini file
    call get_command_argument(3, fname_ux)
    call get_command_argument(4, fname_uy)
    call get_command_argument(5, fname_uz)
    call get_command_argument(6, fname_p)
    fname_ini = ""
  endif

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "Computing pressure from"
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(fname_p))
    write(*,*) "Parameter file (may be empty):"
    write(*,*) trim(adjustl(fname_ini))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if (fname_ini /= "") then
    call get_params(fname_ini, Insect, .true.)
  endif

  ! if the simulation was done using free-flight, we'd need to read rigidsolidsolver.t
  ! and interpolate the state vector etc (as done in --dry-run). TODO.
  if (Insect%BodyMotion=="free_flight") then
    call abort(1239, "This routine needs the mask/us array but cannot construct it in free-flight.")
  endif

  ! note this is actually redundant except for the TIME attribute
  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  neq = 3
  nd = 3
  nrw = 0
  ncw = 1 ! used in "pressure_from_uk_use_existing_mask"
  ga = ra ! since we do not use ghosts nodes
  gb = rb

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw))
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw))
  allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  allocate(press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
  allocate(mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  mask_color=0

  call read_single_file( fname_ux, u(:,:,:,1) )
  call read_single_file( fname_uy, u(:,:,:,2) )
  call read_single_file( fname_uz, u(:,:,:,3) )
  ! to Fourier space
  call fft3 (inx=u, outk=uk)

  ! read in turbulent inlet fields
  if (use_turbulent_inlet=="yes") then
    call init_turbulent_inlet ( )
  endif

  ! create mask
  call create_mask( time, Insect, beams )

  ! compute pressure
  call pressure_from_uk_use_existing_mask(time,u,uk,nlk,vort,workr,workc,press,Insect)

  ! now press contains the pressure (and we do not use ghost nodes, so ga=ra and gb=rb)
  call save_field_hdf5 ( time, fname_p, press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )

  ! clean up, bye bye captain.
  deallocate(u,uk,mask,us,workr,workc, vort,nlk,press)
  call fft_free()

end subroutine convert_pressure
