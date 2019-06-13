!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! can be done in parallel. the flag --second order can be used for filtering
subroutine convert_pressure(help)
  use vars
  use p3dfft_wrapper
  use module_helpers
  use basic_operators
  use penalization
  use module_insects
  use solid_model
  use turbulent_inlet_module
  use hdf5_wrapper

  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, fname_p, fname_ini
  character(len=strlen) :: fname_mask, fname_usx, fname_usy, fname_usz
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk,nlk
  complex(kind=pr),dimension(:,:,:,:),allocatable :: workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, workr, vort
  real(kind=pr),dimension(:,:,:),allocatable :: press
  real(kind=pr) :: time, tmp(1)
  type(diptera) :: Insect
  type(solid), dimension(1) :: beams

  method="fsi"
  allocate(lin(1)) ! Set up the linear term

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --pressure ux uy uz p"
    write(*,*) "./flusi -p --pressure ux uy uz p mask usx usy usz"
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
    fname_mask = ""
    fname_usx = ""
    fname_usy = ""
    fname_usz = ""
  else
    ! first argument is NOT an ini file
    call get_command_argument(3, fname_ux)
    call get_command_argument(4, fname_uy)
    call get_command_argument(5, fname_uz)
    call get_command_argument(6, fname_p)
    call get_command_argument(7, fname_mask)
    call get_command_argument(8, fname_usx)
    call get_command_argument(9, fname_usy)
    call get_command_argument(10, fname_usz)
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
    write(*,*) trim(adjustl(fname_mask))
    write(*,*) trim(adjustl(fname_usx))
    write(*,*) trim(adjustl(fname_usy))
    write(*,*) trim(adjustl(fname_usz))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(fname_p))
    write(*,*) "Parameter file (may be empty):"
    write(*,*) trim(adjustl(fname_ini))
    write(*,'(80("-"))')
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if (fname_mask /= "") then
    call check_file_exists( fname_mask )
    call check_file_exists( fname_usx )
    call check_file_exists( fname_usy )
    call check_file_exists( fname_usz )
  endif

  if (fname_ini /= "") then
      call get_params(fname_ini, Insect, .true.)
      if (iMask == "Insect") then
          call insect_init(time, fname_ini, Insect, .false., "", (/xl,yl,zl/), nu, &
          dx, periodic=periodic)
          ! max color
          if (Insect%second_wing_pair) then
            endcolor = 5
          else
            endcolor = 3
          endif
      endif
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
  nrw = 1
  ncw = 1 ! used in "pressure_from_uk_use_existing_mask"
  if (iVorticitySponge=="yes")  ncw=3
  ga = ra ! since we do not use ghosts nodes
  gb = rb

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  u = 0.0d0
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  uk = 0.0d0
  allocate(mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  mask = 0.0d0
  allocate(us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  us = 0.0d0
  allocate(workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw))
  workr = 0.0d0
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw))
  workc = 0.0d0
  allocate(vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd))
  vort = 0.0d0
  allocate(nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  nlk  = 0.0d0
  allocate(press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
  press = 0.0d0
  allocate(mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  mask_color=0

  call read_single_file( fname_ux, u(:,:,:,1) )
  call read_single_file( fname_uy, u(:,:,:,2) )
  call read_single_file( fname_uz, u(:,:,:,3) )
  ! to Fourier space
  call fft3 (inx=u, outk=uk)

  ! create mask only if possible (i.e. if parameters are given.)
  if (fname_ini /= "") then
    ! read in turbulent inlet fields
    if (use_turbulent_inlet=="yes") then
      call init_turbulent_inlet ( )
    endif

    ! create mask
    call create_mask( time, Insect, beams )
  endif

  ! if files are specified, we read the mask function from file, as well as the solid
  ! velocity field.
  if (fname_mask /= "") then
    call read_single_file( fname_mask, mask )
    call read_single_file( fname_usx, us(:,:,:,1) )
    call read_single_file( fname_usy, us(:,:,:,2) )
    call read_single_file( fname_usz, us(:,:,:,3) )

    ! fetch value of penalization parameter from mask file. It should be stored
    ! there for every version of flusi. Note we do not need to divide by eps here, since
    ! this is done in RHS. It is the first time we need
    ! the epsi parameter since 2014 (now is 2018), so it was not included in
    ! fetch_attributes
    call read_attribute( fname_mask, "mask", "epsi", tmp)
    eps = tmp(1)
    if (root) write(*,'("Value of eps=C_eta=",g15.8)') eps
  endif

  ! compute pressure
  call pressure_from_uk_use_existing_mask(time,u,uk,nlk,vort,workr,workc,press,Insect)

  ! remove mean
  press = press - mpisum( sum(press) )/(dble(nx)*dble(ny)*dble(nz))

  ! now press contains the pressure (and we do not use ghost nodes, so ga=ra and gb=rb)
  call save_field_hdf5 ( time, fname_p, press(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )

  ! clean up, bye bye captain.
  deallocate(u,uk,mask,us,workr,workc, vort,nlk,press)
  call fft_free()

end subroutine convert_pressure
