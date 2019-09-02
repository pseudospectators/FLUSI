subroutine convert_to_wing_system(help)
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use module_insects
  use ghosts
  use interpolation

  implicit none
  logical, intent(in) :: help

  real(kind=pr)          :: t1,t2
  real(kind=pr)          :: time
  integer                :: ix, iy, iz, mpicode, i, nargs, nx1, nx2, ny1, ny2, nz1, nz2
  character (len=strlen) :: infile, outfile, whichwing
  ! arrays needed for interpolation
  real(kind=pr),dimension(:,:,:,:),allocatable :: u_org,u_interp
  real(kind=pr),dimension(:,:),allocatable :: u_buffer, u_buffer2
  ! this is the insect we're using (object oriented)
  type(diptera) :: Insect
  ! this is the solid model beams:
  type(solid), dimension(1:nBeams) :: beams
  real(kind=pr) :: x_wing(1:3),x_glob(1:3),M_wing_r(1:3,1:3),M_wing_l(1:3,1:3),M_body(1:3,1:3)
  real(kind=pr) :: M_wing_r2(1:3,1:3),M_wing_l2(1:3,1:3)
  real(kind=pr) :: u1,u2, x_body(1:3), u(1:3), rot_rel_wing_w(1:3)
  logical :: wing_system=.false., vector=.false., relative_velocity = .false.
  integer(kind=2) :: idw

  ! Set method information in vars module.
  method="fsi" ! We are doing fluid-structure interactions
  nf=1    ! We are evolving one field (that means 1 integrating factor)
  nd=3*nf ! The one field has three components.
  neq=nd  ! number of equations, can be higher than 3 if using passive scalar
  nrw=1   ! number of real valued work arrays
  ncw=1   ! number of complex values work arrays (decide that later)
  nrhs=2  ! number of right-hand side registers

  nargs = command_argument_count()

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p [--convert-to-wing-system|--convert-to-body-system] [--scalar|--vector] &
    &PARAMS.ini input_0000.h5 output_0000.h5 [--left|--right] [-x=100:200 -y=256:512 -z30:990] [--relative-velocity]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --convert-to-wing-system --scalar PARAMS.ini input_0000.h5 output_0000.h5 [--left|--right]"
    write(*,*) "./flusi -p --convert-to-body-system --scalar PARAMS.ini input_0000.h5 output_0000.h5"
    write(*,*) "./flusi -p --convert-to-wing-system --vector PARAMS.ini inx_0000.h5 iny_0000.h5 &
                &inz_0000.h5 outx_0000.h5 outy_0000.h5 outz_0000.h5  [--left|--right]"
    write(*,*) "./flusi -p --convert-to-body-system --vector PARAMS.ini inx_0000.h5 iny_0000.h5 &
                &inz_0000.h5 outx_0000.h5 outy_0000.h5 outz_0000.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! Convert a file to the (left) wings coordinate system"
    write(*,*) "! We read a field from file, which we assume to be a scalar field (or vector component)"
    write(*,*) "! given as q(x,y,z) where x,y,z are understood in the global coordinate system"
    write(*,*) "! "
    write(*,*) "! We now want q(xw,yw,zw) in the wing coordinate system. The output field has the same"
    write(*,*) "! spacing dx,dx,dz as the input field, and the coordinates are "
    write(*,*) "! -xl/2 <= xw <= xl/2"
    write(*,*) "! -yl/4 <= yw <= 3*yl/4"
    write(*,*) "! -zl/2 <= zw <= zl/2"
    write(*,*) "! "
    write(*,*) "! --right --left chooses left or right wing, if --convert-to-wing-system. default is left wing"
    write(*,*) "! "
    write(*,*) "! if you pass one of -x100:200, only a subset with the bounds you specify will be computed."
    write(*,*) "! This can speed up the computation significantly, provided you're interested only in narrow region"
    write(*,*) "! around the wing. "
    write(*,*) "! "
    write(*,*) "! For each of these points in the wing system, we compute the corresponding global coordinate"
    write(*,*) "! and interpolate the input field at this point. Note we require the PARAMS file to know what"
    write(*,*) "! motion protocoll to use for wings AND body. In the free-flight case, we require rigidsolidsolver.t to be present."
    write(*,*) "! The time is read from the input file."
    write(*,*) "! "
    write(*,*) "! Note the wing span axis is y, while the chord is x. z is wing normal."
    write(*,*) "! "
    write(*,*) "! Note: this routine can be used for vectors and scalars. The vectors are rotated to the new system as"
    write(*,*) "! well"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  !-----------------------------------------------------------------------------
  ! Read input parameters
  !-----------------------------------------------------------------------------
  allocate(lin(nf)) ! Set up the linear term
  if (root) write(*,'(A)') '*** info: Reading input data...'

  ! do we deal with a vector or a scalar?
  call get_command_argument(3,infile)

  if (infile == "--vector") then
      vector = .true.
      ! we allocate (and read) three components:
      nd = 3
      if (root) write(*,'(A)') 'we run in vector mode (3 components)'

  elseif (infile == "--scalar") then
      vector = .false.
      ! routine applied to scalar, one field only
      nd = 1
      if (root) write(*,'(A)') 'we run in scalar mode'

  else
      call abort(4556, "flag is unkown. --scalar or --vector")

  endif

  ! get filename of PARAMS file from command line
  call get_command_argument(4,infile)
  ! read all parameters from that file
  call get_params(infile,Insect,.true.)
  if (root) write(*,'("The resolution given by params file is: ",3(i4,1x))') nx,ny,nz

  !-----------------------------------------------------------------------------
  ! ghost points. we need that for interpolation
  !-----------------------------------------------------------------------------
  ng=1 ! one ghost points
  if (root) write(*,'("Set up ng=",i1," ghost points")') ng
  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()
  ! Setup communicators used for ghost point update
  call setup_cart_groups

  !-----------------------------------------------------------------------------
  ! Allocate memory:
  !-----------------------------------------------------------------------------
  allocate( u_org(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) )
  ! interpolated field has ghost nodes
  allocate( u_interp(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nd) )
  allocate( u_buffer(0:nx-1,1:nd),u_buffer2(0:nx-1,1:nd) )

  if (iMask/="Insect") then
    call abort(7774,"transforming to wing system makes no sense if not applied to an insect")
  endif

  ! initialize insect
  call insect_init(time, infile, Insect, .false., "", (/xl,yl,zl/), nu, dx, periodic=periodic)
  ! max color
  if (Insect%second_wing_pair) then
    endcolor = 5
  else
    endcolor = 3
  endif

  !*****************************************************************************
  ! main (active) part of this postprocessing tool
  !*****************************************************************************
  ! read in the input file(s) to be transformed
  do i = 1, nd
    call get_command_argument(i+4, infile)
    call read_single_file ( infile, u_org(:,:,:,i) )
  enddo
  ! call fetch attributes (even though we knew the resolution in advance) since we
  ! want to know what time the data is at
  call fetch_attributes(infile, nx, ny, nz, xl, yl, zl, time, nu, origin)

  ! synchronize ghosts (required for interpolation)
  u_interp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),:) = u_org
  call synchronize_ghosts( u_interp, nd )

  !-----------------------------------------------------------------------------
  ! fetch current motion state of the insect
  !-----------------------------------------------------------------------------
  ! fetch insect state
  call Update_Insect( time, Insect )

  if (Insect%BodyMotion == "free_flight") then
      if (root) write(*,*) "The insects body motion is free_flight, therefore we try to read the"
      if (root) write(*,*) "insect state vector from the file rigidsolidsolver.t"
      if (root) write(*,*) "from the state vector we can construct the body rotation matrix"
      call read_insect_STATE_from_file(time, Insect, 'rigidsolidsolver.t', .true.)
  endif

  call BodyMotion (time, Insect)
  idw = 1
  call FlappingMotionWrap (time, Insect, idw)
  idw = 2
  call FlappingMotionWrap (time, Insect, idw)
  if (Insect%second_wing_pair) then
    idw = 3
    call FlappingMotionWrap (time, Insect, idw)
    idw = 4
    call FlappingMotionWrap (time, Insect, idw)
  endif
  call StrokePlane (time, Insect)
  call body_rotation_matrix( Insect, M_body )
  call wing_right_rotation_matrix( Insect, M_wing_r )
  call wing_left_rotation_matrix( Insect, M_wing_l )
  if (Insect%second_wing_pair) then
    call wing_right2_rotation_matrix( Insect, M_wing_r2 )
    call wing_left2_rotation_matrix( Insect, M_wing_l2 )
  endif

  u_org = -9.9d10

  nz1 = 0
  nz2 = nz-1
  ny1 = 0
  ny2 = ny-1
  nx1 = 0
  nx2 = nx-1

  ! look for -x=A:B -y=A:B -z=A:B arguments (they can be all there or none)
  do i = 1, nargs
      call get_command_argument(i, infile)
      if (index(infile,'-x=') /= 0) then
          infile = infile( index(infile,'=')+1:len_trim(infile) )
          read (infile(1:index(infile,':')-1) ,*) nx1
          read (infile(index(infile,':',.true.)+1:len_trim(infile)),*) nx2
      endif

      if (index(infile,'-y=') /= 0) then
          infile = infile( index(infile,'=')+1:len_trim(infile) )
          read (infile(1:index(infile,':')-1) ,*) ny1
          read (infile(index(infile,':',.true.)+1:len_trim(infile)),*) ny2
      endif

      if (index(infile,'-z=') /= 0) then
          infile = infile( index(infile,'=')+1:len_trim(infile) )
          read (infile(1:index(infile,':')-1) ,*) nz1
          read (infile(index(infile,':',.true.)+1:len_trim(infile)),*) nz2
      endif

      if (index(infile,"--relative-velocity") /= 0) relative_velocity = .true.

      if (index(infile,"--convert-to-wing-system") /= 0) wing_system = .true.

      if (index(infile,"--convert-to-body-system") /= 0) wing_system = .false.

      if (index(infile,"--left") /= 0 .or. index(infile,"--right") /= 0 .or. &
          index(infile,"--left2") /= 0 .or. index(infile,"--right2") /= 0) then
          whichwing = infile
      endif
  enddo

  if (root) then
      write(*,'(40("~+"))')
      write(*,'("OPERATING AT TIME=",es15.8)') time
      write(*,*) "Computing on x=", nx1, nx2
      write(*,*) "Computing on y=", ny1, ny2
      write(*,*) "Computing on z=", nz1, nz2
      write(*,*) "Relative velocity: ", relative_velocity
      if (wing_system) then
          write(*,*) "we convert to wing system"
          write(*,*) "Chosen wing is: "//trim(adjustl(whichwing))
      else
          write(*,*) "we convert to body system"
      endif
      write(*,'(40("~+"))')
  endif

  do iz = nz1, nz2
      do iy = ny1, ny2
          do ix = nx1, nx2
              if (wing_system) then
                  ! define the position in the wing coordinate system (we look for u in this
                  ! coordinate system, so our output matrix is to be understood in this)
                  x_wing = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /) - (/xl/2.0d0,yl/4.0d0,zl/2.0d0/)
                  ! x_wing = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /) - (/xl/2.0d0,yl/2.0d0,zl/2.0d0/)

                  ! compute global coordinate
                  if (whichwing == "--right2") then
                      x_glob = Insect%xc_body_g+ matmul(transpose(M_body), &
                                (matmul(transpose(M_wing_r2),x_wing)+Insect%x_pivot_r2_b) )
                  elseif (whichwing == "--left2") then
                      x_glob = Insect%xc_body_g+ matmul(transpose(M_body), &
                                (matmul(transpose(M_wing_l2),x_wing)+Insect%x_pivot_l2_b) )
                  elseif (whichwing == "--right") then
                      x_glob = Insect%xc_body_g+ matmul(transpose(M_body), &
                                (matmul(transpose(M_wing_r),x_wing)+Insect%x_pivot_r_b) )
                  else
                      x_glob = Insect%xc_body_g+ matmul(transpose(M_body), &
                                (matmul(transpose(M_wing_l),x_wing)+Insect%x_pivot_l_b) )
                  endif
              else
                  ! transform to body system
                  x_body = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /) - (/xl,yl,zl/)/2.d0
                  x_glob = Insect%xc_body_g + matmul(transpose(M_body), x_body)
              endif

              ! periodization (note we cannot use periodize_coordinate as it is centeres around the midpoint)
              if (x_glob(1)<0.d0) x_glob(1)=x_glob(1)+xl
              if (x_glob(2)<0.d0) x_glob(2)=x_glob(2)+yl
              if (x_glob(3)<0.d0) x_glob(3)=x_glob(3)+zl

              if (x_glob(1)>=xl-dx) x_glob(1)=x_glob(1)-xl
              if (x_glob(2)>=yl-dy) x_glob(2)=x_glob(2)-yl
              if (x_glob(3)>=zl-dz) x_glob(3)=x_glob(3)-zl

              ! interpolate the value
              do i = 1,nd
                  u_buffer(ix,i) = trilinear_interp( (/dble(ga(1))*dx, dble(ga(2))*dy, dble(ga(3))*dz/), &
                  (/dx,dy,dz/),u_interp(:,:,:,i), x_glob, periodic=.false.  )
              enddo
          enddo

          ! at this point, all procs have interpolated a value in u_buffer. most of them
          ! are -9.9E10, which means the CPU does not have the data required for interpolation
          call MPI_ALLREDUCE(u_buffer,u_buffer2,nx*nd,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)

          ! now, all CPU have the line in x-direction with the properly interpolated values. however,
          ! only one of them actually holds this data. NOw we note that in FLUSI, we never split
          ! the x-coordinate, so this direction is ALWAYS contiguous
          if ( on_proc( (/0,iy,iz/) ) ) then
              ! only the responsible rank does this
              u_org(:,iy,iz,1:nd) = u_buffer2
          end if
      enddo
  enddo

  ! for some points, no value can be interpolated, since they do not exist ( happens for example at the corners)
  ! set 0 at these points:
  where (u_org < -9.0d10)
      u_org = 0.d0
  end where

  ! for now, interpolation has been done, i.e. we rotated the camera. now, if we are dealing with a
  ! vector, this is not enough, since we do not only shift, but also rotate the individual vectors.
  ! this is done now:

  if (vector) then
      do iz = ra(3), rb(3)
          do iy = ra(2), rb(2)
              do ix = ra(1), rb(1)
                  ! velocity in global coordinate system
                  u = (/u_org(ix,iy,iz,1), u_org(ix,iy,iz,2), u_org(ix,iy,iz,3)/)
                  x_wing = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /) - (/xl/2.0d0,yl/4.0d0,zl/2.0d0/)

                  if (wing_system) then
                      ! go to wing system
                      if (whichwing == "--right2") then
                          u = matmul(M_wing_r2, matmul(M_body,u))
                          rot_rel_wing_w = Insect%rot_rel_wing_r2_w
                      elseif (whichwing == "--left2") then
                          u = matmul(M_wing_l2, matmul(M_body,u))
                          rot_rel_wing_w = Insect%rot_rel_wing_l2_w 
                      elseif (whichwing == "--right") then
                          u = matmul(M_wing_r, matmul(M_body,u))
                          rot_rel_wing_w = Insect%rot_rel_wing_r_w
                      else
                          u = matmul(M_wing_l, matmul(M_body,u))
                          rot_rel_wing_w = Insect%rot_rel_wing_l_w
                      endif

                      if (relative_velocity) then
                          ! if applied to us, i.e. the input vector is the solid body velociy,
                          ! then this line kills all velocity (plausibility test)
                          u_org(ix,iy,iz,1:3) = u - cross(rot_rel_wing_w, x_wing)
                      else
                          u_org(ix,iy,iz,1:3) = u
                      endif
                  else
                      ! go to body system
                      u_org(ix,iy,iz,1:3) = matmul(M_body, u)
                  endif

              enddo
          enddo
      enddo
  endif


  !-----------------------------------------------------------------------------
  ! done, save output
  !-----------------------------------------------------------------------------
  if (vector .eqv. .true.) then
    ! vector case
    do i = 1, nd
      call get_command_argument(i+7,outfile)
      if (root) write(*,*) "vector output will be written to "//outfile
      call save_field_hdf5(time, outfile, u_org(:,:,:,i))
    enddo
  else
    ! scalar case
    call get_command_argument(6,outfile)
    if (root) write(*,*) "scalar output will be written to "//trim(adjustl(outfile))
    call save_field_hdf5(time, trim(adjustl(outfile)), u_org(:,:,:,1))
  endif

  !-----------------------------------------------------------------------------
  ! Deallocate memory
  !-----------------------------------------------------------------------------
  deallocate(u_org, u_interp, u_buffer, u_buffer2)
  ! Clean insect (the globally stored arrays for Fourier coeffs etc..)
  call insect_clean(Insect)
end subroutine convert_to_wing_system
