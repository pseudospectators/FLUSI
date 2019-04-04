!-------------------------------------------------------------------------------
! Wrapper for different postprocessing tools
!-------------------------------------------------------------------------------
subroutine postprocessing()
  use vars
  use mpi
  implicit none
  character(len=strlen) :: postprocessing_mode, filename, key1,key2, mode
  logical :: help
  real(kind=pr) :: t1
  help = .false.
  t1=MPI_wtime()

  ! the first argument is either "-p" or "-h"
  call get_command_argument(1,mode)

  ! the second argument tells us what to do
  call get_command_argument(2,postprocessing_mode)

  ! if mode is help, then we'll call the routines with help=.true. and skip
  ! all other output
  if ((mode=="-h").or.(mode=="--help")) then
    ! we'll just show the help for the routine (we even skip the header)
    help=.true.
  else
    ! we'll actually do something, postprocessing
    if (mpirank==0) then
      write(*,'(80("~"))')
      write(*,'(A)') "~~~ FLUSI is running in postprocessing mode"
      write(*,'(80("~"))')
    endif
    ! show the call from the command line in output
    call postprocessing_ascii_header( 6 )
  endif

  !-----------------
  ! check what to do
  !-----------------
  select case (postprocessing_mode)
  case ("--force")
      call post_force(help)
  case ("--POD")
      call POD(help)
  case ("--divergence")
    call post_div(help)
  case ("--flexible-wing-mask")
    call flexible_wing_mask(help)
  case ("--remove-mean")
    call remove_mean(help)
  case ("--coherent-vortex-extraction","--coherent-scalar-extraction","--CVE","--CSE")
    call post_CVE(help)
  case ("--readwrite")
    call readwrite(help)
  case ("--crop")
    call crop(help)
  case ("--uCT-assemble")
    call uCT_assemble_HDF5(help)
  case ("--transpose-test")
    call transpose_test(help)
  case ("--pressure-force")
    call pressure_force(help)
  case ("--convert-to-wing-system","--convert-to-body-system")
    call convert_to_wing_system(help)
  case ("--zero-padd")
    call post_zero_padd(help)
  case ("--force-decomp")
    call force_decomposition(help)
  case ("--stl2dist")
    call stl2dist(help)
  case ("--dist2mask","--dist2chi")
    call dist2chi(help)
  case ("--simple-field-operation")
    call simple_field_operation(help)
  case ("--cp")
    call copy_hdf_file(help)
  case ("--keyvalues")
    call get_command_argument(3,filename)
    call keyvalues (filename)
  case ("--compare-keys")
    call get_command_argument(3,key1)
    call get_command_argument(4,key2)
    call compare_key (key1,key2)
  case ("--compare-timeseries")
    call compare_timeseries(help)
  case ("--vorticity","--vor-abs","--vorticity-FD","--vor-abs-FD","--Q","--Q-FD")
    call convert_vorticity(help)
  case ("--vor2u")
    call convert_velocity(help)
  case ("--hdf2bin")
    call convert_hdf2bin(help)
  case ("--bin2hdf")
    call convert_bin2hdf(help)
  case ("--p2Q")
    call pressure_to_Qcriterion(help)
  case ("--extract-subset")
    call extract_subset(help)
  case ("--time-avg")
    call time_avg_HDF5(help)
  case ("--upsample")
    call upsample(help)
  case ("--spectrum")
    call post_spectrum(help)
  case ("--turbulence-analysis")
    call turbulence_analysis(help)
  case ("--field-analysis")
    call field_analysis(help)
  case ("--pressure")
    call convert_pressure(help)
  case ("--TKE-mean")
    call tke_mean(help)
  case ("--pointcloud2mask")
    call pointcloud2mask(help)
  case ("--max-over-x")
    call max_over_x(help)
  case ("--mean-over-x-subdomain")
    call mean_over_x_subdomain(help)
  case ("--mean-2D")
    call mean_2D(help)
  case ("--set-hdf5-attribute")
    call set_hdf5_attribute(help)
  case ("-ux-from-uyuz")
    call ux_from_uyuz(help)
  case ("--check-params-file")
    call check_params_file(help)
  case ("--magnitude")
    call magnitude_post(help)
  case ("--energy")
    call energy_post(help)
  case ("--helicity")
    call post_helicity(help)
  case ("--smooth-inverse-mask")
    call post_smooth_mask(help)
  case ("--gradient")
    call post_grad(help)
  case default
    if (root) then
      write(*,*) "Available Postprocessing tools are:"
      write(*,'(80("~"))')
      write(*,*) "--bin2hdf"
      write(*,*) "--check-params-file"
      write(*,*) "--coherent-vortex-extraction --coherent-scalar-extraction --CVE --CSE"
      write(*,*) "--compare-keys"
      write(*,*) "--compare-timeseries"
      write(*,*) "--convert-to-wing-system  // --convert-to-body-system"
      write(*,*) "--cp"
      write(*,*) "--crop"
      write(*,*) "--dist2mask   --dist2chi"
      write(*,*) "--divergence"
      write(*,*) "--energy"
      write(*,*) "--extract-subset"
      write(*,*) "--field-analysis"
      write(*,*) "--force-decomp"
      write(*,*) "--force"
      write(*,*) "--flexible-wing-mask"
      write(*,*) "--gradient"
      write(*,*) "--hdf2bin"
      write(*,*) "--helicity"
      write(*,*) "--keyvalues"
      write(*,*) "--magnitude"
      write(*,*) "--max-over-x"
      write(*,*) "--mean-2D"
      write(*,*) "--mean-over-x-subdomain"
      write(*,*) "--p2Q"
      write(*,*) "--POD"
      write(*,*) "--pointcloud2mask"
      write(*,*) "--pressure-force"
      write(*,*) "--pressure"
      write(*,*) "--readwrite"
      write(*,*) "--remove-mean"
      write(*,*) "--set-hdf5-attribute"
      write(*,*) "--simple-field-operation"
      write(*,*) "--smooth-inverse-mask"
      write(*,*) "--spectrum"
      write(*,*) "--stl2dist"
      write(*,*) "--time-avg"
      write(*,*) "--TKE-mean"
      write(*,*) "--tranpose-test"
      write(*,*) "--turbulence-analysis"
      write(*,*) "--upsample"
      write(*,*) "--uCT-assemble"
      write(*,*) "--ux-from-uyuz"
      write(*,*) "--vor_abs --vor-abs"
      write(*,*) "--vor2u"
      write(*,*) "--vorticity","--vor-abs","--vorticity-FD","--vor-abs-FD"
      write(*,*) "--Q   --Q-FD"
      write(*,*) "--zero-padd"
      write(*,'(80("~"))')
      write(*,*) "Postprocessing option is "// trim(adjustl(postprocessing_mode))
      write(*,*) "But I don't know what to do with that"
    endif
  end select

  if ((mpirank==0).and.(help.eqv..false.)) then
    write(*,'("Elapsed time=",es12.4)') MPI_wtime()-t1
  endif
end subroutine postprocessing




!-------------------------------------------------------------------------------
! write the call to flusi (i.e. the command line arguments) to an ascii file
! use this in conjunction with small ascii output files to document a little
! bit what you have been doing.
!-------------------------------------------------------------------------------
subroutine postprocessing_ascii_header( io_stream )
  use vars
  implicit none
  integer, intent(in) :: io_stream
  integer :: i
  character(len=strlen) :: arg

  if (mpirank /= 0) return

  ! MATLAB comment character:
  write(io_stream,'(A,1x)',advance='no') "% CALL: ./flusi "

  arg = "-p"
  i=1
  ! loop over command line arguments and dump them to the file:
  do while ( arg /= "" )
    call get_command_argument(i,arg)
    write(io_stream,'(A,1x)',advance='no') trim(adjustl(arg))
    i=i+1
  end do

  write(io_stream,'(A,1x)',advance='yes') "%"
end subroutine postprocessing_ascii_header



!-------------------------------------------------------------------------------
! ./flusi -p --max-over-x tkeavg_000.h5 outfile.dat
! This function reads in the specified *.h5 file and outputs the maximum value
! max_yz(x) into the specified ascii-outfile
! It may be used rarely, but we needed it for turbulent bumblebees.
! Can be done in parallel.
!-------------------------------------------------------------------------------
subroutine max_over_x(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers, only : check_file_exists
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ekin, outfile
  real(kind=pr),dimension(:,:,:),allocatable :: u
  real(kind=pr), dimension(:), allocatable :: umaxx,umaxx_loc
  real(kind=pr) :: time
  integer :: ix, mpicode

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --max-over-x tkeavg_000.h5 outfile.dat"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! This function reads in the specified *.h5 file and outputs the maximum value"
    write(*,*) "! max_yz(x) into the specified ascii-outfile"
    write(*,*) "! It may be used rarely, but we needed it for turbulent bumblebees."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,outfile)
  call check_file_exists( fname_ekin )

  call fetch_attributes( fname_ekin, nx, ny, nz, xl, yl, zl, time, nu, origin )

  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( umaxx(0:nx-1) )
  allocate( umaxx_loc(0:nx-1) )

  call read_single_file ( fname_ekin, u )

  do ix=0,nx-1
    umaxx_loc(ix)=maxval(u(ix,:,:))
  enddo

  call MPI_ALLREDUCE(umaxx_loc,umaxx,nx,MPI_DOUBLE_PRECISION,MPI_MAX,&
  MPI_COMM_WORLD,mpicode)

  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    write(17,'(A)') "%-----------------------------------"
    write(17,'(A)') "%FLUSI max-over-x file="//trim(adjustl(fname_ekin))
    write(17,'(A)') "%-----------------------------------"
    do ix=0,nx-1
      write(17,'(es15.8)') umaxx(ix)
    enddo
    close(17)
  endif



  deallocate (u,umaxx,umaxx_loc)
  call fft_free()

end subroutine max_over_x








!-------------------------------------------------------------------------------
! ./flusi -p --mean_over_x_subdomain tkeavg_000.h5 outfile.dat
! Compute the avg value as a function of x for a subdomain [-1.3,1.3]x[-1.3,1.3]
! in the y-z plane
! Can be done in parallel.
!-------------------------------------------------------------------------------
subroutine mean_over_x_subdomain(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers, only : check_file_exists
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ekin, outfile
  real(kind=pr),dimension(:,:,:),allocatable :: u, mask
  real(kind=pr), dimension(:), allocatable :: umaxx,umaxx_loc
  real(kind=pr) :: time,x,y,z,points,allpoints
  integer :: ix,iy,iz, mpicode

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --mean_over_x_subdomain tkeavg_000.h5 outfile.dat "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Compute the avg value as a function of x for a subdomain [-1.3,1.3]x[-1.3,1.3] in the y-z plane"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,outfile)
  call check_file_exists( fname_ekin )

  call fetch_attributes( fname_ekin, nx, ny, nz, xl, yl, zl, time, nu, origin )

  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( umaxx(0:nx-1) )
  allocate( umaxx_loc(0:nx-1) )

  call read_single_file ( fname_ekin, u )

  ! set a 1/0 mask to cancel values out of the bounds we want
  mask =0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x=dble(ix)*dx-2.0
        ! note we center y,z in the middle
        y=dble(iy)*dy-0.5*yl
        z=dble(iz)*dz-0.5*zl
        ! is the point inside the valid bounds?
        if ( abs(y)<=1.3d0 .and. abs(z)<=1.3d0 ) then
          mask(ix, iy, iz) = 1.d0
        endif
      enddo
    enddo
  enddo

  ! check if we loaded bullshit
  call checknan(mask,"mask")
  call checknan(u,"energy")

  if (minval(u)<0.d0) then
    write(*,*) "Warning, E<0 will produce NaN...correct that"
    where (u<0.d0)
      u=0.d0
    end where
  endif

  ! make TU intensity out of e
  u = dsqrt(2.d0*u/3.d0) * mask

  umaxx_loc=0.d0
  umaxx=0.d0

  call checknan(u,"urms")

  do ix=0,nx-1
    if (sum(mask(ix,:,:))>0.d0) then
      umaxx_loc(ix)=sum(u(ix,:,:))
    else
      umaxx_loc(ix)=0.0
    endif
  enddo

  ! get results from all CPUs, gather on all ranks
  call MPI_ALLREDUCE(umaxx_loc,umaxx,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

  ! get total number of points from all CPUs
  points = sum(mask)/dble(nx)
  call MPI_ALLREDUCE(points,allpoints,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)


  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    write(17,'(A)') "%-----------------------------------"
    write(17,'(A)') "%FLUSI mean over x (boxed) file="//trim(adjustl(fname_ekin))
    write(17,'(A)') "%-----------------------------------"
    do ix=0,nx-1
      write(17,'(es15.8)') umaxx(ix) / allpoints
    enddo
    close(17)
  endif



  deallocate (u,umaxx,umaxx_loc,mask)
  call fft_free()

end subroutine mean_over_x_subdomain



!-------------------------------------------------------------------------------
! ./flusi --postprocess --ux-from-uyuz ux_00000.h5 uy_00000.h5 uz_00000.h5
!-------------------------------------------------------------------------------
! compute missing ux component from given uy,uz components assuming incompressibility
subroutine ux_from_uyuz(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time
  integer :: ix,iy,iz
  real(kind=pr) :: kx,ky,kz

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --ux-from-uyuz ux_00000.h5 uy_00000.h5 uz_00000.h5 "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " not working"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)

  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )


  call fetch_attributes( fname_uy, nx, ny, nz, xl, yl, zl, time, nu, origin )

  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  u =0.d0
  uk = dcmplx(0.d0,0.d0)

  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft(u(:,:,:,2), uk(:,:,:,2))
  call fft(u(:,:,:,3), uk(:,:,:,3))

  call dealias(uk(:,:,:,2))
  call dealias(uk(:,:,:,3))

  do iz=ca(1),cb(1)
    !-- wavenumber in z-direction
    kz = wave_z(iz)
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)
        if (kx > 1.0d-12) then
          uk(iz,iy,ix,1) = -(ky/kx)*uk(iz,iy,ix,2) -(kz/kx)*uk(iz,iy,ix,3)
        else
          uk(iz,iy,ix,1) = dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo
  call dealias(uk(:,:,:,1))
  call ifft(uk(:,:,:,1), u(:,:,:,1))

  ! call fft (inx=u(:,:,:,1),outk=uk(:,:,:,1))
  ! call ifft (u(:,:,:,1),uk(:,:,:,1))

  call save_field_hdf5 ( time,fname_ux,u(:,:,:,1))
  if (mpirank==0) write(*,*) "Wrote vorx to "//trim(fname_ux)

  deallocate (u)
  deallocate (uk)
  call fft_free()

end subroutine


!-------------------------------------------------------------------------------
! ./flusi --postprocess --check-params-file PARAMS.ini
!-------------------------------------------------------------------------------
! load a parameter file and check for a bunch of common mistakes/typos
! you tend to make, in order to help preventing stupid mistakes
subroutine check_params_file(help)
  use vars
  use solid_model
  use module_insects
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: infile
  type(diptera) :: Insect

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --check-params-file PARAMS.ini template.ini"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Compare the given params.ini file with a reference template file."
    write(*,*) " If no reference file is given, we will perform a couple of tests only"
    write(*,*) " If a reference file is given, we check what entries are "
    write(*,*) " * present in both files"
    write(*,*) " * have been removed from template (thus deprecated)"
    write(*,*) " * have been added to template"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: no"
    return
  endif


  method="fsi"
  nf = 1
  nd = 3
  allocate(lin(1))

  call get_command_argument(3,infile)
  call check_file_exists(infile)
  call get_params(infile,Insect,.false.)

  call get_command_argument(4,infile)
  if (infile /= "") then
    ! we have a reference file given
    call check_file_exists(infile)
  endif


  ! now we have the parameters and perform the tests
  if ((dx /= dy).or.(dx /= dz).or.(dz /=dy)) then
    write(*,*) "The resolution is NOT equidistant ", dx,dy,dz
  endif

  if (iMask=="insect") then
    if ((Insect%BodyMotion=="free_flight").and.(iTimeMethodFluid/="AB2_rigid_solid")) then
      write(*,*) "Insect%BodyMotion==free_flight but iTimeMethodFluid="//trim(adjustl(iTimeMethodFluid))
    endif
  endif

  if ((method=="fsi").and.(use_slicing=="yes")) then
    if (maxval(slices_to_save)>nx-1) then
      write(*,*) "Slicing is ON but at least one index is out of bounds: ", slices_to_save
    endif
  endif

  if ((use_solid_model=="yes").and.(iMask /= "Flexibility")) then
    write(*,*) "we use the solid model but the mask is wrongly set"
  endif

  write(*,'("Penalization parameter C_eta=",es12.4," and K_eta=",es12.4)') eps, &
  sqrt(nu*eps)/dx

  write(*,'("This simulation will produce ",f5.1,"GB HDD output, if it runs until the end")') &
  dble(iSaveVelocity*3+iSavePress+iSaveVorticity*3+iSaveMask+iSaveSolidVelocity*3) &
  *((dble(nx)*dble(ny)*dble(nz))*4.d0/1000.d0**3)*tmax/tsave &
  +dble(iDoBackup*2)*18.d0*((dble(nx)*dble(ny)*dble(nz))*4.d0/1000.d0**3)*2.d0 !two backup files à 18 fields double precision each.


end subroutine



!-------------------------------------------------------------------------------
!./flusi -p --smooth-inverse-mask mask_000.h5 smooth_000.h5
!-------------------------------------------------------------------------------
subroutine post_smooth_mask(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_ini_files_parser_mpi
  use mpi
  use module_helpers, only : check_file_exists
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: infile, outfile
  complex(kind=pr),dimension(:,:,:),allocatable :: uk,rhs
  real(kind=pr),dimension(:,:,:),allocatable :: u1,u2,tmp
  real(kind=pr),dimension(:,:,:,:),allocatable :: expvis
  real(kind=pr) :: time, dt
  integer :: mpicode, it
  type(inifile) :: params

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --smooth-inverse-mask mask_000.h5 smooth_000.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Reads in the mask function and solves the penalized heat equation, i.e. it diffuses the"
    write(*,*) "mask while trying to keep 1 inside the solid. then, we save (1-mask) to disk, so that "
    write(*,*) "one can multiply a field with it in order to smoothly cut out the region near the "
    write(*,*) "penalized domain. The equation is:"
    write(*,*) " d chi / dt = nu*laplace(chi) - (chi0/eta)*(chi-1.0)"
    write(*,*) " where chi is the dilluted mask and chi0 is the initial mask (sharpened)"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "The code reads in a small smoothing.ini file which contains the required parameters"
    write(*,*) "contents of this file:"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "[smoothing]"
    write(*,*) "nu=1.0e-1; viscosity"
    write(*,*) "eps=1e-2; penalization"
    write(*,*) "tmax=0.2; final time"
    write(*,*) "dt=0.95e-2; time step"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  call get_command_argument(3,infile)
  call get_command_argument(4,outfile)

  call check_file_exists( infile )
  call fetch_attributes( infile, nx, ny, nz, xl, yl, zl, time, nu, origin )

  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  nf = 1

  ! read in parameters
  call read_ini_file_mpi( params, "smoothing.ini", .true.)
  call read_param_mpi( params,"smoothing","nu",nu,1.0d-1 )
  call read_param_mpi( params,"smoothing","eps",eps,1.0d-2 )
  call read_param_mpi( params,"smoothing","dt",dt,0.5d-2 )
  call read_param_mpi( params,"smoothing","tmax",tmax,1.d0 )
  call clean_ini_file_mpi( params )

  allocate(lin(1))
  lin(1) = nu

  call fft_initialize() ! also initializes the domain decomp

  allocate(u1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(u2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  allocate(rhs(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  allocate(expvis(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf))

  call read_single_file ( infile, u1 )

  !-----------------------------------------------------------------------------
  ! we use u2 for forcing:
  u2 = u1
  if (root) write(*,*) "Removing smoothing layer from source term..."
  ! wherever the mask is present, we will force to one, in order to exclude the
  ! smoothing layer. That way, we are hopeful that the resulting function will be
  ! zero at the first onset on smoothing
  where (u2>1.0d-10)
    u2 = 1.d0
  end where

  !-----------------------------------------------------------------------------
  ! the initial condition is what we read from file (stepping in Fourier space)
  call fft ( inx=u1, outk=uk )

  ! compute number of time steps (approx.) to reach tmax from ini file
  nt = nint(tmax/dt)
  if (root) then
    write(*,*) "nt=",nt
    write(*,'("Penalization parameter C_eta=",es12.4," and K_eta=",es12.4)') eps, &
    sqrt(nu*eps)/dx
  endif

  ! compute integrating factor for diffusive term
  call cal_vis( dt, expvis )

  !-----------------------------------------------------------------------------
  ! loop over time steps for penalized heat equation problem
  !-----------------------------------------------------------------------------
  do it=1, nt
    if (root) then
      write(*,*) "step",it,"of",nt
    endif

    ! this is the penalization term:
    tmp = -u2/eps*(u1-1.d0)
    ! to Fourier space
    call fft( inx=tmp, outk=rhs )
    ! advance in time (euler), viscosity is treated with integrating factor:
    uk = (uk + dt*rhs)*expvis(:,:,:,1)
    ! for the next time step we need phys. space to compute penalization term
    call ifft ( ink=uk, outx=u1)
  enddo

  ! save inverse of smoothed mask to disk
  tmp = 1.d0 - u1
  call save_field_hdf5 ( time, outfile, tmp )

  ! done
  deallocate (u1,u2,uk,expvis,rhs,tmp)
  call fft_free()

end subroutine post_smooth_mask



subroutine readwrite(help)
  use vars
  use p3dfft_wrapper
  use module_helpers, only : check_file_exists

  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --readwrite ux_00.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Read a file and directly write it to disk. Useful for converting precision."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes, sure"
    return
  endif

  call get_command_argument(3,fname)
  call check_file_exists( fname )
  call fetch_attributes( fname, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call read_single_file( fname, work )
  call save_field_hdf5( time, fname, work )
  deallocate(work)

end subroutine readwrite


subroutine remove_mean(help)
  use vars
  use p3dfft_wrapper
  use basic_operators, only : fieldmean
  use module_helpers, only : check_file_exists

  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname
  real(kind=pr),dimension(:,:,:),allocatable :: work
  real(kind=pr) :: time, mean

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --remove-mean ux_00.h5 output_00.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Read a file, remove the mean value from the field, save."
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes, sure"
    return
  endif

  call get_command_argument(3,fname)
  call check_file_exists( fname )
  call fetch_attributes( fname, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize domain decomposition, but do not use FFTs
  call decomposition_initialize()

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call read_single_file( fname, work )

  ! compute mean value
  mean  = fieldmean(work)
  ! remove it from datafield. it now has zero mean.
  work = work - mean
  call get_command_argument(4,fname)
  call save_field_hdf5( time, fname, work )
  deallocate(work)

end subroutine remove_mean
