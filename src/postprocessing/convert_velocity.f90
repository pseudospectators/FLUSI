!-------------------------------------------------------------------------------
! Read in vorticity fields and compute the corresponding velocity (Biot-Savart Law)
! Some checks are performed on the result, namely divergence(u) and we see
! if the curl(result) is the same as the input values.
!
! you can optionally specify the list of outfiles:
! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5 --outfiles outx_00.h5 outy_00.h5 outz_00.h5
!
! or write to the default file names: (ux_00.h5 in the example)
! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5
!
! or specify the basename prefix: (writes to outx_00.h5 in the example:)
! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5 --outputprefix out
!
!-------------------------------------------------------------------------------
subroutine convert_velocity(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, outfiles_given
  character(len=strlen) :: fname_outx, fname_outy, fname_outz, prefix
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk, workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u,workr
  real(kind=pr) :: time, divu_max, errx, erry, errz
  integer :: ix,iy,iz, mpicode

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --vor2u"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! Read in vorticity fields and compute the corresponding velocity (Biot-Savart Law)"
    write(*,*) "! Some checks are performed on the result, namely divergence(u) and we see"
    write(*,*) "! if the curl(result) is the same as the input values."
    write(*,*) "!"
    write(*,*) "! you can optionally specify the list of outfiles:"
    write(*,*) "! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5 --outfiles outx_00.h5 outy_00.h5 outz_00.h5"
    write(*,*) "!"
    write(*,*) "! or write to the default file names: (ux_00.h5 in the example)"
    write(*,*) "! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5"
    write(*,*) "!"
    write(*,*) "! or specify the basename prefix: (writes to outx_00.h5 in the example:)"
    write(*,*) "! ./flusi -p --vor2u vorx_00.h5 vory_00.h5 vorz_00.h5 --outputprefix out"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,outfiles_given)

  if ( outfiles_given == "--outfiles" ) then
    call get_command_argument(7,fname_outx)
    call get_command_argument(8,fname_outy)
    call get_command_argument(9,fname_outz)
  elseif ( outfiles_given == "--outputprefix" ) then
    ! if not specified, use standard output names:
    call get_command_argument(7,prefix)
    fname_outx=trim(adjustl(prefix))//"x"//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
    fname_outy=trim(adjustl(prefix))//"y"//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
    fname_outz=trim(adjustl(prefix))//"z"//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)
  else
    ! if not specified, use standard output names:
    fname_outx='ux'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
    fname_outy='uy'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
    fname_outz='uz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)
  endif

  ! header and information
  if (mpirank==0) then
    write(*,'(80("-"))')
    write(*,*) "vor2u (Biot-Savart)"
    write(*,*) "Computing velocity from vorticity given in these files: "
    write(*,'(80("-"))')
    write(*,*) trim(adjustl(fname_ux))
    write(*,*) trim(adjustl(fname_uy))
    write(*,*) trim(adjustl(fname_uz))
    write(*,*) "Writing to:"
    write(*,*) trim(adjustl(fname_outx))
    write(*,*) trim(adjustl(fname_outy))
    write(*,*) trim(adjustl(fname_outz))
    write(*,'(80("-"))')
    if ((fname_ux(1:4).ne."vorx").or.(fname_uy(1:4).ne."vory").or.(fname_uz(1:4).ne."vorz")) then
      write (*,*) "WARNING in arguments, files do not start with vorx vory and vorz"
      write (*,*) "note files have to be in the right order"
    endif
  endif

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  call fetch_attributes( fname_ux, nx, ny, nz, xl, yl, zl, time, nu, origin )

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  neq=3
  nd=3
  ncw=3
  nrw=3

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  ! read vorticity from files to u
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  ! compute divergence of input fields and show the maximum value, this is interesting
  ! to know since actually, div(curl(u)) should be identically zero. however, in
  ! a discrete setting, this may not strictly be the case
  call fft3( inx=u,outk=uk ) ! uk is now vork
  call divergence( uk, workc(:,:,:,1) )
  call ifft( ink=workc(:,:,:,1), outx=workr(:,:,:,1) )
  divu_max = fieldmax(workr(:,:,:,1))
  if(mpirank==0) write(*,'("maximum divergence in input field=",es12.4)') divu_max

  ! copy original vorticity to workr (for comparison later, error checks)
  workr = u

  !-----------------------------------------------------------------------------
  ! compute velocity from vorticity
  !-----------------------------------------------------------------------------
  ! uk: vorticity; workc: velocity
  call Vorticity2Velocity(uk(:,:,:,1:3),workc(:,:,:,1:3))
  call ifft3 (ink=workc, outx=u)
  Uxmean=0.0_pr
  Uymean=0.0_pr
  Uzmean=0.0_pr
  call set_mean_flow(workc, time)
  uk = workc ! uk: velocity in F-space

  ! now u contains the velocity in physical space
  call save_field_hdf5 ( time,fname_outx,u(:,:,:,1) )
  call save_field_hdf5 ( time,fname_outy,u(:,:,:,2) )
  call save_field_hdf5 ( time,fname_outz,u(:,:,:,3) )

  if (mpirank==0) then
    write(*,*) "Done writing output!"
    write(*,'(80("-"))')
    write(*,*) "Done computing, performing some analysis of the result..."
  endif

  !-----------------------------------------------------------------------------
  ! check divergence of new field
  !-----------------------------------------------------------------------------
  call divergence(uk,workc(:,:,:,1))
  call ifft(ink=workc(:,:,:,1),outx=u(:,:,:,1))
  divu_max = fieldmax(u(:,:,:,1))
  if (mpirank==0) write(*,*) "max(div(u))", divu_max
  divu_max = fieldmin(u(:,:,:,1))
  if (mpirank==0) write(*,*) "min(div(u))", divu_max

  deallocate(workc)

  !-----------------------------------------------------------------------------
  ! check if the curl of computed velocity is the same as input values
  !-----------------------------------------------------------------------------
  call curl3_inplace(uk)
  call ifft3(ink=uk, outx=u)

  ! difference
  u(:,:,:,1)=u(:,:,:,1)-workr(:,:,:,1)
  u(:,:,:,2)=u(:,:,:,2)-workr(:,:,:,2)
  u(:,:,:,3)=u(:,:,:,3)-workr(:,:,:,3)

  errx = fieldmax(u(:,:,:,1))
  erry = fieldmax(u(:,:,:,2))
  errz = fieldmax(u(:,:,:,3))

  if (mpirank==0) then
    write(*,*) "max diff vor_original-curl(result)", errx
    write(*,*) "max diff vor_original-curl(result)", erry
    write(*,*) "max diff vor_original-curl(result)", errz
  endif

  ! check relative difference in mag(vor):
  errx=0.d0
  erry=0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        errx = errx + sqrt(u(ix,iy,iz,1)**2 + u(ix,iy,iz,2)**2 + u(ix,iy,iz,3)**2)
        erry = erry + sqrt(workr(ix,iy,iz,1)**2 + workr(ix,iy,iz,2)**2 + workr(ix,iy,iz,3)**2)
      enddo
    enddo
  enddo

  call MPI_ALLREDUCE (errx,errz,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE (erry,divu_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode)
  ! relative error:
  if (mpirank==0) write(*,*) "relative error in magnitude:", errz/divu_max

  deallocate (u)
  deallocate (uk)
  deallocate (workr)
  call fft_free()

end subroutine convert_velocity
