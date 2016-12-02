!-------------------------------------------------------------------------------
! Upsampling from a source resolution to a target resolution
! ./flusi -p --upsample source.h5 target.h5 256 256 526
!-------------------------------------------------------------------------------
! We first read in the original field from the source file, with it's resolution
! and domain size and timestamp.
!-------------------------------------------------------------------------------
subroutine upsample(help)
  use vars
  use p3dfft_wrapper
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_in, fname_out, tmp
  integer :: nx_new, ny_new, nz_new
  integer :: nx_org, ny_org, nz_org, ix_org,iy_org,iz_org, ix_new,iy_new,iz_new
  integer :: i,j,k
  real(kind=pr) :: time, kx_org,ky_org,kz_org, kx_new,ky_new,kz_new
  complex(kind=pr),dimension(:,:,:),allocatable :: uk_org, uk_new
  real(kind=pr),dimension(:,:,:),allocatable :: u_org, u_new
  integer, dimension(1:3) :: ra_org,rb_org,ca_org,cb_org,ra_new,rb_new,ca_new,cb_new

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --upsample source.h5 target.h5 256 256 526"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Upsampling from a source resolution to a target resolution"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: no"
    return
  endif



  if (mpisize/=1) then
    write(*,*) "./flusi --postprocess --upsample is a SERIAL routine, use 1CPU only"
    call abort()
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  call get_command_argument(4,fname_out)

  ! read target resolution from command line
  call get_command_argument(5,tmp)
  read (tmp,*) nx_new
  call get_command_argument(6,tmp)
  read (tmp,*) ny_new
  call get_command_argument(7,tmp)
  read (tmp,*) nz_new

  write(*,'("Target resolution= ",3(i4,1x))') nx_new, ny_new, nz_new

  call fetch_attributes( fname_in, nx_org, ny_org, nz_org, xl, yl, zl, time, nu )
  write(*,'("Origin resolution= ",3(i4,1x))') nx_org,ny_org,nz_org

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl


  !-----------------------
  write(*,*) "Initializing small FFT and transforming source field to k-space"
  nx = nx_org
  ny = ny_org
  nz = nz_org
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  call fft_initialize
  ra_org = ra
  rb_org = rb
  ca_org = ca
  cb_org = cb
  !-----------------------

  allocate(u_org(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk_org(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  call fft_unit_test( u_org, uk_org )

  write(*,*) "Reading file "//trim(adjustl(fname_in))
  call read_single_file(fname_in,u_org)

  call fft(inx=u_org, outk=uk_org)

  deallocate(u_org)

  call fft_free
  !-----------------------
  write(*,*) "Initializing big FFT and copying source Fourier coefficients to &
  & target field in k-space"

  nx=nx_new
  ny=ny_new
  nz=nz_new
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  call fft_initialize
  ra_new = ra
  rb_new = rb
  ca_new = ca
  cb_new = cb
  !-----------------------

  allocate(u_new(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk_new(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  call fft_unit_test( u_new, uk_new )
  uk_new = dcmplx(0.d0,0.d0)

  !------------------------------------------------------------------
  do iz_org=ca_org(1),cb_org(1)
    nx=nx_org; ny=ny_org; nz=nz_org
    kz_org=wave_z(iz_org)
    ! find corresponding iz_new
    do iz_new=ca_new(1),cb_new(1)
      nx=nx_new; ny=ny_new; nz=nz_new
      kz_new = wave_z(iz_new)
      if (kz_new==kz_org) exit
    enddo
    ! we now have the pair (iz_org, iz_new)
    !------------------------------------------------------------------
    do iy_org=ca_org(2),cb_org(2)
      nx=nx_org; ny=ny_org; nz=nz_org
      ky_org=wave_y(iy_org)
      ! find corresponding iz_new
      do iy_new=ca_new(2),cb_new(2)
        nx=nx_new; ny=ny_new; nz=nz_new
        ky_new = wave_y(iy_new)
        if (ky_new==ky_org) exit
      enddo
      ! we now have the pair (iy_org, iy_new)
      !------------------------------------------------------------------
      do ix_org=ca_org(3),cb_org(3)
        nx=nx_org; ny=ny_org; nz=nz_org
        kx_org=wave_x(ix_org)
        ! find corresponding iz_new
        do ix_new=ca_new(3),cb_new(3)
          nx=nx_new; ny=ny_new; nz=nz_new
          kx_new = wave_x(ix_new)
          if (kx_new==kx_org) exit
        enddo
        ! we now have the pair (ix_org, ix_new)

        ! copy the old Fourier coefficients to the new field
        uk_new(iz_new,iy_new,ix_new) = uk_org(iz_org,iy_org,ix_org)
      enddo
    enddo
  enddo

  deallocate( uk_org )

  ! transform the zero-padded Fourier coefficients back to physical space. this
  ! is the upsampled (=interpolated) field.
  write(*,*) "transforming zero-padded Fourier coefficients back to x-space"
  call ifft(ink=uk_new,outx=u_new)

  deallocate( uk_new )

  ! save the final result to the specified file
  write(*,*) "Saving upsampled field to " // trim(adjustl(fname_out))
  call save_field_hdf5(time,fname_out,u_new)

  deallocate( u_new )
end subroutine upsample
