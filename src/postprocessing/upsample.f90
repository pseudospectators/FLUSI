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
  use hdf5_wrapper
  use module_helpers
  use slicing
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname_in, fname_out, tmp
  integer :: nx_new, ny_new, nz_new
  integer :: nx_org, ny_org, nz_org, ix_org,iy_org,iz_org, ix_new,iy_new,iz_new
  integer :: i,j,k, mpicode, q
  real(kind=pr) :: time, kx_org,ky_org,kz_org, kx_new,ky_new,kz_new, E, E2
  complex(kind=pr),dimension(:,:,:),allocatable :: uk_org, uk_new, uk
  real(kind=pr),dimension(:,:,:),allocatable :: u_org, u_new
  integer, dimension(1:3) :: ra_org,rb_org,ca_org,cb_org,ra_new,rb_new,ca_new,cb_new

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --upsample source.h5 target.h5 256 256 526"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Upsampling from a source resolution to a target resolution by zero-padding in k-space"
    write(*,*) ""
    write(*,*) "The routine is half-parallel: we read the source field without parallelization, i.e. "
    write(*,*) "each proc contains the entire source field. This works ok if your source field is not"
    write(*,*) "too large. Then, the target field is allocated in a distributed fashion and the source"
    write(*,*) "data is zero padded. It is suprisingly difficult to zero-padd the data completely in parallel"
    write(*,*) "and we use upsampling so incredibly rarely (namely for initial conditions) that is it not"
    write(*,*) "a performance critical part of the code. optimizing it is not worth the pain, as long as"
    write(*,*) "the memory restriction is not limiting us."
    write(*,*) ""
    write(*,*) " Step 1: read data [PARALLEL]"
    write(*,*) " Step 2: FFT (small) data [PARALLEL]"
    write(*,*) " Step 3: gather & zero-padd [SERIAL]"
    write(*,*) " Step 4: iFFT (big) data [PARALLEL]"
    write(*,*) " Step 5: save (big) data [PARALLEL]"
    write(*,*) ""
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: in part"
    return
  endif


  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call get_command_argument(4,fname_out)
  call check_file_exists( fname_in )

  ! read target resolution from command line
  call get_command_argument(5,tmp)
  read (tmp,*) nx_new
  call get_command_argument(6,tmp)
  read (tmp,*) ny_new
  call get_command_argument(7,tmp)
  read (tmp,*) nz_new

  call fetch_attributes( fname_in, nx_org, ny_org, nz_org, xl, yl, zl, time, nu, origin )

  if (mpirank==0) then
    write(*,'("Target resolution= ",3(i4,1x))') nx_new, ny_new, nz_new
    write(*,'("Origin resolution= ",3(i4,1x))') nx_org, ny_org, nz_org
  endif


  !-----------------------
  if (mpirank==0) write(*,*) "Initializing small FFT and transforming source field to k-space"
  nx = nx_org;  ny = ny_org;  nz = nz_org
  call fft_initialize()
  ra_org = ra;  rb_org = rb
  ca_org = ca;  cb_org = cb
  !-----------------------
  ! array for input array, x and k-space
  allocate(u_org(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk_org(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  ! it is always a good idea to do this test
  call fft_unit_test(u_org, uk_org)

  ! read file from disk
  call read_single_file(fname_in, u_org)

  E = sum(u_org**2) / 2.d0
  E = mpisum(E) / (dble(nx_org)*dble(ny_org)*dble(nz_org))
  if (mpirank==0) write(*,*) "Input energy, x-space", E

  ! transform to fourier space
  call fft(inx=u_org, outk=uk_org)

  call compute_energies1_k(uk_org,E)
  if (mpirank==0) write(*,*) "Input energy, k-space", E / (xl*yl*zl)


  ! we don't need the original array anymore
  deallocate(u_org)

  ! array that hold the entire, non-mpi-distributed fourier transform of input data
  allocate(uk(minval(ca_table(1,:)):maxval(cb_table(1,:)),&
              minval(ca_table(2,:)):maxval(cb_table(2,:)),&
              minval(ca_table(3,:)):maxval(cb_table(3,:))))
  ! gather the entire 3D source (Fourier-transformed) field on root
  ! Note this limits parallelism im memory: This routine works only if you can fit
  ! the complete source data on each proc.
  call gather_all( uk_org, uk )

  ! now that we have a copy of the Fourier coefficients, delete the (small) part:
  deallocate( uk_org )

  ! and distribute it to all procs
  call MPI_BCAST(uk, size(uk,1)*size(uk,2)*size(uk,3), MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, mpicode)

  ! free fft, in order to prepare the larger one
  call fft_free()
  !-----------------------
  if (mpirank == 0) then
    write(*,*) "Initializing big FFT and copying source Fourier coefficients to target field in k-space"
  endif
  nx = nx_new;  ny = ny_new;  nz = nz_new
  call fft_initialize()
  ra_new = ra;  rb_new = rb
  ca_new = ca;  cb_new = cb
  !-----------------------


  allocate(uk_new(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  uk_new = dcmplx(0.d0,0.d0)

  !------------------------------------------------------------------
  do iz_org = lbound(uk,1), ubound(uk,1)
    nx=nx_org; ny=ny_org; nz=nz_org
    kz_org=wave_z(iz_org)
    ! find corresponding iz_new
    do iz_new=ca_new(1),cb_new(1)
      nx=nx_new; ny=ny_new; nz=nz_new
      kz_new = wave_z(iz_new)
      if ( abs(kz_new-kz_org) <= 1.0d-7 ) exit
    enddo
    ! we now have the pair (iz_org, iz_new)
    !------------------------------------------------------------------
    do iy_org = lbound(uk,2), ubound(uk,2)
      nx=nx_org; ny=ny_org; nz=nz_org
      ky_org=wave_y(iy_org)
      ! find corresponding iz_new
      do iy_new=ca_new(2),cb_new(2)
        nx=nx_new; ny=ny_new; nz=nz_new
        ky_new = wave_y(iy_new)
        if (abs(ky_new-ky_org) <= 1.0d-7) exit
      enddo
      ! we now have the pair (iy_org, iy_new)
      !------------------------------------------------------------------
      do ix_org = lbound(uk,3), ubound(uk,3)
        nx=nx_org; ny=ny_org; nz=nz_org
        kx_org=wave_x(ix_org)
        ! find corresponding iz_new
        do ix_new=ca_new(3),cb_new(3)
          nx=nx_new; ny=ny_new; nz=nz_new
          kx_new = wave_x(ix_new)
          if (abs(kx_new-kx_org) <= 1.0d-7) exit
        enddo
        ! we now have the pair (ix_org, ix_new)

        ! copy the old Fourier coefficients to the new field
        if ( iz_new >= ca_new(1) .and. iz_new <= cb_new(1)) then
          if ( iy_new >= ca_new(2) .and. iy_new <= cb_new(2)) then
            if ( ix_new >= ca_new(3) .and. ix_new <= cb_new(3)) then
              uk_new(iz_new,iy_new,ix_new) = uk(iz_org,iy_org,ix_org)
            endif
          endif
        endif

      enddo
    enddo
  enddo

  call compute_energies1_k(uk_new,E2)
  if (mpirank==0) write(*,*) "Output energy, k-space", E2 / (xl*yl*zl)

  ! as we have now zero-padded the original data, we can free that very large array:
  deallocate( uk )

  ! this will be our final result in x-space:
  allocate(u_new(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  ! transform the zero-padded Fourier coefficients back to physical space. this
  ! is the upsampled (=interpolated) field.
  if (mpirank==0) write(*,*) "transforming zero-padded Fourier coefficients back to x-space"
  call ifft(ink=uk_new,outx=u_new)

  E2 = sum(u_new**2) / 2.d0
  E2 = mpisum(E2) / (dble(nx_new)*dble(ny_new)*dble(nz_new))
  if (mpirank==0) write(*,*) "Output energy, x-space", E2
  call compute_energies1(u_new,E2)
  if (mpirank==0) write(*,*) "Output energy, x-space", E2 / (xl*yl*zl)

  if ( abs(E2-E) > 1.0e-7) then
    ! if this error occurs, it is likely a bug in gather_all. it happens if there
    ! is some slight load imbalancing (eg n=64 ncpu=3)
    ! I'm also not very sure what happens if nx /= ny /= nz
    ! call abort(123,"During upsampling, energy changed, which is not right.")
  endif

  call fft_free()
  deallocate( uk_new )

  ! save the final result to the specified file
  call save_field_hdf5(time, fname_out, u_new)

  deallocate( u_new )
end subroutine upsample
