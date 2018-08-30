!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
subroutine transpose_test(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use hdf5_wrapper
  use mpi
  use flusi_wavelet_lib
  implicit none
  logical, intent(in) :: help
  integer, dimension(1:3) :: ka,kb,ks, kat, kbt, kst
  integer :: ix, iy, iz, idir, mpicode
  real(kind=pr) :: x,y,z
  real(kind=pr), dimension(:,:,:,:), allocatable :: work1, work2, work3
  type(orth_wavelet) :: wavelet

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --transpose-test"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! "
    write(*,*) "! "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  !-----------------------------------------------------------------------------
  ! Initializations
  !-----------------------------------------------------------------------------
  nx = 256; ny = 256; nz = 256
  xl = 2.d0*3.141592653589793d+00; yl=xl; zl=xl

  ! initialize code and domain decomposition, but do not use FFTs
  call decomposition_initialize()
  ! call fft_initialize()
  call setup_cart_groups()

  ! ============================================================================
  ! test transpose y
  ! ============================================================================

  allocate( work1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )
  allocate( work2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )
  allocate( work3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3) )

  work1 = 0.d0
  work2 = 0.d0
  work3 = 0.d0

  ! do iz=ra(3),rb(3)
  !   do iy=ra(2),rb(2)
  !     do ix=ra(1),rb(1)
  !       x = dble(ix)*dx
  !       y = dble(iy)*dy
  !       z = dble(iz)*dz
  !       ! work1(ix,iy,iz) = dsin((3.d0*x)**2)! 3.d0*sin(x) + sin(5.d0*x) +0.7*sin(10.d0*x)
  !       ! work1(ix,iy,iz) = 10.d0*dsin(x)*dsin(y)*dsin(z) + dsin(20.d0*x)*dsin(20.d0*y)*dsin(20.d0*z) + rand_nbr()
  !       work1(ix,iy,iz) = 10.d0*dsin(x)*dsin(y)*dsin(z) + rand_nbr()
  !     enddo
  !   enddo
  ! enddo

  ! call read_single_file( 'mask_00.h5', work1)
  !
  call read_single_file( 'vorx_0000.h5', work1(:,:,:,1))
  call read_single_file( 'vory_0000.h5', work1(:,:,:,2))
  call read_single_file( 'vorz_0000.h5', work1(:,:,:,3))


  work2 = work1
  ! call save_field_hdf5(0.d0,'utot_00.h5',work1)

  call setup_coiflet_coefs( 1, wavelet )

  ! call coherent_scalar_extraction( wavelet, 10, work1, work2, work3 )
  ! call save_field_hdf5(0.d0,'ucoh_00.h5',work2)
  ! call save_field_hdf5(0.d0,'uinc_00.h5',work3)

  call coherent_vortex_extraction( wavelet, 20, work1, work2, work3 )


  call read_single_file( 'vorcohx_0000.h5', work2(:,:,:,3))
  call read_single_file( 'vorcohy_0000.h5', work2(:,:,:,3))
  call read_single_file( 'vorcohz_0000.h5', work2(:,:,:,3))

  ! call FWT3_PO( work1, work3, wavelet, 1)
  ! call save_field_hdf5(0.d0,'wc_00.h5',work3)
  !
  !
  ! call IWT3_PO( work3, work1, wavelet, 1)
  ! call save_field_hdf5(0.d0,'IWT_00.h5',work1)
  ! call save_field_hdf5(0.d0,'error_00.h5',work2-work1)
  !   write(*,*) "error:", maxval(work2-work1)

  !
  !
  ! work3 = work1
  !
  ! ! idir = 1 - permute 1 and 2 indices
  ! ! idir = 2 - permute 1 and 3 indices
  ! idir = 1
  !
  ! call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  ! ! the indices ka, kb are the same as ra-ra=0 and rb-ra, so simply ka=0 always
  ! ! then kb is number of elements-1
  !
  ! allocate( work2(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
  ! work2 = 0.d0
  !
  ! call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, MPI_DOUBLE_PRECISION, work1, work2 )
  ! write(*,'(i3,32(f5.2,1x))') mpirank, work2(:,kat(2),kat(3))
  !
  ! call trextents(idir, (/ny,nx,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  ! call subtr ( idir, (/ny,nx,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, MPI_DOUBLE_PRECISION, work2, work1 )
  !
  ! call MPI_barrier(MPI_COMM_WORLD, mpicode)
  !
  ! write(*,*) "round-trip-error", maxval(work3-work1), mpirank
  ! ! call write_field_hdf5('after_00.h5', 'after', ra, rb, work1)
  ! ! call write_attribute( 'after_00.h5', 'after', "time",(/0.d0/))
  ! ! call write_attribute( 'after_00.h5', 'after', "domain_size",(/xl,yl,zl/))
  ! ! call write_attribute( 'after_00.h5', 'after', "nxyz",(/nx,ny,nz/))
  !
  ! deallocate( work1, work2, work3 )
  !
  ! call MPI_barrier(MPI_COMM_WORLD, mpicode)
  ! ! ============================================================================
  ! ! test transpose z
  ! ! ============================================================================
  !
  ! allocate( work1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  !
  ! do iz=ra(3),rb(3)
  !   do iy=ra(2),rb(2)
  !     do ix=ra(1),rb(1)
  !       x = dble(ix)*dx
  !       y = dble(iy)*dy
  !       z = dble(iz)*dz
  !       work1(ix,iy,iz) = z
  !     enddo
  !   enddo
  ! enddo
  !
  ! allocate( work3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  ! work3 = work1
  !
  ! ! idir = 1 - permute 1 and 2 indices
  ! ! idir = 2 - permute 1 and 3 indices
  ! idir = 2
  !
  ! call trextents(idir, (/nx,ny,nz/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  ! ! the indices ka, kb are the same as ra-ra=0 and rb-ra, so simply ka=0 always
  ! ! then kb is number of elements-1
  !
  ! allocate( work2(kat(1):kbt(1),kat(2):kbt(2),kat(3):kbt(3)) )
  ! work2 = 0.d0
  !
  ! call subtr ( idir, (/nx,ny,nz/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, MPI_DOUBLE_PRECISION, work1, work2 )
  !
  ! write(*,'(i3,32(f5.2,1x))') mpirank, work2(:,kat(2),kat(3))
  !
  ! call trextents(idir, (/nz,ny,nx/), mpidims, mpicoords, ka, kb, ks, kat, kbt, kst )
  ! call subtr ( idir, (/nz,ny,nx/), ka, kb, ks, kat, kbt, kst, mpidims, mpicoords, mpicommslab, MPI_DOUBLE_PRECISION, work2, work1 )
  !
  ! call MPI_barrier(MPI_COMM_WORLD, mpicode)
  !
  ! write(*,*) "round-trip-error", maxval(work3-work1), mpirank
  deallocate( work1, work2, work3 )
end subroutine transpose_test
