!-------------------------------------------------------------------------------
! pressure_to_Qcriterion()
! converts a given pressure field and outputs the Q-criterion, computed with
! second order and periodic boundary conditions. Alternatively, one
! may use distint postprocessing tools, such as paraview, and compute the Q-crit
! there, but the accuracy may be different.
! note the precision is reduced to second order (by using the effective
! wavenumber), since spurious oscillations appear when computing it with
! spectral precision
!-------------------------------------------------------------------------------
! call:
! ./flusi --postprocess --p2Q p_00000.h5 Q_00000.h5
!-------------------------------------------------------------------------------
subroutine pressure_to_Qcriterion(help)
  use mpi
  use vars
  use basic_operators
  use p3dfft_wrapper
  use module_helpers

  implicit none
  character(len=strlen) :: fname_p, fname_Q
  complex(kind=pr),dimension(:,:,:),allocatable :: pk
  real(kind=pr),dimension(:,:,:),allocatable :: p
  real(kind=pr)::time,maxi,mini
  logical, intent(in) :: help

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) ""
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif


  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_p)
  call check_file_exists( fname_p )
  ! get filename to save Q criterion to
  call get_command_argument(4,fname_Q)

  ! read in information from the file
  call fetch_attributes( fname_p, nx, ny, nz, xl, yl, zl, time, nu, origin )

  if (mpirank==0) then
    write(*,'("Computing Q criterion from  file ",A," saving to ",&
    & A," nx=",i4," ny=",i4," nz=",i4, &
    &"xl=",es12.4," yl=",es12.4," zl=",es12.4 )') &
    trim(fname_p), trim(fname_Q), nx,nx,nz,xl,yl,zl
  endif

  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  allocate(p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  call read_single_file(fname_p, p)
  call fft(inx=p, outk=pk)
  ! Q-criterion is 0.5*laplace(P)
  ! see http://books.google.de/books?id=FWmsrNv3BYoC&pg=PA23&lpg=PA23&dq=q-criterion&source=bl&ots=CW-1NPn9p0&sig=Zi_Z2iw-ZuDqJqYctM9OUrb5WMA&hl=de&sa=X&ei=ZBCXU9nHGInVPPTRgagO&ved=0CCgQ6AEwADgK#v=onepage&q=q-criterion&f=false
  call laplacien_inplace_filtered(pk)
  call ifft(ink=pk, outx=p)
  p=0.5d0*p

  call save_field_hdf5(time, fname_Q, p)

  maxi = fieldmax(P)
  mini = fieldmin(P)

  if (mpirank==0) then
    write(*,'("Q-criterion 2nd order maxval=",es12.4," minval=",es12.4)') maxi,mini
  endif

  deallocate(p,pk)
  call fft_free()
end subroutine pressure_to_Qcriterion
