subroutine io_test
  use vars
  use p3dfft_wrapper
  use module_insects

  implicit none

  real(kind=pr),dimension(:,:,:),allocatable :: work
  type(diptera) :: Insect
  character(len=strlen)  :: infile
  integer :: ix,iy,iz

  ! Set method information in vars module.
  method="fsi" ! We are doing fluid-structure interactions
  nf=1    ! We are evolving one field (that means 1 integrating factor)
  nd=3*nf ! The one field has three components.
  neq=nd  ! number of equations, can be higher than 3 if using passive scalar

  !-----------------------------------------------------------------------------
  ! Read input parameters
  !-----------------------------------------------------------------------------
  allocate(lin(nf)) ! Set up the linear term
  if (root) write(*,'(A)') '*** info: Reading input data...'
  ! get filename of PARAMS file from command line
  call get_command_argument(2,infile)
  ! read all parameters from that file
  call get_params(infile,Insect,.true.)

  call decomposition_initialize()

  allocate(work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  do iz=ra(3), rb(3)
    do iy=ra(2), rb(2)
       do ix=ra(1), rb(1)
          work(ix,iy,iz)=rand_nbr()
       end do
    end do
  end do


  call save_field_hdf5(0.d0,'test_00.h5',work)
  call read_single_file('test_00.h5',work)

  call save_field_hdf5(0.d0,'test_00.h5',work)
  call read_single_file('test_00.h5',work)

  call save_field_hdf5(0.d0,'test_00.h5',work)
  call read_single_file('test_00.h5',work)

  deallocate(work)

end subroutine
