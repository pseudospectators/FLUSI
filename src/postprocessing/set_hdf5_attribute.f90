!-------------------------------------------------------------------------------
! Add or modifiy an attribute to the dataset stored in a HDF5 file.
! This is useful for example if one creates stroke-averaged fields
! that would all have the same time 0.0.
! We also want to use this to add new attributes to our data, namely the
! viscosity
!-------------------------------------------------------------------------------
! ./flusi -p --set-hdf5-attrbute [FILE] [ATTRIBUTE_NAME] [ATTRIBUTE_VALUE(S)]
!-------------------------------------------------------------------------------
subroutine set_hdf5_attribute(help)
  use mpi
  use hdf5
  use vars
  use helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, attribute_name, dsetname,tmp
  integer :: error

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: attr_id       ! attribute identifier
  integer(hid_t) :: aspace_id     ! attribute dataspace identifier
  integer(hsize_t), dimension(1) :: data_dims
  real(kind=pr) ::  attr_data  ! attribute data
  real(kind=pr) :: new_value

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --set-hdf5-attrbute [FILE] [ATTRIBUTE_NAME] [ATTRIBUTE_VALUE(S)]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! Add or modifiy an attribute to the dataset stored in a HDF5 file."
    write(*,*) "! This is useful for example if one creates stroke-averaged fields"
    write(*,*) "! that would all have the same time 0.0."
    write(*,*) "! We also want to use this to add new attributes to our data, namely the"
    write(*,*) "! viscosity"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: no"
    return
  endif


  data_dims(1) = 1

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname)
  call check_file_exists( fname )
  dsetname = get_dsetname(fname)

  ! get filename to save file to
  call get_command_argument(4,attribute_name)
  call get_command_argument(5,tmp)
  read(tmp,*) new_value

  ! Initialize FORTRAN interface.
  CALL h5open_f(error)
  ! Open an existing file.
  CALL h5fopen_f (fname, H5F_ACC_RDWR_F, file_id, error)
  ! Open an existing dataset.
  CALL h5dopen_f(file_id, dsetname, dset_id, error)

  if(mpirank==0) then
    write(*,'(80("-"))')
    write (*,*) "FLUSI tries to open the attribute in the file..."
  endif

  ! try to open the attribute
  CALL h5aopen_f(dset_id, trim(adjustl(attribute_name)), attr_id, error)

  if (error == 0) then
    if(mpirank==0) then
      write (*,*) "succesful! The file already contains the attribute"
      write(*,'(80("-"))')
    endif
    ! the attribute is already a part of the HDF5 file
    ! Get dataspace and read
    CALL h5aget_space_f(attr_id, aspace_id, error)
    ! read
    CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)
    if(mpirank==0) then
      write(*,'("old value of ",A," is ",es15.8)') trim(adjustl(attribute_name)), attr_data
    endif
    ! write
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, (/new_value/), data_dims, error)
    ! Close the attribute.
    CALL h5aclose_f(attr_id, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(aspace_id, error)
  else
    ! the attribute is NOT a part of the file yet, we add it now
    if(mpirank==0) then
      write(*,*) "fail: the file did not yet contain this attribute, adding it!"
      write(*,'(80("-"))')
    endif
    call write_attribute_dble(data_dims,trim(adjustl(attribute_name)),(/new_value/),1,dset_id)
  endif

  ! check if everthing worked
  if(mpirank==0) write(*,*) "checking if the operation worked..."
  ! try to open the attribute
  CALL h5aopen_f(dset_id, trim(adjustl(attribute_name)), attr_id, error)
  CALL h5aread_f( attr_id, H5T_NATIVE_DOUBLE, attr_data, data_dims, error)

  if(mpirank==0) write(*,'(A,"=",g12.4)') trim(adjustl(attribute_name)), attr_data

  ! finalize HDF5
  CALL h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
  CALL h5fclose_f(file_id, error) ! Close the file.
  CALL h5close_f(error)  ! Close FORTRAN interface.

end subroutine set_hdf5_attribute
