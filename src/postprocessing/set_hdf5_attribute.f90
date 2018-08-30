!-------------------------------------------------------------------------------
! Add or modifiy an attribute to the dataset stored in a HDF5 file.
! This is useful for example if one creates stroke-averaged fields
! that would all have the same time 0.0.
! We also want to use this to add new attributes to our data, namely the
! viscosity
!-------------------------------------------------------------------------------
! ./flusi -p --set-hdf5-attrbute [FILE] [ATTRIBUTE_NAME] [N] [ATTRIBUTE_VALUE(S)]
!-------------------------------------------------------------------------------
subroutine set_hdf5_attribute(help)
  use mpi
  use vars
  use hdf5_wrapper
  use module_helpers
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, attribute_name, dsetname,tmp
  real(kind=pr),dimension(:),allocatable :: new_value

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --set-hdf5-attribute [FILE] [ATTRIBUTE_NAME] [N] [ATTRIBUTE_VALUE(S)]"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! Add or modifiy an attribute to the dataset stored in a HDF5 file."
    write(*,*) "! This is useful for example if one creates stroke-averaged fields"
    write(*,*) "! that would all have the same time 0.0."
    write(*,*) "! We also want to use this to add new attributes to our data, namely the"
    write(*,*) "! viscosity"
    write(*,*) "! N: length, can be 1 for scalar attributes"
    write(*,*) "! "
    write(*,*) "! example: (3 component attribute)"
    write(*,*) "! ./flusi -p --set-hdf5-attribute ux_0000.h5 origin 3 2.0,1.0,84.9"
    write(*,*) "! "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: no"
    return
  endif

  call get_command_argument(3,fname)
  call check_file_exists( fname )
  dsetname = get_dsetname(fname)

  call get_command_argument(4,attribute_name)
  call get_command_argument(5,tmp)
  read(tmp,*) nx
  allocate(new_value(1:nx))

  call get_command_argument(6,tmp)
  read(tmp,*) new_value

  call write_attribute( fname, dsetname, attribute_name, new_value )
  deallocate(new_value)
end subroutine set_hdf5_attribute
