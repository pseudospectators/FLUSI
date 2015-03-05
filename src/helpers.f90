module helpers
  use fsi_vars
  implicit none

  ! module global variables


contains

  !-----------------------------------------------------------------------------
  ! This function returns, to a given filename, the corresponding dataset name
  ! in the hdf5 file, following flusi conventions (folder/ux_0000.h5 -> "ux")
  !-----------------------------------------------------------------------------
  character(len=strlen)  function get_dsetname(fname)
    implicit none
    character(len=*), intent(in) :: fname
    ! extract dsetname (from "/" until "_", excluding both)
    get_dsetname  = fname  ( index(fname,'/',.true.)+1:index( fname, '_' )-1 )
    return
  end function get_dsetname

end module helpers
