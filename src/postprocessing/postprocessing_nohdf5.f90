!-------------------------------------------------------------------------------
! This is a stub. Postprocessing tools are not available without HDF5 support.
!-------------------------------------------------------------------------------
subroutine postprocessing()
  use vars
  call abort(76181,'Postprocessing without HDF5 not implemented.')
end subroutine
