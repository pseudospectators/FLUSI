! MHD wrapper for different (possibly time-dependend) mask functions 
subroutine Create_Mask_mhd(time)
  use mpi_header
  use mhd_vars
  implicit none

  real(kind=pr), intent(in) :: time

  ! Actual mask functions:
  select case (iMask)
  case default    
     if (mpirank == 0) then
        write (*,*) "iMask not properly set. Suicide!"
        stop
     endif
  end select
end subroutine Create_Mask_mhd
