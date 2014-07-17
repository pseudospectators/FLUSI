!-------------------------------------------------------------------------------
! runtime control routines
! flusi regularily reads from a file runtime_control.ini if it should do some-
! thing, such as abort, reload_params or save data.
!-------------------------------------------------------------------------------
subroutine Initialize_runtime_control_file()
  ! overwrites the file again with the standard runtime_control file
  use vars
  use mpi
  implicit none
  
  open  (14,file='runtime_control.ini',status='replace')
  write (14,'(A)') "# This is flusi's runtime control file"
  write (14,'(A)') "# You can influence a running simulation with this file"
  write (14,'(A)') "# flusi regularily reads it and checks what to do"
  write (14,'(A)') "# "
  write (14,'(A)') "# "
  write (14,'(A)') "# The following is currently available:"
  write (14,'(A)') "# "
  write (14,'(A)') "# ----------------------------------"
  write (14,'(A)') "# Re-read PARAMS.ini file from disk"
  write (14,'(A)') "# ----------------------------------"
  write (14,'(A)') "# This is useful to correct some parameters. "
  write (14,'(A)') "# Be careful: the entire file is reloaded. Some parameters"
  write (14,'(A)') "# can easily be changed, for example the penalization param"
  write (14,'(A)') "# eps or saving flags."
  write (14,'(A)') "# Others cause an immediate crash (eg changing nx,ny,nz)"
  write (14,'(A)') "# runtime_control=reload_params;"
  write (14,'(A)') "# "
  write (14,'(A)') "# ----------------------------------"
  write (14,'(A)') "# Safe stop (dump backup then stop)"
  write (14,'(A)') "# ----------------------------------"  
  write (14,'(A)') "# Stop the run but makes a backup first"
  write (14,'(A)') "# Memory is properly dealloacted, unlike KILL"
  write (14,'(A)') "# runtime_control=save_stop;"
  write (14,'(A)') ""
  write (14,'(A)') "[runtime_control]"
  write (14,'(A)') "runtime_control=nothing;"
  close (14)        
      
end subroutine Initialize_runtime_control_file



subroutine runtime_control_command( command )
  ! reads runtime control command
  use vars
  use mpi
  implicit none
  character(len=strlen), intent(out)     :: command
  character(len=strlen) :: file
  character(len=strlen), dimension(1:nlines) :: CTRL_FILE
  integer :: i

  file ="runtime_control.ini"
  ! root reads in the control file
  ! and fetched the command 
  if (mpirank ==0) then
    call read_params_file( CTRL_FILE, i, file, .false. ) ! false = non-verbose
  endif
  command = "nothing"
  call param_str(CTRL_FILE, i, "runtime_control","runtime_control", &
                        command, command)
      
  ! all MPI procs now know what to do.
  return
end subroutine runtime_control_command


