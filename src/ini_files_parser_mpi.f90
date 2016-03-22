! Ini files parser for MPI parallel codes
! This module is just a layer that includes the serial ini_files_parser
! if you read a parameter, only root rank extracts it from the file (also onl
! read by root) and then broadcasts it to all other processes
module ini_files_parser_mpi
  ! use the serial ini files parser. this module is just a wrapper for the
  ! mpi versions of our codes
  use ini_files_parser
  use mpi
  use vars, only : root, mpirank, pr, strlen


  ! the generic call "read_param" redirects to these routines, depending on the data
  ! type and the dimensionality. vectors can be read without setting a default.
  interface read_param_mpi
    module procedure param_dbl_mpi, param_int_mpi, param_vct_mpi, param_str_mpi, param_vct_nodefault_mpi
  end interface


!!!!!!!!
contains
!!!!!!!!


  !-------------------------------------------------------------------------------
  ! clean a previously read ini file, deallocate its string array, and reset
  ! verbosity to .true. (as a matter of precaution)
  !-------------------------------------------------------------------------------
  subroutine clean_ini_file_mpi(PARAMS)
    implicit none
    type(inifile), intent(inout) :: PARAMS

    if (root) then
      call clean_ini_file(PARAMS)
    endif
  end subroutine clean_ini_file_mpi

  !-------------------------------------------------------------------------------
  ! Read the file paramsfile, count the lines and put the
  ! text in PARAMS.
  !-------------------------------------------------------------------------------
  subroutine read_ini_file_mpi(PARAMS, file, verbose)
    implicit none

    type(inifile), intent(inout) :: PARAMS
    character(len=*) :: file ! this is the file we read the PARAMS from
    logical, intent(in) :: verbose

    if (root) then
      ! check if the specified file exists
      call check_file_exists( file )
      call read_ini_file( PARAMS, file, verbose )
    endif
  end subroutine read_ini_file_mpi




  !-------------------------------------------------------------------------------
  ! Fetches a REAL VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find the parameter, we return this and warn
  ! Output:
  !       params_real: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_dbl_mpi (PARAMS, section, keyword, params_real, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real (kind=pr) :: params_real, defaultvalue
    integer :: mpicode

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call read_param (PARAMS, section, keyword, params_real, defaultvalue)
    endif

    ! And then broadcast
    call MPI_BCAST( params_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_dbl_mpi




  !-------------------------------------------------------------------------------
  ! Fetches a STRING VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find the parameter, we return this and warn
  ! Output:
  !       params_string: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_str_mpi (PARAMS, section, keyword, params_string, defaultvalue)
    implicit none

    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! what section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=strlen), intent (inout) :: params_string
    character(len=*), intent (in) :: defaultvalue
    integer :: mpicode

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call read_param (PARAMS, section, keyword, params_string, defaultvalue)
    endif

    ! And then broadcast
    call MPI_BCAST( params_string, strlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpicode)
  end subroutine param_str_mpi



  !-------------------------------------------------------------------------------
  ! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find a vector, we return this and warn
  !       n: length of vector
  ! Output:
  !       params_vector: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_vct_mpi (PARAMS, section, keyword, params_vector, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=pr), intent(inout) :: params_vector(1:)
    real(kind=pr), intent(in) :: defaultvalue(1:)

    integer :: n
    integer :: mpicode

    n = size(params_vector,1)

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call read_param (PARAMS, section, keyword, params_vector, defaultvalue)
    endif

    ! And then broadcast
    call MPI_BCAST( params_vector, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_vct_mpi


  !-------------------------------------------------------------------------------
  ! Fetches a VECTOR VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find a vector, we return this and warn
  !       n: length of vector
  ! Output:
  !       params_vector: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_vct_nodefault_mpi (PARAMS, section, keyword, params_vector)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=pr) :: params_vector(1:)

    integer :: n
    integer :: mpicode

    n = size(params_vector,1)

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call read_param (PARAMS, section, keyword, params_vector)
    endif

    ! And then broadcast
    call MPI_BCAST( params_vector, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_vct_nodefault_mpi



  !-------------------------------------------------------------------------------
  ! Fetches a INTEGER VALUED parameter from the PARAMS.ini file.
  ! Displays what it does on stdout (so you can see whats going on)
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  !       defaultvalue: if the we can't find the parameter, we return this and warn
  ! Output:
  !       params_int: this is the parameter you were looking for
  !-------------------------------------------------------------------------------
  subroutine param_int_mpi(PARAMS, section, keyword, params_int, defaultvalue)
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    integer :: params_int, defaultvalue
    integer :: mpicode

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call read_param(PARAMS, section, keyword, params_int, defaultvalue)
    endif

    ! And then broadcast
    call MPI_BCAST( params_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_int_mpi
end module
