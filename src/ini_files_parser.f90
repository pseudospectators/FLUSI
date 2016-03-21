module ini_files_parser
  use vars

  integer, parameter :: maxcolumns=1024

  ! is set to true, we'll produce some output on the screen (for documentation of runs)
  ! the flag is set with the read_ini_file routine
  logical, save :: verbosity = .true.

  ! type definition of an inifile type. it contains an allocatable string array
  ! of maxcolumns length. we will allocate the array, then fill it with what we
  ! read from the file.
  type inifile
    ! string array that contains the text file:
    character(len=maxcolumns),allocatable,dimension(:) :: PARAMS
    ! number of lines in file:
    integer :: nlines
  end type

  ! the generic call "read_param" redirects to these routines, depending on the data
  ! type and the dimensionality. vectors can be read without setting a default.
  interface read_param
    module procedure param_dbl, param_int, param_vct, param_str, param_vct_nodefault
  end interface


!!!!!!!!
contains
!!!!!!!!


  !-------------------------------------------------------------------------------
  ! clean a previously read ini file, deallocate its string array, and reset
  ! verbosity to .true. (as a matter of precaution)
  !-------------------------------------------------------------------------------
  subroutine clean_ini_file(PARAMS)
    implicit none
    type(inifile), intent(inout) :: PARAMS

    if (allocated(PARAMS%PARAMS)) deallocate(PARAMS%PARAMS)
    verbosity = .true.
  end subroutine clean_ini_file

  !-------------------------------------------------------------------------------
  ! Read the file paramsfile, count the lines and put the
  ! text in PARAMS.
  !-------------------------------------------------------------------------------
  subroutine read_ini_file(PARAMS, file, verbose)
    use vars
    implicit none

    type(inifile), intent(inout) :: PARAMS
    character(len=*) :: file ! this is the file we read the PARAMS from
    character(len=maxcolumns) :: dummy, line
    logical, intent(in) :: verbose
    integer :: io_error, i

    ! check if the specified file exists
    call check_file_exists( file )

    ! we set the module-global variable verbosity. if set to false, all routines
    ! will perform their task quietly.
    verbosity = verbose

    !-----------------------------------------------------------------------------
    ! determine number of lines in file
    !-----------------------------------------------------------------------------
    io_error=0
    if (verbosity) then
      write (*,*) "*************************************************"
      write (*,'(A,i3)') " *** info: reading params from "//trim(file)//" rank=", mpirank
      write (*,*) "*************************************************"
    endif
    i = 0

    open(unit=14,file=trim(adjustl(file)),action='read',status='old')
    do while ((io_error==0).and.(i<=nlines))
      read (14,'(A)',iostat=io_error) dummy
      i = i+1
    enddo
    close (14)
    PARAMS%nlines = i

    allocate( PARAMS%PARAMS(1:PARAMS%nlines) )

    !-----------------------------------------------------------------------------
    ! read actual file
    !-----------------------------------------------------------------------------
    io_error=0
    i=1
    if (mpirank==0) then
      open(unit=14,file=trim(adjustl(file)),action='read',status='old')
      do while (io_error==0)
        read (14,'(A)',iostat=io_error) line
        PARAMS%PARAMS(i) = adjustl(line)
        i=i+1
      enddo
      close (14)
    endif
  end subroutine read_ini_file






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
  subroutine param_dbl (PARAMS, section, keyword, params_real, defaultvalue)
    use vars
    use mpi
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) :: value    ! returns the value
    real (kind=pr) :: params_real, defaultvalue
    integer :: mpicode

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then

      call GetValue(PARAMS, section, keyword, value)

      if (value .ne. '') then
        if ( index(value,'*dx') == 0 ) then
          !-- value to be read is in absolute form (i.e. we just read the value)
          read (value, *) params_real
          write (value,'(g10.3)') params_real
          if (verbosity) then
            write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
          endif
        else
          !-- the value is given in gridpoints (e.g. thickness=5*dx)
          read (value(1:index(value,'*dx')-1),*) params_real
          params_real = params_real*max(dy,dz)
          write (value,'(g10.3,"(=",g10.3,"*dx)")') params_real, params_real/max(dy,dz)
          if (verbosity) then
            write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
          endif
        endif
      else
        write (value,'(g10.3)') defaultvalue
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
        endif
        params_real = defaultvalue
      endif
    endif

    ! And then broadcast
    call MPI_BCAST( params_real, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_dbl




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
  subroutine param_str (PARAMS, section, keyword, params_string, defaultvalue)
    use vars
    use mpi
    implicit none

    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! what section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns)  value    ! returns the value
    character(len=strlen), intent (inout) :: params_string
    character(len=*), intent (in) :: defaultvalue
    integer :: mpicode

    !------------------
    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    !------------------
    if (mpirank==0) then
      call GetValue(PARAMS, section, keyword, value)
      if (value .ne. '') then
        params_string = value
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
      else
        value = defaultvalue
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
        endif
        params_string = defaultvalue
      endif
    endif

    ! And then broadcast
    call MPI_BCAST( params_string, strlen, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpicode)
  end subroutine param_str



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
  subroutine param_vct (PARAMS, section, keyword, params_vector, defaultvalue, n)
    use vars
    use mpi
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=pr) :: params_vector(1:n)
    real(kind=pr) :: defaultvalue(1:n)
    integer, intent(in) :: n

    character(len=maxcolumns) :: value
    character(len=14)::formatstring
    integer :: mpicode

    if(n==0) return
    write(formatstring,'("(",i2.2,"(g10.3,1x))")') n

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call GetValue(PARAMS, section, keyword, value)
      if (value .ne. '') then
        ! read the three values from the vector string
        read (value, *) params_vector
        write (value,formatstring) params_vector
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
      else
        write (value,formatstring) defaultvalue
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
        endif
        params_vector = defaultvalue
      endif
    endif

    ! And then broadcast
    call MPI_BCAST( params_vector, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_vct


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
  subroutine param_vct_nodefault (PARAMS, section, keyword, params_vector, n)
    use vars
    use mpi
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    real(kind=pr) :: params_vector(1:n)
    real(kind=pr) :: defaultvalue(1:n)
    integer, intent(in) :: n

    character(len=maxcolumns) :: value
    character(len=14)::formatstring

    integer :: mpicode
    defaultvalue=0.0;

    if(n==0) return
    write(formatstring,'("(",i2.2,"(g10.3,1x))")') n

    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    if (mpirank==0) then
      call GetValue(PARAMS, section, keyword, value)
      if (value .ne. '') then
        ! read the three values from the vector string
        read (value, *) params_vector
        write (value,formatstring) params_vector
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
      else
        write (value,formatstring) defaultvalue
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
        endif
        params_vector = defaultvalue
      endif
    endif

    ! And then broadcast
    call MPI_BCAST( params_vector, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_vct_nodefault



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
  subroutine param_int(PARAMS, section, keyword, params_int, defaultvalue)
    use mpi
    use vars
    implicit none
    ! Contains the ascii-params file
    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword ! what keyword do you look for? for example nx=128
    character(len=maxcolumns) ::  value    ! returns the value
    integer :: params_int, defaultvalue
    integer :: mpicode

    !------------------
    ! Root rank fetches value from PARAMS.ini file (which is in PARAMS)
    !------------------
    if (mpirank==0) then
      call GetValue(PARAMS, section, keyword, value)
      if (value .ne. '') then
        read (value, *) params_int
        write (value,'(i7)') params_int
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))
        endif
      else
        write (value,'(i7)') defaultvalue
        if (verbosity) then
          write (*,*) "read "//trim(section)//"::"//trim(keyword)//" = "//adjustl(trim(value))//&
          " (THIS IS THE DEFAULT VALUE!)"
        endif
        params_int = defaultvalue
      endif
    endif

    ! And then broadcast
    call MPI_BCAST( params_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )
  end subroutine param_int





  !-------------------------------------------------------------------------------
  ! Extracts a value from the PARAMS.ini file, which is in "section" and
  ! which is named "keyword"
  ! Input:
  !       PARAMS: the complete *.ini file
  !       section: the section we're looking for
  !       keyword: the keyword we're looking for
  ! Output:
  !       value: is a string contaiing everything between '=' and ';'
  !              to be processed further, depending on the expected type
  !              of variable (e.g. you read an integer from this string)
  !-------------------------------------------------------------------------------
  subroutine GetValue (PARAMS, section, keyword, value)
    use vars
    implicit none

    type(inifile), intent(inout) :: PARAMS
    character(len=*), intent(in) :: section ! What section do you look for? for example [Resolution]
    character(len=*), intent(in) :: keyword   ! what keyword do you look for? for example nx=128
    character(len=*), intent(inout) :: value   ! returns the value
    integer :: i
    integer :: index1,index2
    logical :: foundsection

    foundsection = .false.
    value = ''

    !-- loop over the lines of PARAMS.ini file
    do i=1, PARAMS%nlines
      !-- ignore commented lines completely
      if ((PARAMS%PARAMS(i)(1:1).ne.'#').and.(PARAMS%PARAMS(i)(1:1).ne.';').and.&
      (PARAMS%PARAMS(i)(1:1).ne.'!').and.(PARAMS%PARAMS(i)(1:1).ne.'%')) then

      !-- does this line contain the "[section]" statement?
      if (index(PARAMS%PARAMS(i),'['//section//']')==1) then
        ! yes, it does
        foundsection = .true.
      elseif (PARAMS%PARAMS(i)(1:1) == '[') then
        ! we're already at the next section mark, so we leave the section we
        ! were looking for again
        foundsection = .false.
      endif

      !-- we're inside the section we want
      if (foundsection) then
        if (index(PARAMS%PARAMS(i),keyword//'=')==1) then
          if (index(PARAMS%PARAMS(i),';')/=0) then
            ! found "keyword=" in this line, as well as the delimiter ";"
            index1 = index(PARAMS%PARAMS(i),'=')+1
            index2 = index(PARAMS%PARAMS(i),';')-1
            value = PARAMS%PARAMS(i)(index1:index2)
            exit ! leave do loop
          else
            ! found "keyword=" in this line, but delimiter ";" is missing.
            ! proceed, but this may fail (e.g. value=999commentary or
            ! value=999\newline fails)
            if (verbosity) then
              write (*,'(A)') "Though found the keyword, I'm unable to find &
                            &value for variable --> "//trim(keyword)//" <-- &
                            &missing delimiter (;)"
            endif
            index1 = index(PARAMS%PARAMS(i),'=')+1
            index2 = index(PARAMS%PARAMS(i),' ')-1
            value = PARAMS%PARAMS(i)(index1:index2)
            if (verbosity) then
              write (*,'(A)') "As you forgot to set ; at the end of your value,&
                          & I try to extract a value nonetheless. Proceed, with&
                          & fingers crossed (there is no guarantee this will work)"
              write (*,'(A)') "Extracted --->"//trim(value)//"<-----"
            endif
            exit
          endif
        endif
      endif

    endif
  enddo
end subroutine getvalue




end module
