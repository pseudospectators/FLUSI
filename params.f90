subroutine get_params (paramsfile)
  use mpi_header ! Module incapsulates mpif.
  use share_vars
  implicit none
  integer                              :: mpicode, io_error,i
  integer, dimension (11)              :: comm_int
  real (kind=pr), dimension (11)       :: comm_real
  character (len=80)                   :: paramsfile ! this is the file we read the PARAMS from
  character (len=80)                   :: tmp
  character PARAMS(nlines)*256	! this array will contain the ascii-params file
 

  if (mpirank==0) then
    !-----------------------------------------------------------
    ! read in the params file (root rank only)
    !-----------------------------------------------------------
    write (*,'(A)') "*** info: reading params from "//trim(paramsfile)
    i = 1
    open ( unit=14, file=paramsfile, action='read', status='old' )    
    do while ((io_error==0).and.(i<=nlines))
      read (14,'(A)',iostat=io_error) PARAMS(i)  
      i = i+1
    enddo      
    close (14)
    i = i-1 ! counted one too far
    
    !-----------------------------------------------------------
    ! set default values (they are overwritten if specified in PARAMS.ini)
    !-----------------------------------------------------------
    nx = 4
    ny = 4
    nz = 4
    nt = 99999999
    Tmax = 9999999.9
    CFL = 0.1
    nu = 1.d-1
    inicond =1
    iDealias=0
    iPenalization=0
    iMoving=0
    iMask=0
    eps=0.0
    xl=1.0
    yl=1.0
    zl=1.0
    x0=0.0
    y0=0.0
    z0=0.0
    iMeanFlow=1
    tsave=999999999.99
    tdrag=999999999.99
    iDoBackup=1
    iSaveVelocity=0
    iSavePress=0
    iSaveVorticity=0
    iSaveMask=0
    iSaveSolidVelocity=0
    
    !-----------------------------------------------------------
    ! fetch all the parameters
    !-----------------------------------------------------------
    call GetValue(PARAMS,i,"Resolution","nx",tmp)
    if (tmp .ne. '') read (tmp, *) nx
    call GetValue(PARAMS,i,"Resolution","ny",tmp)
    if (tmp .ne. '') read (tmp, *) ny
    call GetValue(PARAMS,i,"Resolution","nz",tmp)
    if (tmp .ne. '') read (tmp, *) nz
    call GetValue(PARAMS,i,"Time","Tmax",tmp)
    if (tmp .ne. '') read (tmp, *) Tmax
    call GetValue(PARAMS,i,"Time","nt",tmp)
    if (tmp .ne. '') read (tmp, *) nt   
    call GetValue(PARAMS,i,"Time","iTimeMethodFluid",tmp)
    if (tmp .ne. '') read (tmp, *) iTimeMethodFluid   
    call GetValue(PARAMS,i,"Time","CFL",tmp)
    if (tmp .ne. '') read (tmp, *) CFL       
    call GetValue(PARAMS,i,"ReynoldsNumber","nu",tmp)
    if (tmp .ne. '') read (tmp, *) nu
    call GetValue(PARAMS,i,"InitialCondition","inicond",tmp)
    if (tmp .ne. '') read (tmp, *) nu    
    call GetValue(PARAMS,i,"Dealiasing","iDealias",tmp)
    if (tmp .ne. '') read (tmp, *) iDealias      
    call GetValue(PARAMS,i,"Penalization","iPenalization",tmp)
    if (tmp .ne. '') read (tmp, *) iPenalization        
    call GetValue(PARAMS,i,"Penalization","iMoving",tmp)
    if (tmp .ne. '') read (tmp, *) iMoving    
    call GetValue(PARAMS,i,"Penalization","iMask",tmp)
    if (tmp .ne. '') read (tmp, *) iMask    
    call GetValue(PARAMS,i,"Penalization","eps",tmp)
    if (tmp .ne. '') read (tmp, *) eps
    call GetValue(PARAMS,i,"Geometry","xl",tmp)
    if (tmp .ne. '') read (tmp, *) xl
    call GetValue(PARAMS,i,"Geometry","yl",tmp)
    if (tmp .ne. '') read (tmp, *) yl
    call GetValue(PARAMS,i,"Geometry","zl",tmp)
    if (tmp .ne. '') read (tmp, *) zl
    call GetValue(PARAMS,i,"Geometry","x0",tmp)
    if (tmp .ne. '') read (tmp, *) x0
    call GetValue(PARAMS,i,"Geometry","y0",tmp)
    if (tmp .ne. '') read (tmp, *) y0
    call GetValue(PARAMS,i,"Geometry","z0",tmp)
    if (tmp .ne. '') read (tmp, *) z0   
    call GetValue(PARAMS,i,"MeanFlow","iMeanFlow",tmp)
    if (tmp .ne. '') read (tmp, *) iMeanFlow 
    call GetValue(PARAMS,i,"Saving","tsave",tmp)
    if (tmp .ne. '') read (tmp, *) tsave     
    call GetValue(PARAMS,i,"Saving","tdrag",tmp)
    if (tmp .ne. '') read (tmp, *) tdrag     
    call GetValue(PARAMS,i,"Saving","iDoBackup",tmp)
    if (tmp .ne. '') read (tmp, *) iDoBackup     
    call GetValue(PARAMS,i,"Saving","iSaveVelocity",tmp)
    if (tmp .ne. '') read (tmp, *) iSaveVelocity     
    call GetValue(PARAMS,i,"Saving","iSavePress",tmp)
    if (tmp .ne. '') read (tmp, *) iSavePress     
    call GetValue(PARAMS,i,"Saving","iSaveVorticity",tmp)
    if (tmp .ne. '') read (tmp, *) iSaveVorticity     
    call GetValue(PARAMS,i,"Saving","iSaveMask",tmp)
    if (tmp .ne. '') read (tmp, *) iSaveMask     
    call GetValue(PARAMS,i,"Saving","iSaveSolidVelocity",tmp)
    if (tmp .ne. '') read (tmp, *) iSaveSolidVelocity         
  endif
  
  !-------------------------------------------------------
  ! broadcast all parameters
  !-------------------------------------------------------
  call MPI_BCAST( nx, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( ny, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( nz, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( nt, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( inicond, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iTimeMethodFluid, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iDealias, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iMoving, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iMask, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iMeanFlow, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iDoBackup, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iSaveVelocity, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iSavePress, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iSaveVorticity, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iSaveMask, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( iSaveSolidVelocity, 1, mpiinteger, 0, MPI_COMM_WORLD, mpicode )
  
  call MPI_BCAST( Tmax, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( CFL, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( nu, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( eps, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( xl, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( yl, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( zl, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( x0, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( y0, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( z0, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( tsave, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( tdrag, 1, mpireal, 0, MPI_COMM_WORLD, mpicode )
  
  !-------------------------------------------------------
  ! set other parameters
  !-------------------------------------------------------  
  pi     = 4.d0 * datan (1.d0)
  scalex = 2.d0*pi / xl
  scaley = 2.d0*pi / yl
  scalez = 2.d0*pi / zl  
  dx = xl / dble (nx)
  dy = yl / dble (ny)
  dz = zl / dble (nz) 
  tstart = 0.d0
  
  
  !-- Stop if the sizes are odd or smaller than 4
  if ( nx<4 .or. ny<4 .or. nz<4 .or. modulo(nx,2)==1 .or. modulo(ny,2)==1 .or. modulo(nz,2)==1 ) then
    if ( mpirank == 0 ) then
      print *, 'nx, ny, nz must be even and not smaller than 4'
    endif
    stop
  endif
  
  !-------------------------------------------------------
  ! initialize FFT
  !-------------------------------------------------------
  call fft_initialize 
  
  !--------------------------------------------------------
  ! HEADER
  !--------------------------------------------------------
  if ( mpirank == 0 ) then
     write (*,*) '--------------------------------------------'
     write (*, '(" nx = ",i4, 1x, "ny = ", i4, 1x, "nz = ", i4)') nx, ny, nz
     write (*, '(" xl = ", f6.2, 1x, "yl = ", f6.2, 1x, "zl = ", f6.2)') xl, yl, zl
     write (*, '(" nt = ", i6, 1x)') nt
     write (*, '(" tmax = ", f6.2)') tmax
     write (*, '(" viscosity = ", es11.4)') nu
     write (*,*) '--------------------------------------------'
     write (*, '(" calculate drag every", es11.4, " time steps")') tdrag
     write (*, '(" save fields every   ", es11.4)') tsave
     if (iSavePress == 1) write (*,*) 'save pressure'
     if (iSaveVelocity == 1) write (*,*) 'save velocity'
     if (iSaveVorticity == 1) write (*,*) 'save vorticity'
     if (iSaveMask == 1) write (*,*) 'save mask'
     if (iSaveSolidVelocity == 1) write (*,*) 'save mask velocity'
     write (*,*) '--------------------------------------------'
     if (iPenalization > 0) then 
        write (*,*) 'with obstacle'
        write (*, '(" x0 = ", f6.2, 1x, "y0 = ", f6.2, 1x, "z0 = ", f6.2, 1x, &
            & "size = ", f6.2)') x0, y0, z0, size
        write (*, '(" epsilon_eff   = ", es11.4)') eps * xl / real (nx)
        write (*, '(" epsilon       = ", es11.4)') eps
        write (*, '(" Ux = ", f6.2, " Uy = ", f6.2, " Uz = ", f6.2)') Ux, Uy, Uz
        write (*, '(" Ax = ", f6.2, " Ay = ", f6.2, " Az = ", f6.2)') Ax, Ay, Az
        write (*, '(" Re = ", es11.4)') sqrt(Ux**2+Uy**2+Uz**2) * size / nu
     else
        write (*,*) 'without obstacle'
     end if
     write (*,*) ' '
     write (*,*) '--------------------------------------------'
     write (*,*) ' '
  endif
  
  stop "hhaha"
  
end subroutine get_params



subroutine GetValue (PARAMS, actual_lines, section, keyword, value)
  use share_vars
  implicit none
  character section*(*)			! what section do you look for? for example [Resolution]
  character keyword*(*)			! what keyword do you look for? for example nx=128
  character value*(*)			! returns the value
  character PARAMS(nlines)*256		! this is the complete PARAMS.ini file
  integer actual_lines			! how many lines did you actually read?  
  integer :: maxline = 256			! how many characters per line?
  integer i, j,k
  character line*256
  logical foundsection
  
  foundsection = .false.
  value = ''

  do i=1, actual_lines					! loop over the lines of PARAMS.ini file
  if (PARAMS(i)(1:1) .ne. '#') then			! ignore commented lines compleetly
  
  
      if (PARAMS(i)(1:1) == '[') then			! the first char would have to be '['
	  do j = 2, maxline				! then we look fot the corrresponding ']'
	    if (PARAMS(i)(j:j) == ']') then		! we found it
	      if (section == PARAMS(i)(2:j-1)) then	! is this the section we"re looking for?
		foundsection = .true.			! yes, it is
		exit
	      endif
	    endif
	  enddo	
      else      
	  if (foundsection .eqv. .true.) then		! yes we found the section, now we're looking for the keyword
	    do j=1, maxline				! scan the line
	      if (PARAMS(i)(j:j) == '=') then		! found the '='
	      if (keyword == PARAMS(i)(1:j-1)) then	! is this the keyword you're looking for?
		do k = j+1, maxline			! everything behind the '=' and before ';' is the value
		if (PARAMS(i)(k:k) == ';') then		! found the delimiter
		    value = PARAMS(i)(j+1:k-1)		! value is between '=', and ';'   
		    exit
		endif
		enddo
		if (value == '') then
		  write (*,'(A)') "??? Though found the keyword, I'm unable to find value for variable --> "& 
		               //trim(keyword)//" <-- maybe missing delimiter (;)?"		  
		endif
	      endif
	      endif
	    enddo
	  endif
      endif
      
      
  endif
  enddo
  
  if (value =='') then
    write(*,'(A)') "??? WARNING: No value found for "//trim(keyword)//" in section "//trim(section)//" ----> using default value!"
  endif


end subroutine getvalue
