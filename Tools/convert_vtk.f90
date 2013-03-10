!=========================================================
! this file converts several *.mpiio files to one ascii *.vtk file
!
! run it from command line: ./convert_vtk outfile.vtk 64 64 64 ux.mpiio uy.mpiio uz_mpiio [other files]
!
! note program checks if files are present, if not, they are skipped (yet still program continues!)
! 
! note you can't store velocity or vorticity if one of the 3 components is missing.
!
! files must begin with their variable identifier: ux_500.mpiio is valid, mpiio_ux_500 is not!
! order of arguments is arbitrary (besides outfile, nx,ny,nz)
!
!=========================================================
! THOMAS ENGELS, Aix-Marseille Université and Technische Universität Berlin
! MARCH 2013
!=========================================================

program convert_mpiio
  use mpi
  implicit none

  integer, parameter :: pr_in = 8, pr_out = 4 ! precision for input and output (INPUT PRECISION MUST BE KNOWN!!!) 
  integer :: nx, ny, nz, mpicode, nunit
  integer :: ix, iy, iz, iSaveAscii, iArg
  real (kind=pr_out), dimension (:,:,:), allocatable :: field_out
  real (kind=pr_out), dimension (:,:,:), allocatable :: ux,uy,uz
  character (len=128) :: fname, nx_str, ny_str, nz_str, fname_out
  character (len=128) :: fname_ux,fname_uy,fname_uz,fname_vorx,fname_vory,fname_vorz, fname_mask, fname_p
  logical :: bool_ux=.false., bool_uy=.false., bool_uz=.false., bool_vorx=.false.
  logical :: bool_vorz=.false., bool_vory=.false., bool_mask=.false., bool_p=.false., file_exists, check

  !--------------------------------------------------
  ! read in command line arguments=.false.
  !--------------------------------------------------  
  call get_command_argument(1, fname_out)
  if  (len_trim(nx_str)==0) then
    write (*,*) "I'm confused... whats the output filename?"
    call Help
    stop    
  endif
  
  
  call get_command_argument(2, nx_str)
  if  (len_trim(nx_str).ne.0) then
    read (nx_str, *) nx
  else
    write (*,*) "I'm confused... whats nx?"
    call Help
    stop    
  endif
  
  call get_command_argument(3, ny_str)
  if  (len_trim(ny_str).ne.0) then
    read (ny_str, *) ny
  else
    write (*,*) "I'm confused... whats ny?"
    call Help
    stop    
  endif
  
  call get_command_argument(4, nz_str)
  if  (len_trim(nz_str).ne.0) then
    read (nz_str, *) nz
  else
    write (*,*) "I'm confused... whats nz?"
    call Help
    stop    
  endif
  
  !----------------------------------------------------------
  ! retrieve all the filenames (their order is arbitrary, also their presence)
  !----------------------------------------------------------
  do iArg = 5,12
    call get_command_argument(iArg, fname)
    if (len_trim(fname).ne.0) then
      select case (fname(1:3))
	  case ('ux_')
	      write(*,*) 'identified ', trim(fname), ' as ux'
	      bool_ux = check(fname)
	      fname_ux = trim(fname)
	  case ('uy_')
	      write(*,*) 'identified ', trim(fname), ' as uy'
	      bool_uy = check(fname)
	      fname_uy = trim(fname)
	  case ('uz_')
	      write(*,*) 'identified ', trim(fname), ' as uz'
	      bool_uz = check(fname)
	      fname_uz = trim(fname)
      end select
      
      select case (fname(1:4))
	  case ('vorx')
	      write(*,*) 'identified ', trim(fname), ' as vorx'
	      bool_vorx = check(fname)
	      fname_vorx = trim(fname)
	  case ('vory')
	      write(*,*) 'identified ', trim(fname), ' as vory'
	      bool_vory = check(fname)
	      fname_vory = trim(fname)
	  case ('vorz')
	      write(*,*) 'identified ', trim(fname), ' as vorz'
	      bool_vorz = check(fname)
	      fname_vorz = trim(fname)
	  case ('mask')
	      write(*,*) 'identified ', trim(fname), ' as mask'
	      bool_mask = check(fname)
	      fname_mask = trim(fname)
      end select
      
      if (fname(1:2)=="p_") then
	  write(*,*) 'identified ', trim(fname), ' as p'
	  bool_p = check(fname)
	  fname_p = trim(fname)
      endif      
    endif  
  enddo
  
  call MPI_INIT (mpicode)
  allocate ( field_out(0:nx-1, 0:ny-1, 0:nz-1) )
  
  !---------------------------------------------------------------
  ! read in the files and write the VTK file
  !---------------------------------------------------------------
  
  nunit = 11
  open (nunit, file=fname_out, form='formatted', status='replace')

  ! write some header stuff
  write(unit=nunit,fmt="('# vtk DataFile Version 2.0')")
  write(unit=nunit,fmt="('flusi-data')")
  write(unit=nunit,fmt="('ASCII')")
  write(unit=nunit,fmt="('DATASET STRUCTURED_GRID')")
  write(unit=nunit,fmt="('DIMENSIONS ',3I8)")nx,ny,nz
  write(unit=nunit,fmt="('POINTS',I8,A6)")nx*ny*nz,"float"
  ! save the grid (I know it sucks)
  do iz=0,nz-1
      do iy=0,ny-1
	do ix=0,nx-1
	    write(unit=nunit,fmt='(3G16.6)') real(ix),real(iy),real(iz)
	end do
      end do
  end do
  write(unit=nunit,fmt="('POINT_DATA',I8)")nx*ny*nz
  
  ! pressure
  if (bool_p) then
    call ReadMPIIO (nx, ny, nz, mpicode, fname_p, field_out)
    write(unit=nunit,fmt="('SCALARS p FLOAT')")
    write(unit=nunit,fmt="('LOOKUP_TABLE default')")
    do iz=0,nz-1
	do iy=0,ny-1
	  do ix=0,nx-1
	      write(unit=nunit,fmt='(G13.6)') field_out(ix,iy,iz)
	  end do
	end do
    end do  
  endif
  
  ! mask
  if (bool_mask) then
    call ReadMPIIO (nx, ny, nz, mpicode, fname_mask, field_out)
    write(unit=nunit,fmt="('SCALARS mask FLOAT')")
    write(unit=nunit,fmt="('LOOKUP_TABLE default')")
    do iz=0,nz-1
	do iy=0,ny-1
	  do ix=0,nx-1
	      write(unit=nunit,fmt='(G13.6)') field_out(ix,iy,iz)
	  end do
	end do
    end do  
  endif
  
  ! velocity
  if ((bool_ux).and.(bool_uy).and.(bool_uz)) then
      allocate ( ux(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uy(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uz(0:nx-1, 0:ny-1, 0:nz-1) )
      call ReadMPIIO (nx, ny, nz, mpicode, fname_ux, ux)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_uy, uy)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_uz, uz)
      
      write(unit=nunit,fmt="('VECTORS u FLOAT')")
      do iz=0,nz-1
	  do iy=0,ny-1
	    do ix=0,nx-1
		write(unit=nunit,fmt='(3G16.6)')&
		ux(ix,iy,iz),uy(ix,iy,iz),uz(ix,iy,iz)
	    end do
	  end do
      end do      
      
      deallocate (ux,uy,uz)
  endif

  ! vorticity
  if ((bool_vorx).and.(bool_vory).and.(bool_vorz)) then
      allocate ( ux(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uy(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uz(0:nx-1, 0:ny-1, 0:nz-1) )
      call ReadMPIIO (nx, ny, nz, mpicode, fname_vorx, ux)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_vory, uy)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_vorz, uz)
      
      write(unit=nunit,fmt="('VECTORS vor FLOAT')")
      do iz=0,nz-1
	  do iy=0,ny-1
	    do ix=0,nx-1
		write(unit=nunit,fmt='(3G16.6)')&
		ux(ix,iy,iz),uy(ix,iy,iz),uz(ix,iy,iz)
	    end do
	  end do
      end do      
      
      deallocate (ux,uy,uz)
  endif

  close(unit=nunit)
  
  deallocate ( field_out )  
  call MPI_FINALIZE (mpicode)

end program convert_mpiio


subroutine help
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "		converter "
  write (*,*) " *.mpiio -> *.vtk (for PARAVIEW)"
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "  usage: ./convert_vtk.out file_out nx ny nz file1 [file2]...[file11]"
  write (*,*) ""
  write (*,*) " for example:"
  write (*,*) " ./convert_vtk.out test.vtk 64 128 64 ux_001.00.mpiio uy_001.00.mpiio p_001.00.mpiio"
  write (*,*) "--------------------------------------------------------------"

end subroutine



subroutine ReadMPIIO (nx, ny, nz, mpicode, fname, field_out)
  use mpi
  implicit none
  character(len=128), intent (in) :: fname 
  integer, intent(in) :: nx,ny,nz
  integer, intent(inout) :: mpicode
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus
  integer :: filedesc
  real(kind=4), intent(out) :: field_out (0:nx-1, 0:ny-1, 0:nz-1)
  real(kind=8), dimension(:,:,:), allocatable :: field_in
  integer, parameter :: mpireal = MPI_DOUBLE_PRECISION  ! double precision array for input
   
  
  !--------------------------------------------------
  ! read MPI data
  !-------------------------------------------------- 
  allocate ( field_in(0:nx-1, 0:ny-1, 0:nz-1) )
  call MPI_FILE_OPEN (MPI_COMM_WORLD,trim(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,filedesc,mpicode)
  call MPI_FILE_READ_ORDERED (filedesc,field_in,nx*ny*nz,mpireal,mpistatus,mpicode)
  call MPI_FILE_CLOSE (filedesc,mpicode)

  write(*,'("Reading ",A,"... Min:Max=",es12.4,":",es12.4,1x,"nx:ny:nz=",i3,":",i3,":",i3)') &
  trim(fname),minval (field_in), maxval (field_in),nx,ny,nz
  
  
  field_out = real(field_in, kind=4)  ! convert field to output precision  
  
  deallocate (field_in)
end subroutine 



!-------------------------------------------
! check if file exists (returns true if so)
!-------------------------------------------
logical function check(fname)
  implicit none
  character(len=128), intent (in) :: fname
  logical :: file_exists
  
  inquire( file=trim(fname),exist=file_exists )   
  if (file_exists .eqv. .false.) then
    write (*,'("ERROR!! File ",A," NOT found -------> skipping")') trim(fname)
  endif	  
  
end function