!=========================================================
! converts *.mpiio files to paraview readable *.vtk format
!
! can use either scalar or vector mode
!
! SCALAR:
! 	./convert_vtk.out nx ny nz -scalar mask_00010.mpiio
!	this creates the mask_00010.vtk file
!
! VECTOR:
!	./convert_vtk.out nx ny nz -vector ux_00010.mpiio uy_00010.mpiio uz_00010.mpiio
!	this creates ONE file u_00010.vtk
!
!=========================================================
! THOMAS ENGELS, Aix-Marseille Université and Technische Universität Berlin
! APRIL 2013
!=========================================================

program convert_mpiio
  use mpi
  implicit none

  integer, parameter :: pr_in = 8, pr_out = 4 ! precision for input and output (INPUT PRECISION MUST BE KNOWN!!!) 
  integer :: nx, ny, nz, mpicode, nunit
  integer :: ix, iy, iz, iSaveAscii, iArg, ii, j1,j2
  real (kind=pr_out), dimension (:,:,:), allocatable :: field_out
  real (kind=pr_out), dimension (:,:,:), allocatable :: ux,uy,uz
  real (kind=pr_out), dimension (:,:), allocatable :: tmp
  character (len=128) :: fname, nx_str, ny_str, nz_str, fname_out, modus
  character (len=128) :: fname_x, fname_y, fname_z
  character (len=128) :: cbuffer, frmt 

  
  !--------------------------------------------------
  ! read in command line arguments
  !--------------------------------------------------  
  call get_command_argument(1, nx_str)
  if  (len_trim(nx_str).ne.0) then
    read (nx_str, *) nx
  else
    call help
    stop "nx not specified!"
  endif
  
  
  call get_command_argument(2, ny_str)
  if  (len_trim(ny_str).ne.0) then
    read (ny_str, *) ny
  else
    call help
    stop "ny not specified!"
  endif
  
  
  call get_command_argument(3, nz_str)
  if  (len_trim(nz_str).ne.0) then
    read (nz_str, *) nz
  else
    call help
    stop "nz not specified!" 
  endif
  
  
  call get_command_argument(4, modus)
  if  (len_trim(modus)==0) then
    call help
    stop "you forgot to set the modus, either -scalar or -vector?"
  endif
  
  
  call MPI_INIT (mpicode)
  
  !-----------------------------------------
  ! scalar modus
  !-----------------------------------------  
  if ( modus=='-scalar' ) then
  
      !-------------------------------------------
      ! get filename (only one in this case)
      !-------------------------------------------
      call get_command_argument( 5, fname )        
      if (len_trim(fname)==0) stop ("no file given")
      call check(fname)
      
      ! fetch variable name "mask_00010.mpiio"
      j1 = index ( fname,"_" ) ! name is fname(1:j1-1)     "mask"
      j2 = index ( fname,"." ) ! filename is fname(1:j2-1) "mask_00010"
      
      !---------------------------------------------------------------
      ! read in the files and write the VTK file
      !---------------------------------------------------------------
      
      nunit = 11
      open (nunit, file=fname(1:j2-1)//".vtk", status='replace', access='STREAM',convert='big_endian')

      !----------------------------
      ! write VTK file
      !----------------------------  
      write(nunit) "# vtk DataFile Version 3.0"//char(10)		! standard header
      write(nunit) "Fluid Field"//char(10)				! description
      write(nunit) "BINARY"//char(10)					! it's a binary file
      
      write(nunit) "DATASET STRUCTURED_POINTS"//char(10)		! we use uniform grids
      write(cbuffer,fmt="('DIMENSIONS ',3(i5,1x))") nx,ny,nz		! this is the grid dimensions
      write(nunit) cbuffer//char(10)
      write(cbuffer,fmt="('ORIGIN ',3(i1,1x))") 0, 0, 0 		! our origin is 0,0,0 
      write(nunit) cbuffer//char(10) 
      write(cbuffer,fmt="('SPACING ',3(i1,1x))") 1, 1, 1		! our spacing is 1,1,1 (note it's possible to use the real spacing with the domain size)
      write(nunit) cbuffer//char(10) 
      
      write(cbuffer,fmt="('POINT_DATA',I10)") nx*ny*nz			! tell paraview to expect point data, not cell data
      write(nunit) cbuffer//char(10)

      write(cbuffer,fmt="('SCALARS ',A10,' float 1')") fname(1:j1-1)	! save ONE scalar field, named as the variable name (the last "1" is the dim of scalar field)
      write(nunit) cbuffer//char(10)
      write(cbuffer,fmt="('LOOKUP_TABLE default')")			! standard lookup table (don't touch)
      write(nunit) cbuffer//char(10)
      
      
      allocate ( field_out(0:nx-1, 0:ny-1, 0:nz-1) )			! alloc output field
      
      call ReadMPIIO (nx, ny, nz, mpicode, fname, field_out)		! read in file (field_out is single precision)
      
      write(nunit) field_out						! dump to vtk file
      
      close(unit=nunit)							! close vtk file
    
      deallocate ( field_out )  
      
      
      
  !-----------------------------------------
  ! vector modus
  !-----------------------------------------  
  elseif ( modus=='-vector') then
  
      !-------------------------------------------
      ! get filenames (three files in this case)
      !-------------------------------------------
      call get_command_argument( 5, fname_x )        
      if (len_trim(fname_x)==0) stop ("no x-file given")
      call check(fname_x)
      
      call get_command_argument( 6, fname_y )        
      if (len_trim(fname_y)==0) stop ("no y-file given")
      call check(fname_x)
      
      call get_command_argument( 7, fname_z )        
      if (len_trim(fname_z)==0) stop ("no z-file given")
      call check(fname_x)
      
      ! fetch variable name "vorx_00010.mpiio" "vory_00010.mpiio" "vorz_00010.mpiio"
      j1 = index ( fname_x,"x" ) ! name is fname(1:j1-1)     "vor"
      ii = index ( fname_x,"_" )
      j2 = index ( fname_x,"." ) ! filename is fname(1:j2-1) "vorx_00010"
      
      fname_out = fname_x(1:j1-1)//fname_x(ii:j2-1)//".vtk" ! out name is "vor_00010.vtk"
      
      !---------------------------------------------------------------
      ! read in the files and write the VTK file
      !---------------------------------------------------------------
      
      nunit = 11
      open (nunit, file=fname_out, status='replace', access='STREAM',convert='big_endian')

      !----------------------------
      ! write VTK file
      !----------------------------  
      write(nunit) "# vtk DataFile Version 3.0"//char(10)		! standard header
      write(nunit) "Fluid Field"//char(10)				! description
      write(nunit) "BINARY"//char(10)					! it's a binary file
      
      write(nunit) "DATASET STRUCTURED_POINTS"//char(10)		! we use uniform grids
      write(cbuffer,fmt="('DIMENSIONS ',3(i5,1x))") nx,ny,nz		! this is the grid dimensions
      write(nunit) cbuffer//char(10)
      write(cbuffer,fmt="('ORIGIN ',3(i1,1x))") 0, 0, 0 		! our origin is 0,0,0 
      write(nunit) cbuffer//char(10) 
      write(cbuffer,fmt="('SPACING ',3(i1,1x))") 1, 1, 1		! our spacing is 1,1,1 (note it's possible to use the real spacing with the domain size)
      write(nunit) cbuffer//char(10) 
      
      write(cbuffer,fmt="('POINT_DATA',I10)") nx*ny*nz			! tell paraview to expect point data, not cell data
      write(nunit) cbuffer//char(10)

      write(cbuffer,fmt="('VECTORS ',A10,' float')") fname_x(1:j1-1)	! save ONE scalar field, named as the variable name
      write(nunit) cbuffer//char(10)
      
      
      allocate ( ux(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uy(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( uz(0:nx-1, 0:ny-1, 0:nz-1) )
      allocate ( tmp(3,nx*ny*nz) )
      
      call ReadMPIIO (nx, ny, nz, mpicode, fname_x, ux)			! read in X-file (single precision)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_y, uy)			! read in Y-file (single precision)
      call ReadMPIIO (nx, ny, nz, mpicode, fname_z, uz)			! read in Z-file (single precision)
      
      ii = 0 ! counter
      do iz=0,nz-1
	  do iy=0,ny-1
	    do ix=0,nx-1
		ii = ii + 1 
		tmp(1,ii) = ux(ix,iy,iz)
		tmp(2,ii) = uy(ix,iy,iz)
		tmp(3,ii) = uz(ix,iy,iz)    
	    end do
	  end do
      end do
      
      write (nunit) tmp 						! dump to disk
      
      close(unit=nunit)							! close vtk file
    
      deallocate (ux,uy,uz,tmp)
  
  
  else
  
    write (*,*) "wrong modus, choose -scalar or -vector"
    stop  
  
  endif
  
  
  
  call MPI_FINALIZE (mpicode)

end program convert_mpiio


subroutine help
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "		converter "
  write (*,*) " *.mpiio -> *.vtk (for PARAVIEW)"
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "  usage: ./convert_vtk.out  nx ny nz -modus file1 [file2] [file3]"
  write (*,*) ""
  write (*,*) " for example (SCALAR):"
  write (*,*) " ./convert_vtk.out 64 128 64 -scalar mask_00010.mpiio"
  write (*,*) ""
  write (*,*) " for example (VECTOR):"
  write (*,*) " ./convert_vtk.out 64 128 64 -vector vorx_00010.mpiio vory_00010.mpiio vorz_00010.mpiio"
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
subroutine check(fname)
  implicit none
  character(len=128), intent (in) :: fname
  logical :: file_exists
  
  inquire( file=trim(fname),exist=file_exists )   
  if (file_exists .eqv. .false.) then
    write (*,'("ERROR!! File ",A," NOT found.")') trim(fname)
    stop
  endif	  
  
  
end subroutine