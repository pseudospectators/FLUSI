!-------------------------------------------------------------------------------
! Wrapper for different postprocessing tools
!-------------------------------------------------------------------------------
subroutine postprocessing()
  use fsi_vars
  use mpi_header
  implicit none
  character (len=80)     :: postprocessing_mode, filename, key1,key2
  
  if (mpirank==0) write (*,*) "*** FLUSI is running in postprocessing mode ***"
  
  
  ! the second argument tells us what to do with the file
  call get_command_argument(2,postprocessing_mode)
  ! it then depends on the second argument what follows
  
  !-----------------
  ! check what to do
  !-----------------     
  select case (postprocessing_mode)
  case ("--keyvalues")
    call get_command_argument(3,filename)
    call Keyvalues (filename)      
  case ("--compare-keys")
    call get_command_argument(3,key1)
    call get_command_argument(4,key2)
    call Compare_key (key1,key2)         
  case ("--vorticity")
    call Convert_vorticity()
  case ("--vor_abs")
    call Convert_abs_vorticity()    
  case ("--hdf2bin")
    call convert_hdf2bin()
  end select
      
  if (mpirank==0) write (*,*) "*** bye bye ***"   
end subroutine postprocessing




!-------------------------------------------------------------------------------
! ./flusi --postprocessing --hdf2bin ux_00000.h5
!-------------------------------------------------------------------------------
! converts the *.h5 file to an ordinairy binary file
subroutine convert_hdf2bin()
  use fsi_vars
  use mpi_header
  implicit none
  character(len=80) :: fname, dsetname  
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer, parameter :: pr_out = 4 
  integer :: ix, iy ,iz
  real(kind=pr_out), dimension(:,:,:), allocatable :: field_out ! single precision
  real(kind=pr) :: time 
  logical :: exist1
  call get_command_argument(3,fname)
  
  ! check if input file exists
  call check_file_exists ( fname )
  
  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return 
  endif    
  
  dsetname = fname ( 1:index( fname, '_' )-1 )
  call Fetch_attributes( fname, dsetname, nx, ny, nz, xl, yl, zl, time )
  
  write (*,'("Converting ",A," to ",A,".binary. Resolution is" 3(i4,1x))') &
        trim(fname), trim(fname), nx,ny,nz
  write (*,'("time=",es12.4," xl=",es12.4," yl=",es12.4," zl=",es12.4)') &
        time, xl, yl, zl
      
  allocate ( field(0:nx-1,0:ny-1,0:nz-1),field_out(0:nx-1,0:ny-1,0:nz-1) )
  ! read field from hdf file
  call Read_Single_File_serial (fname, field)
  ! convert to single precision
  field_out = real(field, kind=pr_out)
  
  write (*,'("maxval=",es12.4," minval=",es12.4)') &
        maxval(field_out),minval(field_out)
  
  ! dump binary file (this file will be called ux_00100.h5.binary)
  open (12, file = trim(fname)//".binary", form='unformatted', status='replace')
  write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
  close (12)
  
  deallocate (field, field_out) 
end subroutine convert_hdf2bin




!-------------------------------------------------------------------------------
! ./flusi --postprocessing --vor_abs ux_00000.h5 uy_00000.h5 uz_00000.h5
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! directly compute the absolute value of vorticity, do not save components
! can be done in parallel
subroutine Convert_abs_vorticity()
  use fsi_vars
  use mpi_header
  implicit none
  character(len=80) :: fname_ux, fname_uy, fname_uz, dsetname
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time 
  logical :: exist1,exist2,exist3
  
  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  
  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )
    
  if (mpirank == 0) then
    write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
  endif
  
  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
     write (*,*) "Error in arguments, files do not start with ux uy and uz"
     write (*,*) "note files have to be in the right order"
     stop
  endif
  
  
  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call Fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )
  
  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl  
    
  call fft_initialize() ! also initializes the domain decomp
  
  if (mpirank==0) write (*,*) "Done fft_initialize"
  
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  
  if (mpirank==0) write (*,*) "Allocated memory"
  
  call Read_Single_File ( fname_ux, u(:,:,:,1) )
  call Read_Single_File ( fname_uy, u(:,:,:,2) )
  call Read_Single_File ( fname_uz, u(:,:,:,3) )
    
  call FFT (uk(:,:,:,1),u(:,:,:,1))
  call FFT (uk(:,:,:,2),u(:,:,:,2))
  call FFT (uk(:,:,:,3),u(:,:,:,3))
  
  call curl_inplace(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  
  call IFFT (u(:,:,:,1),uk(:,:,:,1))
  call IFFT (u(:,:,:,2),uk(:,:,:,2))
  call IFFT (u(:,:,:,3),uk(:,:,:,3))
  
  ! now u contains the vorticity in physical space
  fname_ux='vor_abs'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
  
  if (mpirank == 0) then
    write (*,'("Writing to file",A)') trim(fname_ux)
  endif
  
  ! compute absolute vorticity:
  u(:,:,:,1) = dsqrt(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)
    
  call Save_Field_HDF5 ( time,fname_ux,u(:,:,:,1),"vor_abs")
  
  deallocate (u)
  deallocate (uk)
  call fft_free()
  
end subroutine Convert_abs_vorticity



!-------------------------------------------------------------------------------
! ./flusi --postprocessing --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! can be done in parallel
subroutine Convert_vorticity()
  use fsi_vars
  use mpi_header
  implicit none
  character(len=80) :: fname_ux, fname_uy, fname_uz, dsetname
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time 
  logical :: exist1,exist2,exist3
  
  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  
  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  
  if (mpirank == 0) then
    write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
  endif

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
     write (*,*) "Error in arguments, files do not start with ux uy and uz"
     write (*,*) "note files have to be in the right order"
     stop
  endif
  
  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call Fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )
  
  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl  
    
  call fft_initialize() ! also initializes the domain decomp
  
  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  
  call Read_Single_File ( fname_ux, u(:,:,:,1) )
  call Read_Single_File ( fname_uy, u(:,:,:,2) )
  call Read_Single_File ( fname_uz, u(:,:,:,3) )
  
  call FFT (uk(:,:,:,1),u(:,:,:,1))
  call FFT (uk(:,:,:,2),u(:,:,:,2))
  call FFT (uk(:,:,:,3),u(:,:,:,3))
  
  call curl_inplace(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  
  call IFFT (u(:,:,:,1),uk(:,:,:,1))
  call IFFT (u(:,:,:,2),uk(:,:,:,2))
  call IFFT (u(:,:,:,3),uk(:,:,:,3))
  
  ! now u contains the vorticity in physical space
  fname_ux='vorx'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
  fname_uy='vory'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
  fname_uz='vorz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)
  
  if (mpirank == 0) then
    write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
  endif
    
  call Save_Field_HDF5 ( time,fname_ux,u(:,:,:,1),"vorx")
  call Save_Field_HDF5 ( time,fname_uy,u(:,:,:,2),"vory")
  call Save_Field_HDF5 ( time,fname_uz,u(:,:,:,3),"vorz")
  
  deallocate (u)
  deallocate (uk)
  call fft_free()
  
end subroutine Convert_vorticity



!-------------------------------------------------------------------------------
! ./flusi --postprocessing --keyvalues mask_00000.h5
!-------------------------------------------------------------------------------
! load the specified *.h5 file and creates a *.key file that contains
! min / max / mean / L2 norm of the field data. This is used for unit testing
! so that we don't need to store entire fields but rather the *.key only
subroutine Keyvalues(filename)
  use fsi_vars
  use mpi_header
  implicit none
  character(len=*), intent(in) :: filename
  character(len=80) :: dsetname
  integer :: nx_file,ny_file,nz_file
  real(kind=pr) :: xl_file, yl_file, zl_file, time
  real(kind=pr), dimension(:,:,:), allocatable :: field
  logical :: exist1
  
  if (mpisize>1) then
    write (*,*) "--keyvalues is currently a serial version only, run it on 1CPU"
    stop 
  endif  
  
  call check_file_exists( filename )
  
  write (*,*) "analyzing file "//trim(adjustl(filename))//" for keyvalues"  
  
  !---------------------------------------------------------
  ! in a first step, we fetch the attributes from the dataset
  ! namely the resolution is whats important
  ! this routine was created in the mpi2vis repo -> convert_hdf2xmf.f90
  !---------------------------------------------------------
  ! the dataset is named the same way as the file:
  dsetname = filename ( 1:index( filename, '_' )-1 )
  write (*,*) "dsetname: ", dsetname
  call Fetch_attributes( filename, dsetname, nx, ny, nz, xl, yl, zl, time )
  allocate ( field(0:nx-1,0:ny-1,0:nz-1) )
  
  call Read_Single_File_serial (filename, field)
  
  open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
  write (14,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
  write (*,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
  close (14)  
  
  deallocate (field)  
end subroutine Keyvalues


!-------------------------------------------------------------------------------
! ./flusi --postprocessing --compare-keys mask_00000.key saved.key
!-------------------------------------------------------------------------------
! compares to *.key files if they're equal
subroutine Compare_key(key1,key2)
  use fsi_vars
  use mpi_header
  implicit none
  character(len=*), intent(in) :: key1,key2
  real(kind=pr) :: a1,a2,b1,b2,c1,c2,d1,d2
  logical :: exist1, exist2
  
  inquire ( file=key1, exist=exist1 )
  inquire ( file=key2, exist=exist2 )
  
  if ( exist1.eqv..false. .or. exist2.eqv..false. ) then
    write (*,*) "Input file not found..."
    stop
  endif  
  
  open  (14, file = key1, status = 'unknown', action='read')
  read (14,'(4(es17.10,1x))') a1,b1,c1,d1  
  close (14)   
  
  open  (14, file = key2, status = 'unknown', action='read')
  read (14,'(4(es17.10,1x))') a2,b2,c2,d2  
  close (14)   
  
  write (*,'("present  : max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
        a1,b1,c1,d1
  write (*,'("reference: max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
        a2,b2,c2,d2
  write (*,'("differenc: max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
        a2-a1,b2-b1,c2-c1,d2-d1
  
  if ((a1-a2<1.d-12) .and. (b2-b1<1.d-12) .and. (c2-c1<1.e-12) .and. (d2-d1<1.d-12)) then
    ! all cool
    write (*,*) "SUCCES"
    call exit(0)            
  elseif ((a1-a2<1.d-09) .and. (b2-b1<1.d-09) .and. (c2-c1<1.e-09) .and. (d2-d1<1.d-09)) then
    ! not so cool but no catastrophy
    write (*,*) "WARNING"
    call exit(2)            
  else
    ! very bad
    write (*,*) "ERROR"
    call exit(1)
  endif 
end subroutine Compare_key