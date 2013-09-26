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
  end select
      
  if (mpirank==0) write (*,*) "*** bye bye ***"   
end subroutine postprocessing







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
  
  inquire ( file=fname_ux, exist=exist1 )
  inquire ( file=fname_uy, exist=exist2 )
  inquire ( file=fname_uz, exist=exist3 )
  
  if ( exist1.eqv..false. .or. exist2.eqv..false. .or. exist3.eqv..false. ) then
    write (*,*) "Input file not found..."
    return
  endif
  
  if (mpirank == 0) then
    write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
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
  
  inquire ( file=filename, exist=exist1 )
  
  if (exist1.eqv..false.) then
    write(*,*) "Input file not found..."
    return
  endif
  
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