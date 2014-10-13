!-------------------------------------------------------------------------------
! Wrapper for different postprocessing tools
!-------------------------------------------------------------------------------
subroutine postprocessing()
  use vars
  use mpi
  implicit none
  character(len=strlen)     :: postprocessing_mode, filename, key1,key2
  real(kind=pr) :: t1
  
  t1=MPI_wtime()
  
  if (mpirank==0) write (*,*) "*** FLUSI is running in postprocessing mode ***"
  
  ! the second argument tells us what to do with the file
  call get_command_argument(2,postprocessing_mode)
  ! it then depends on the second argument what follows
!   
!   !-----------------
!   ! check what to do
!   !-----------------     
  select case (postprocessing_mode)
  case ('cp')
    call copy_hdf5file()
!   case ("--keyvalues")
!     call get_command_argument(3,filename)
!     call keyvalues (filename)      
!   case ("--compare-keys")
!     call get_command_argument(3,key1)
!     call get_command_argument(4,key2)
!     call compare_key (key1,key2)     
!   case ("--compare-timeseries")
!     call compare_timeseries() 
!   case ("--vorticity")
!     call convert_vorticity()
!   case ("--vor_abs")
!     call convert_abs_vorticity()    
!   case ("--hdf2bin")
!     call convert_hdf2bin()
!   case ("--bin2hdf")
!     call convert_bin2hdf()
!   case ("--p2Q")
!     call pressure_to_Qcriterion()
!   case ("--extract-subset")
!     call extract_subset()
  end select
      
  if (mpirank==0) write(*,'("Elapsed time=",es12.4)') MPI_wtime()-t1    
end subroutine postprocessing


subroutine copy_hdf5file
  use vars
  implicit none
  
  character(len=strlen) :: file_src, file_dst, dset_src, dset_dst
  real(kind=pr),dimension(:,:,:),allocatable::field
  real(kind=pr)::time
  
  call get_command_argument(3,file_src)
  call get_command_argument(4,file_dst)
  
  dset_dst = file_dst(1:index(file_dst,'_')-1)
  dset_src = file_src(1:index(file_src,'_')-1)
  
  call check_file_exists ( file_src )
  call fetch_attributes( file_src, dset_src, nx, ny, nz, xl, yl, zl, time )
  allocate ( field(0:nx-1,0:ny-1,0:nz-1))
  call read_single_file_serial(file_src, field)
  ra=0
  rb=(/nx-1,ny-1,nz-1/)
  call save_field_hdf5(time,file_dst(1:index(file_dst,'.')-1),field,dset_dst)
  deallocate(field)
end subroutine

! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --hdf2bin ux_00000.h5 filename.bin
! !-------------------------------------------------------------------------------
! ! converts the *.h5 file to an ordinairy binary file
! subroutine convert_hdf2bin()
!   use vars
!   use mpi
!   use basic_operators
!   implicit none
!   character(len=strlen) :: fname, dsetname  ,fname_bin
!   real(kind=pr), dimension(:,:,:), allocatable :: field
!   integer, parameter :: pr_out = 4 
!   integer :: ix, iy ,iz
!   real(kind=pr_out), dimension(:,:,:), allocatable :: field_out ! single precision
!   real(kind=pr) :: time 
!   call get_command_argument(3,fname)
!   call get_command_argument(4,fname_bin)
!   
!   ! check if input file exists
!   call check_file_exists ( fname )
!   
!   if ( mpisize>1 ) then
!     write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
!     return 
!   endif    
!   
!   dsetname = fname ( 1:index( fname, '_' )-1 )
!   call fetch_attributes( fname, dsetname, nx, ny, nz, xl, yl, zl, time )
!   
!   write (*,'("Converting ",A," to ",A," Resolution is",3(i4,1x))') &
!         trim(fname), trim(fname_bin), nx,ny,nz
!   write (*,'("time=",es12.4," xl=",es12.4," yl=",es12.4," zl=",es12.4)') &
!         time, xl, yl, zl
!       
!   allocate ( field(0:nx-1,0:ny-1,0:nz-1),field_out(0:nx-1,0:ny-1,0:nz-1) )
!   ! read field from hdf file
!   call read_single_file_serial (fname, field)
!   ! convert to single precision
!   field_out = real(field, kind=pr_out)
!   
!   write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field_out),minval(field_out)
!   
!   ! dump binary file (this file will be called ux_00100.h5.binary)
!   open (12, file = trim(fname_bin), form='unformatted', status='replace',&
!       convert="little_endian")
! !   write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
!   write(12) field_out
!   close (12)
!   
!   deallocate (field, field_out) 
! end subroutine convert_hdf2bin
! 
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --bin2hdf [file_bin] [file_hdf5] [nx] [ny] [nz] [xl] [yl] [zl] [time]
! ! ./flusi --postprocessing --bin2hdf ux_file.binary ux_00000.h5 128 128 384 3.5 2.5 10.0 0.0
! !-------------------------------------------------------------------------------
! ! converts the given binary file into an HDF5 file following flusi's conventions
! subroutine convert_bin2hdf()
!   use vars
!   use mpi
!   use basic_operators
!   implicit none
!   character(len=strlen) :: fname_bin,fname_hdf,dsetname,tmp  
!   real, dimension(:,:,:), allocatable :: field
!   integer, parameter :: pr_out = 4 
!   integer :: ix, iy ,iz, record_length
!   real(kind=pr) :: time 
!   
!   if ( mpisize>1 ) then
!     write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
!     return 
!   endif 
!   
!   
!   ! binary file name
!   call get_command_argument(3,fname_bin)
!   ! hdf5 file name
!   call get_command_argument(4,fname_hdf)
!   ! name of field in hdf5 file
!   dsetname=fname_hdf ( 1:index( fname_hdf, '_' )-1 )
!   call get_command_argument(5,tmp)
!   read (tmp,*) nx
!   call get_command_argument(6,tmp)
!   read (tmp,*) ny
!   call get_command_argument(7,tmp)
!   read (tmp,*) nz
!   call get_command_argument(8,tmp)
!   read (tmp,*) xl
!   call get_command_argument(9,tmp)
!   read (tmp,*) yl
!   call get_command_argument(10,tmp)
!   read (tmp,*) zl
!   call get_command_argument(11,tmp)
!   read (tmp,*) time
! 
!   write(*,'("converting ",A," into ",A," resolution: ",3(i4,1x)," box size: ",&
!        &3(es15.8,1x)," time=",es15.8)') trim(adjustl(fname_bin)), &
!        trim(adjustl(fname_hdf)), nx,ny,nz, xl,yl,zl,time
!   
!   !-----------------------------------------------------------------------------
!   ! read in the binary field to be converted
!   !-----------------------------------------------------------------------------
!   allocate ( field(0:nx-1,0:ny-1,0:nz-1) )
!   inquire (iolength=record_length) field
!   open(11, file=fname_bin, form='unformatted', &
!   access='direct', recl=record_length, convert="little_endian")
!   read (11,rec=1) field
!   close (11)  
!   write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)
! 
!   !-----------------------------------------------------------------------------
!   ! write the field data to an HDF file
!   !-----------------------------------------------------------------------------
!   ! initializes serial domain decomposition:
!   ra=(/0, 0, 0/)
!   rb=(/nx-1, ny-1, nz-1/)
!     
!   call save_field_hdf5(time,trim(adjustl(fname_hdf)),dble(field),trim(adjustl(dsetname))) 
!   
!   deallocate (field) 
! end subroutine convert_bin2hdf
! 
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --vor_abs ux_00000.h5 uy_00000.h5 uz_00000.h5
! !-------------------------------------------------------------------------------
! ! load the velocity components from file and compute & save the vorticity
! ! directly compute the absolute value of vorticity, do not save components
! ! can be done in parallel
! subroutine convert_abs_vorticity()
!   use vars
!   use p3dfft_wrapper
!   use mpi
!   use basic_operators
!   implicit none
!   character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname
!   complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
!   real(kind=pr),dimension(:,:,:,:),allocatable :: u
!   real(kind=pr) :: time 
!   
!   call get_command_argument(3,fname_ux)
!   call get_command_argument(4,fname_uy)
!   call get_command_argument(5,fname_uz)
!   
!   call check_file_exists(fname_ux)
!   call check_file_exists(fname_uy)
!   call check_file_exists(fname_uz)
!     
!   if (mpirank == 0) then
!     write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
!   endif
!   
!   if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
!      write (*,*) "Error in arguments, files do not start with ux uy and uz"
!      write (*,*) "note files have to be in the right order"
!      call abort()
!   endif
!   
!   
!   dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
!   call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )
!   
!   pi=4.d0 *datan(1.d0)
!   scalex=2.d0*pi/xl
!   scaley=2.d0*pi/yl
!   scalez=2.d0*pi/zl  
!   dx = xl/dble(nx)
!   dy = yl/dble(ny)
!   dz = zl/dble(nz)
!     
!   call fft_initialize() ! also initializes the domain decomp
!   
!   if (mpirank==0) write (*,*) "Done fft_initialize"
!   
!   allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
!   allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
!   
!   if (mpirank==0) write (*,*) "Allocated memory"
!   
!   call read_single_file ( fname_ux, u(:,:,:,1) )
!   call read_single_file ( fname_uy, u(:,:,:,2) )
!   call read_single_file ( fname_uz, u(:,:,:,3) )
!     
!   call fft (uk(:,:,:,1),u(:,:,:,1))
!   call fft (uk(:,:,:,2),u(:,:,:,2))
!   call fft (uk(:,:,:,3),u(:,:,:,3))
!   
!   call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
!   
!   call ifft (u(:,:,:,1),uk(:,:,:,1))
!   call ifft (u(:,:,:,2),uk(:,:,:,2))
!   call ifft (u(:,:,:,3),uk(:,:,:,3))
!   
!   ! now u contains the vorticity in physical space
!   fname_ux='vor_abs'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
!   
!   if (mpirank == 0) then
!     write (*,'("Writing to file",A)') trim(fname_ux)
!   endif
!   
!   ! compute absolute vorticity:
!   u(:,:,:,1) = dsqrt(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)
!     
!   call save_field_hdf5 ( time,fname_ux,u(:,:,:,1),"vor_abs")
!   
!   deallocate (u)
!   deallocate (uk)
!   call fft_free()
!   
! end subroutine convert_abs_vorticity
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --second-order
! !-------------------------------------------------------------------------------
! ! load the velocity components from file and compute & save the vorticity
! ! can be done in parallel. the flag --second order can be used for filtering
! subroutine convert_vorticity()
!   use vars
!   use p3dfft_wrapper
!   use basic_operators
!   use mpi
!   implicit none
!   character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, order
!   complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
!   real(kind=pr),dimension(:,:,:,:),allocatable :: u
!   real(kind=pr) :: time
!   
!   call get_command_argument(3,fname_ux)
!   call get_command_argument(4,fname_uy)
!   call get_command_argument(5,fname_uz)
!   call get_command_argument(6,order)
!   
!   call check_file_exists( fname_ux )
!   call check_file_exists( fname_uy )
!   call check_file_exists( fname_uz )
! 
!   if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
!      write (*,*) "Error in arguments, files do not start with ux uy and uz"
!      write (*,*) "note files have to be in the right order"
!      call abort()
!   endif
!   
!   dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
!   call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )
!   
!   pi=4.d0 *datan(1.d0)
!   scalex=2.d0*pi/xl
!   scaley=2.d0*pi/yl
!   scalez=2.d0*pi/zl 
!   dx = xl/dble(nx)
!   dy = yl/dble(ny)
!   dz = zl/dble(nz)
!     
!   call fft_initialize() ! also initializes the domain decomp
!   
!   allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
!   allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
!   
!   call read_single_file ( fname_ux, u(:,:,:,1) )
!   call read_single_file ( fname_uy, u(:,:,:,2) )
!   call read_single_file ( fname_uz, u(:,:,:,3) )
!   
!   call fft (uk(:,:,:,1),u(:,:,:,1))
!   call fft (uk(:,:,:,2),u(:,:,:,2))
!   call fft (uk(:,:,:,3),u(:,:,:,3))
!   
!   if (order=="--second-order") then
!     call curl_2nd(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3)) 
!   else
!     call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
!   endif
!   
!   call ifft (u(:,:,:,1),uk(:,:,:,1))
!   call ifft (u(:,:,:,2),uk(:,:,:,2))
!   call ifft (u(:,:,:,3),uk(:,:,:,3))
!   
!   ! now u contains the vorticity in physical space
!   fname_ux='vorx'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
!   fname_uy='vory'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
!   fname_uz='vorz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)
!   
!     
!   call save_field_hdf5 ( time,fname_ux,u(:,:,:,1),"vorx")
!   if (mpirank==0) write(*,*) "Wrote vorx to "//trim(fname_ux)
!   call save_field_hdf5 ( time,fname_uy,u(:,:,:,2),"vory")
!   if (mpirank==0) write(*,*) "Wrote vory to "//trim(fname_uy)
!   call save_field_hdf5 ( time,fname_uz,u(:,:,:,3),"vorz")
!   if (mpirank==0) write(*,*) "Wrote vorz to "//trim(fname_uz)
!   
!  
!   
!   deallocate (u)
!   deallocate (uk)
!   call fft_free()
!   
! end subroutine convert_vorticity
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --keyvalues mask_00000.h5
! !-------------------------------------------------------------------------------
! ! load the specified *.h5 file and creates a *.key file that contains
! ! min / max / mean / L2 norm of the field data. This is used for unit testing
! ! so that we don't need to store entire fields but rather the *.key only
! subroutine keyvalues(filename)
!   use vars
!   use mpi
!   implicit none
!   character(len=*), intent(in) :: filename
!   character(len=strlen) :: dsetname
!   real(kind=pr) :: time
!   real(kind=pr), dimension(:,:,:), allocatable :: field
!   
!   if (mpisize>1) then
!     write (*,*) "--keyvalues is currently a serial version only, run it on 1CPU"
!     call abort() 
!   endif  
!   
!   call check_file_exists( filename )
!   
!   write (*,*) "analyzing file "//trim(adjustl(filename))//" for keyvalues"  
!   
!   !---------------------------------------------------------
!   ! in a first step, we fetch the attributes from the dataset
!   ! namely the resolution is whats important
!   ! this routine was created in the mpi2vis repo -> convert_hdf2xmf.f90
!   !---------------------------------------------------------
!   ! the dataset is named the same way as the file:
!   dsetname = filename ( 1:index( filename, '_' )-1 )
!   call fetch_attributes( filename, dsetname, nx, ny, nz, xl, yl, zl, time )
!   allocate ( field(0:nx-1,0:ny-1,0:nz-1) )
!   
!   call read_single_file_serial (filename, field)
!   
!   open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
!   write (14,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
!   write (*,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
!   close (14)  
!   
!   deallocate (field)  
! end subroutine keyvalues
! 
! 
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --compare-timeseries forces.t ref/forces.t 
! !-------------------------------------------------------------------------------
! subroutine compare_timeseries()
!   use vars
!   implicit none
!   character(len=strlen) :: file1,file2
!   character(len=1024) :: header, line
!   character(len=15) ::format
!   real(kind=pr),dimension(:),allocatable :: values1, values2, error
!   real(kind=pr)::diff
!   integer :: i,columns,io_error,columns2
!   
!   call get_command_argument(3,file1)
!   call get_command_argument(4,file2)
!   
!   call check_file_exists(file1)
!   call check_file_exists(file2)
!   
!   !-----------------------------------------------------------------------------
!   ! how many colums are in the *.t file?
!   !-----------------------------------------------------------------------------
!   open  (14, file = file1, status = 'unknown', action='read')
!   read (14,'(A)') header
!   read (14,'(A)') line
!   columns=1
!   do i=2,len_trim(line)
!     if ((line(i:i)==" ").and.(line(i+1:i+1)/=" ")) columns=columns+1
!   enddo
!   close (14) 
!   
!   !-----------------------------------------------------------------------------
!   ! how many colums are in the second *.t file?
!   !-----------------------------------------------------------------------------
!   open  (14, file = file1, status = 'unknown', action='read')
!   read (14,'(A)') header
!   read (14,'(A)') line
!   columns2=1
!   do i=2,len_trim(line)
!     if ((line(i:i)==" ").and.(line(i+1:i+1)/=" ")) columns2=columns2+1
!   enddo
!   close (14)   
!   
!   if(columns/=columns2) then
!     write(*,*) "trying to compare two t files with different #columns..."
!     call exit(666)
!   endif
! 
!   write(format,'("(",i2.2,"(es15.8,1x))")') columns
! !   write(*,*) format
!   !-----------------------------------------------------------------------------
!   ! alloc arrays, then scan line by line for errors
!   !-----------------------------------------------------------------------------
!   allocate( values1(1:columns), values2(1:columns), error(1:columns) )
!   
!   ! read in files, line by line
!   io_error=0
!   open (20, file = file1, status = 'unknown', action='read')
!   open (30, file = file2, status = 'unknown', action='read')
!   read (20,'(A)') header
!   read (30,'(A)') header
!   do while (io_error==0)
!     ! compare this line
!     read (20,*,iostat=io_error) values1
!     read (30,*,iostat=io_error) values2
!     
!     do i=1,columns
!       diff = values1(i)-values2(i)
!       ! ignore values smaller 1e-4 in ref file
!       if ((dabs(values2(i))>1.d-4).and.(diff>1.d-7)) then
!         error(i) = dabs(diff/values2(i))
!       else
!         error(i) = 0.0
!       endif
!     enddo
!     
!     if (maxval(error)>1.d-4) then 
!       write(*,*) "time series comparison failed..."
!       write(*,format) values1
!       write(*,format) values2
!       write(*,format) error
!       call exit(666)
!     endif
!   enddo
!   close (20)
!   close (30)
!   deallocate (values1,values2,error)
! end subroutine compare_timeseries
! 
! !-------------------------------------------------------------------------------
! ! ./flusi --postprocessing --compare-keys mask_00000.key saved.key
! !-------------------------------------------------------------------------------
! ! compares to *.key files if they're equal
! subroutine compare_key(key1,key2)
!   use vars
!   use mpi
!   implicit none
!   character(len=*), intent(in) :: key1,key2
!   real(kind=pr) :: a1,a2,b1,b2,c1,c2,d1,d2
!   real(kind=pr) :: e1,e2,e3,e4
!   
!   call check_file_exists(key1)
!   call check_file_exists(key2) 
!   
!   open  (14, file = key1, status = 'unknown', action='read')
!   read (14,'(4(es17.10,1x))') a1,b1,c1,d1  
!   close (14)   
!   
!   open  (14, file = key2, status = 'unknown', action='read')
!   read (14,'(4(es17.10,1x))') a2,b2,c2,d2  
!   close (14)   
!   
!   write (*,'("present  : max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
!         a1,b1,c1,d1
!               
!   write (*,'("reference: max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
!         a2,b2,c2,d2
!         
!   ! errors:
!   if (dabs(a2)>=1.0d-7) then
!     e1 = dabs( (a2-a1) / a2 )
!   else
!     e1 = dabs( (a2-a1) )
!   endif  
!   
!   if (dabs(b2)>=1.0d-7) then
!     e2 = dabs( (b2-b1) / b2 )
!   else
!     e2 = dabs( (b2-b1) )
!   endif
!   
!   if (dabs(c2)>=1.0d-7) then
!     e3 = dabs( (c2-c1) / c2 )
!   else
!     e3 = dabs( (c2-c1) ) 
!   endif
!   
!   if (dabs(d2)>=1.0d-7) then
!     e4 = dabs( (d2-d1) / d2 )
!   else
!     e4 = dabs( (d2-d1) )
!   endif
!   
!   write (*,'("err(rel) : max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
!         e1,e2,e3,e4
!   
!   if ((e1<1.d-4) .and. (e2<1.d-4) .and. (e3<1.d-4) .and. (e4<1.d-4)) then
!     ! all cool
!     write (*,*) "OKAY..."
!     call exit(0)            
!   else
!     ! very bad
!     write (*,*) "ERROR"
!     call exit(1)
!   endif 
! end subroutine compare_key
! 
! 
! 
! !-------------------------------------------------------------------------------
! ! pressure_to_Qcriterion()
! ! converts a given pressure field and outputs the Q-criterion, computed with 
! ! second order and periodic boundary conditions. Alternatively, one 
! ! may use distint postprocessing tools, such as paraview, and compute the Q-crit
! ! there, but the accuracy may be different. 
! ! note the precision is reduced to second order (by using the effective 
! ! wavenumber), since spurious oscillations appear when computing it with
! ! spectral precision
! !-------------------------------------------------------------------------------
! ! call:
! ! ./flusi --postprocessing --p2Q p_00000.h5 Q_00000.h5
! !-------------------------------------------------------------------------------
! subroutine pressure_to_Qcriterion()
!   use mpi
!   use vars
!   use basic_operators
!   use p3dfft_wrapper
!   implicit none
!   character(len=strlen) :: fname_p, fname_Q
!   complex(kind=pr),dimension(:,:,:),allocatable :: pk
!   real(kind=pr),dimension(:,:,:),allocatable :: p
!   real(kind=pr)::time,maxi,mini
!   
!   ! get file to read pressure from and check if this is present
!   call get_command_argument(3,fname_p)
!   call check_file_exists( fname_p )
!   ! get filename to save Q criterion to
!   call get_command_argument(4,fname_Q)
! 
!   ! read in information from the file
!   call fetch_attributes( fname_p, "p", nx, ny, nz, xl, yl, zl, time )
!   
!   if (mpirank==0) then
!     write(*,'("Computing Q criterion from  file ",A," saving to ",&
!     & A," nx=",i4," ny=",i4," nz=",i4, &
!     &"xl=",es12.4," yl=",es12.4," zl=",es12.4 )') &
!     trim(fname_p), trim(fname_Q), nx,nx,nz,xl,yl,zl
!   endif
!   
!   pi=4.d0 * datan(1.d0)
!   scalex=2.d0*pi/xl
!   scaley=2.d0*pi/yl
!   scalez=2.d0*pi/zl  
!   
!   dx = xl/dble(nx)
!   dy = yl/dble(ny)
!   dz = zl/dble(nz)
!   
!   call fft_initialize() ! also initializes the domain decomp
!   
!   allocate(p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
!   allocate(pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
!   
!   call read_single_file(fname_p, p)
!   call fft(inx=p, outk=pk)
!   ! Q-criterion is 0.5*laplace(P)
!   ! see http://books.google.de/books?id=FWmsrNv3BYoC&pg=PA23&lpg=PA23&dq=q-criterion&source=bl&ots=CW-1NPn9p0&sig=Zi_Z2iw-ZuDqJqYctM9OUrb5WMA&hl=de&sa=X&ei=ZBCXU9nHGInVPPTRgagO&ved=0CCgQ6AEwADgK#v=onepage&q=q-criterion&f=false
!   call laplacien_inplace_filtered(pk)
!   call ifft(ink=pk, outx=p)
!   p=0.5d0*p
!   
!   call save_field_hdf5(time, trim(fname_Q(1:index(fname_Q,'.h5')-1)), p, "Q")
!   
!   maxi = fieldmax(P)
!   mini = fieldmin(P)
!   
!   if (mpirank==0) then
!     write(*,'("Q-criterion 2nd order maxval=",es12.4," minval=",es12.4)') maxi,mini
!   endif
!   
!   deallocate(p,pk)
!   call fft_free()
! end subroutine pressure_to_Qcriterion
! 
! 
! !-------------------------------------------------------------------------------
! ! extract subset
! ! loads a file to memory and extracts a subset, writing to a different file. We
! ! assume here that you do this for visualization; in this case, one usually keeps
! ! the original, larger files. For simplicity, ensure all files follow FLUSI 
! ! naming convention, it is thus recommended to use a subfolder.
! !-------------------------------------------------------------------------------
! ! call:
! ! ./flusi --postprocessing --extract-subset ux_00000.h5 resized/ux_00000.h5 128:1:256 128:2:1024 1:1:end
! !-------------------------------------------------------------------------------
! subroutine extract_subset()
!   use mpi
!   use vars
!   implicit none
!   character(len=strlen) :: fname_in, fname_out, dsetname_in, dsetname_out
!   character(len=strlen) :: xset,yset,zset
!   real(kind=pr)::time
!   integer :: ix,iy,iz,i
!   ! reduced domain size
!   integer :: nx1,nx2, ny1,ny2, nz1,nz2, nxs,nys,nzs
!   ! sizes of the new array
!   integer :: nx_red, ny_red, nz_red, ix_red, iy_red, iz_red
!   ! reduced domain extends
!   real(kind=pr) :: xl1, yl1, zl1
!   ! original field
!   real(kind=pr), dimension(:,:,:), allocatable :: field_in
!   ! reduced field
!   real(kind=pr), dimension(:,:,:), allocatable :: field_out
! 
!   if (mpisize/=1) then
!     write(*,*) "./flusi --postprocessing --extract-subset is a SERIAL routine, use 1CPU only"
!     call abort()
!   endif
!   
!   ! get file to read pressure from and check if this is present
!   call get_command_argument(3,fname_in)
!   call check_file_exists( fname_in )
!   
!   ! get filename to save Q criterion to
!   call get_command_argument(4,fname_out)  
!   
!   dsetname_in = fname_in ( 1:index( fname_in, '_' )-1 )
!   dsetname_out = fname_out ( 1:index( fname_out, '_' )-1 )
!   
!   call fetch_attributes( fname_in, dsetname_in, nx, ny, nz, xl, yl, zl, time )
!   
!   call get_command_argument(5,xset)
!   call get_command_argument(6,yset)
!   call get_command_argument(7,zset)
!   
!   ! red in subset from command line. it is given in the form 
!   ! ixmin:xspacing:ixmax as a string.
!   read (xset(1:index(xset,':')-1) ,*) nx1
!   read (xset(index(xset,':',.true.)+1:len_trim(xset)),*) nx2
!   read (xset(index(xset,':')+1:index(xset,':',.true.)-1),*) nxs
!   
!   read (yset(1:index(yset,':')-1) ,*) ny1
!   read (yset(index(yset,':',.true.)+1:len_trim(yset)),*) ny2
!   read (yset(index(yset,':')+1:index(yset,':',.true.)-1),*) nys
!   
!   read (zset(1:index(zset,':')-1) ,*) nz1
!   read (zset(index(zset,':',.true.)+1:len_trim(zset)),*) nz2
!   read (zset(index(zset,':')+1:index(zset,':',.true.)-1),*) nzs
!   
!   
!   ! stop if subset exceeds array bounds
!   if ( nx1<0 .or. nx2>nx-1 .or. ny1<0 .or. ny2>ny-1 .or. nz1<0 .or. nz2>nz-1) then
!     write (*,*) "subset indices exceed array bounds....proceed, but correct mistake"
!     nx1 = max(nx1,0)
!     ny1 = max(ny1,0)
!     nz1 = max(nz1,0)
!     nx2 = min(nx-1,nx2)
!     ny2 = min(ny-1,ny2)
!     nz2 = min(nz-1,nz2)
!   endif
!   
!   
!   write(*,'("Cropping field from " &
!   &,"0:",i4," | 0:",i4," | 0:",i4,&
!   &"   to subset   "&
!   &,i4,":",i2,":",i4," | ",i4,":",i2,":",i4," | ",i4,":",i2,":",i4)')&
!   nx-1,ny-1,nz-1,nx1,nxs,nx2,ny1,nys,ny2,nz1,nzs,nz2
!   
!     
!   
!   allocate ( field_in(0:nx-1,0:ny-1,0:nz-1) ) 
!   call read_single_file_serial(fname_in,field_in)
!   
!   dx = xl/dble(nx)
!   dy = yl/dble(ny)
!   dz = zl/dble(nz)
!   
!   
!   !-----------------------------------------------------------------------------
!   ix = nx1
!   nx_red = 1
!   do while (ix<=nx2)
!     nx_red = nx_red + 1
!     ix = ix+nxs
!   enddo
!   ! we counted one too far
!   ix = ix-nxs
!   nx_red = nx_red -1
!   ! warn if spacing does not allow to take last point into account
!   if (ix /= nx2) then
!     write(*,'("Warning, with the spacing ",i2," you provided we extract "&
!     &,i4,":",i4, " rather then what you specified (last point is missing)" )') &
!     nxs,nx1,ix
!   endif
!   !-----------------------------------------------------------------------------
!   iy = ny1
!   ny_red = 1
!   do while (iy<=ny2)
!     ny_red = ny_red + 1
!     iy = iy+nys
!   enddo
!   ! we counted one too far
!   iy = iy-nys
!   ny_red = ny_red -1
!   ! warn if spacing does not allow to take last point into account
!   if (iy /= ny2) then
!     write(*,'("Warning, with the spacing ",i2," you provided we extract "&
!     &,i4,":",i4, " rather then what you specified (last point is missing)" )') &
!     nys,ny1,iy
!   endif
!   !-----------------------------------------------------------------------------
!   iz=nz1
!   nz_red = 1
!   do while (iz<=nz2)
!     iz = iz+nzs
!     nz_red = nz_red+1
!   enddo
!   ! we counted one too far
!   nz_red = nz_red-1
!   iz = iz-nzs
!   ! warn if spacing does not allow to take last point into account
!   if (iz /= nz2) then
!     write(*,'("Warning, with the spacing ",i2," you provided we extract "&
!     &,i4,":",i4, " rather then what you specified (last point is missing)" )') &
!     nzs,nz1,iz
!   endif
!   !-----------------------------------------------------------------------------
!   write (*,'("Size of subset is ",3(i4,1x))') nx_red, ny_red, nz_red
!   !-----------------------------------------------------------------------------
!   
!   ! we figured out how big the subset array is
!   allocate ( field_out(0:nx_red-1,0:ny_red-1,0:nz_red-1) )
!   
!   ! fill the target field
!   do ix_red = 0, nx_red-1
!     do iy_red = 0, ny_red-1
!       do iz_red = 0, nz_red-1
!         ix = nx1 + ix_red*nxs
!         iy = ny1 + iy_red*nys
!         iz = nz1 + iz_red*nzs
!         field_out(ix_red,iy_red,iz_red) = field_in(ix,iy,iz)
!       enddo
!     enddo
!   enddo
!   
!   write(*,*) maxval(field_out), maxval(field_in)
!   
!   ! set up dimensions in global variables, since save_field_hdf5 relies on this
!   ra = 0
!   rb(1) = nx_red-1
!   rb(2) = ny_red-1
!   rb(3) = nz_red-1
!   nx = nx_red
!   ny = ny_red
!   nz = nz_red  
!   xl = dx+dble(nx1 + (nx_red-1)*nxs)*dx - dble(nx1)*dx
!   yl = dy+dble(ny1 + (ny_red-1)*nys)*dy - dble(ny1)*dy
!   zl = dz+dble(nz1 + (nz_red-1)*nzs)*dz - dble(nz1)*dz
!   
!   
!   call save_field_hdf5 ( time, fname_out(1:index(fname_out,'.h5')-1), field_out, dsetname_out )
!   
!   deallocate (field_in, field_out)
! end subroutine extract_subset
