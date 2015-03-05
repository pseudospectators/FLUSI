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

  !-----------------
  ! check what to do
  !-----------------
  select case (postprocessing_mode)
  case ("--cp")
    call copy_hdf_file()
  case ("--keyvalues")
    call get_command_argument(3,filename)
    call keyvalues (filename)
  case ("--compare-keys")
    call get_command_argument(3,key1)
    call get_command_argument(4,key2)
    call compare_key (key1,key2)
  case ("--compare-timeseries")
    call compare_timeseries()
  case ("--vorticity")
    call convert_vorticity()
  case ("--vor2u")
    call convert_velocity()
  case ("--vor_abs")
    call convert_abs_vorticity()
  case ("--hdf2bin")
    call convert_hdf2bin()
  case ("--bin2hdf")
    call convert_bin2hdf()
  case ("--p2Q")
    call pressure_to_Qcriterion()
  case ("--extract-subset")
    call extract_subset()
  case ("--time-avg")
    call time_avg_HDF5()
  case ("--upsample")
    call upsample()
  case ("--spectrum")
    call post_spectrum()
  case ("--turbulence-analysis")
    call turbulence_analysis()
  case ("--TKE-mean")
    call tke_mean()
  case ("--max-over-x")
    call max_over_x()
  case ("--mean-over-x-subdomain")
    call mean_over_x_subdomain()
  case default
    if(mpirank==0) write(*,*) "Postprocessing option is "// postprocessing_mode
    if(mpirank==0) write(*,*) "But I don't know what to do with that"
  end select

  if (mpirank==0) write(*,'("Elapsed time=",es12.4)') MPI_wtime()-t1
end subroutine postprocessing






!-------------------------------------------------------------------------------
! ./flusi --postprocess --hdf2bin ux_00000.h5 filename.bin
!-------------------------------------------------------------------------------
! converts the *.h5 file to an ordinairy binary file
subroutine convert_hdf2bin()
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  character(len=strlen) :: fname, fname_bin
  real(kind=pr), dimension(:,:,:), allocatable :: field
  integer, parameter :: pr_out = 4
  integer :: ix, iy ,iz
  real(kind=pr_out), dimension(:,:,:), allocatable :: field_out ! single precision
  real(kind=pr) :: time
  call get_command_argument(3,fname)
  call get_command_argument(4,fname_bin)

  ! check if input file exists
  call check_file_exists ( fname )

  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return
  endif

  call fetch_attributes( fname, get_dsetname(fname), nx, ny, nz, xl, yl, zl, time )

  write (*,'("Converting ",A," to ",A," Resolution is",3(i4,1x))') &
        trim(fname), trim(fname_bin), nx,ny,nz
  write (*,'("time=",es12.4," xl=",es12.4," yl=",es12.4," zl=",es12.4)') &
        time, xl, yl, zl

  allocate ( field(0:nx-1,0:ny-1,0:nz-1),field_out(0:nx-1,0:ny-1,0:nz-1) )
  ! read field from hdf file
  call read_single_file_serial (fname, field)
  ! convert to single precision
  field_out = real(field, kind=pr_out)

  write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field_out),minval(field_out)

  ! dump binary file (this file will be called ux_00100.h5.binary)
  open (12, file = trim(fname_bin), form='unformatted', status='replace',&
      convert="little_endian")
  write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
!  write(12) field_out
  close (12)

  deallocate (field, field_out)
end subroutine convert_hdf2bin



!-------------------------------------------------------------------------------
! ./flusi --postprocess --time-avg file_list.txt avgx_0000.h5
! Reads in a list of files from a file, then loads one file after the other and
! computes the average field, which is then stored in the specified file.
!-------------------------------------------------------------------------------
subroutine time_avg_HDF5()
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  character(len=strlen) :: fname, fname_bin, fname_avg
  real(kind=pr), dimension(:,:,:), allocatable :: field_avg, field
  integer :: ix, iy ,iz, io_error=0, i=0
  real(kind=pr) :: time
  call get_command_argument(3,fname)
  call get_command_argument(4,fname_avg)

  !-----------------------------------------------------------------------------
  ! check if input file exists, the file contains the list of h5 files to be avg
  !-----------------------------------------------------------------------------
  call check_file_exists ( fname )
  write(*,*) "Reading list of files from "//fname

  if ( mpisize>1 ) then
    write (*,*) "--time-avg is currently a serial version only, run it on 1CPU"
    return
  endif

  !-----------------------------------------------------------------------------
  ! read in the file, loop over lines
  !-----------------------------------------------------------------------------
  open( unit=14,file=fname, action='read', status='old')
  do while (io_error==0)
    ! fetch current filename
    read (14,'(A)', iostat=io_error) fname_bin
    write(*,*) "read "//trim(adjustl(fname_bin))
    if (io_error == 0) then
        write(*,*) "Processing file "//trim(adjustl(fname_bin))

        call check_file_exists ( fname_bin )
        call fetch_attributes( fname_bin, get_dsetname(fname_bin), nx, ny, nz, &
             xl, yl, zl, time )

        ! first time? allocate then.
        if ( .not. allocated(field_avg) ) then
          ra=(/0,0,0/)
          rb=(/nx-1,ny-1,nz-1/)
          allocate(field_avg(0:nx-1,0:ny-1,0:nz-1))
          allocate(field(0:nx-1,0:ny-1,0:nz-1))
          field_avg = 0.d0
        endif

        ! read the field from file
        call read_single_file_serial( fname_bin, field )

        field_avg = field_avg + field

        i = i+1
    endif
  enddo
  close (14)

  field_avg = field_avg / dble(i)

  call save_field_hdf5(0.d0, fname_avg, field_avg, get_dsetname(fname_avg))

  deallocate(field_avg)
  deallocate(field)


end subroutine time_avg_HDF5



!-------------------------------------------------------------------------------
! ./flusi --postprocess --bin2hdf [file_bin] [file_hdf5] [nx] [ny] [nz] [xl] [yl] [zl] [time]
! ./flusi --postprocess --bin2hdf ux_file.binary ux_00000.h5 128 128 384 3.5 2.5 10.0 0.0
!-------------------------------------------------------------------------------
! converts the given binary file into an HDF5 file following flusi's conventions
subroutine convert_bin2hdf()
  use vars
  use mpi
  use basic_operators
  use helpers
  implicit none
  character(len=strlen) :: fname_bin,fname_hdf,tmp
  real, dimension(:,:,:), allocatable :: field
  integer, parameter :: pr_out = 4
  integer :: ix, iy ,iz, i,j,k
  integer(kind=8) :: record_length
  real(kind=pr) :: time

  if ( mpisize>1 ) then
    write (*,*) "--hdf2bin is currently a serial version only, run it on 1CPU"
    return
  endif


  ! binary file name
  call get_command_argument(3,fname_bin)
  ! hdf5 file name
  call get_command_argument(4,fname_hdf)
  call get_command_argument(5,tmp)
  read (tmp,*) nx
  call get_command_argument(6,tmp)
  read (tmp,*) ny
  call get_command_argument(7,tmp)
  read (tmp,*) nz
  call get_command_argument(8,tmp)
  read (tmp,*) xl
  call get_command_argument(9,tmp)
  read (tmp,*) yl
  call get_command_argument(10,tmp)
  read (tmp,*) zl
  call get_command_argument(11,tmp)
  read (tmp,*) time

  write(*,'("converting ",A," into ",A," resolution: ",3(i4,1x)," box size: ",&
       &3(es15.8,1x)," time=",es15.8)') trim(adjustl(fname_bin)), &
       trim(adjustl(fname_hdf)), nx,ny,nz, xl,yl,zl,time

  !-----------------------------------------------------------------------------
  ! read in the binary field to be converted
  !-----------------------------------------------------------------------------
  allocate ( field(0:nx-1,0:ny-1,0:nz-1) )

!   inquire (iolength=record_length) field
!   open(11, file=fname_bin, form='unformatted', &
!   access='direct', recl=record_length, convert="little_endian")
!   read (11,rec=1) field
!   close (11)
!   write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)


  OPEN(10,FILE=fname_bin,FORM='unformatted',STATUS='OLD', convert='LITTLE_ENDIAN')
  read(10) (((field(i,j,k),i=0,nx-1),j=0,ny-1),k=0,nz-1)
  CLOSE(10)
  write (*,'("maxval=",es12.4," minval=",es12.4)') maxval(field),minval(field)
  !-----------------------------------------------------------------------------
  ! write the field data to an HDF file
  !-----------------------------------------------------------------------------
  ! initializes serial domain decomposition:
  ra=(/0, 0, 0/)
  rb=(/nx-1, ny-1, nz-1/)

  call save_field_hdf5(time,trim(adjustl(fname_hdf)),dble(field),get_dsetname(fname_hdf))

  deallocate (field)
end subroutine convert_bin2hdf




!-------------------------------------------------------------------------------
! ./flusi --postprocess --vor_abs ux_00000.h5 uy_00000.h5 uz_00000.h5
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! directly compute the absolute value of vorticity, do not save components
! can be done in parallel
subroutine convert_abs_vorticity()
  use vars
  use p3dfft_wrapper
  use mpi
  use helpers
  use basic_operators
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)

  call check_file_exists(fname_ux)
  call check_file_exists(fname_uy)
  call check_file_exists(fname_uz)

  if (mpirank == 0) then
    write (*,'(3(A,","))') trim(fname_ux), trim(fname_uy), trim(fname_uz)
  endif

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
     write (*,*) "Error in arguments, files do not start with ux uy and uz"
     write (*,*) "note files have to be in the right order"
     call abort()
  endif


  call fetch_attributes( fname_ux, get_dsetname(fname_ux), nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  if (mpirank==0) write (*,*) "Done fft_initialize"

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  if (mpirank==0) write (*,*) "Allocated memory"

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft (uk(:,:,:,1),u(:,:,:,1))
  call fft (uk(:,:,:,2),u(:,:,:,2))
  call fft (uk(:,:,:,3),u(:,:,:,3))

  call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))

  call ifft (u(:,:,:,1),uk(:,:,:,1))
  call ifft (u(:,:,:,2),uk(:,:,:,2))
  call ifft (u(:,:,:,3),uk(:,:,:,3))

  ! now u contains the vorticity in physical space
  fname_ux='vor_abs'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)

  if (mpirank == 0) then
    write (*,'("Writing to file",A)') trim(fname_ux)
  endif

  ! compute absolute vorticity:
  u(:,:,:,1) = dsqrt(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)

  call save_field_hdf5 ( time,fname_ux,u(:,:,:,1),"vor_abs")

  deallocate (u)
  deallocate (uk)
  call fft_free()

end subroutine convert_abs_vorticity



!-------------------------------------------------------------------------------
! ./flusi --postprocess --vorticity ux_00000.h5 uy_00000.h5 uz_00000.h5 --second-order
!-------------------------------------------------------------------------------
! load the velocity components from file and compute & save the vorticity
! can be done in parallel. the flag --second order can be used for filtering
subroutine convert_vorticity()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, order
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,order)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
     write (*,*) "Error in arguments, files do not start with ux uy and uz"
     write (*,*) "note files have to be in the right order"
     call abort()
  endif

  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft (uk(:,:,:,1),u(:,:,:,1))
  call fft (uk(:,:,:,2),u(:,:,:,2))
  call fft (uk(:,:,:,3),u(:,:,:,3))

  if (order=="--second-order") then
    call curl_2nd(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  else
    call curl(uk(:,:,:,1),uk(:,:,:,2),uk(:,:,:,3))
  endif

  call ifft (u(:,:,:,1),uk(:,:,:,1))
  call ifft (u(:,:,:,2),uk(:,:,:,2))
  call ifft (u(:,:,:,3),uk(:,:,:,3))

  ! now u contains the vorticity in physical space
  fname_ux='vorx'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
  fname_uy='vory'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
  fname_uz='vorz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)


  call save_field_hdf5 ( time,fname_ux,u(:,:,:,1),"vorx")
  if (mpirank==0) write(*,*) "Wrote vorx to "//trim(fname_ux)
  call save_field_hdf5 ( time,fname_uy,u(:,:,:,2),"vory")
  if (mpirank==0) write(*,*) "Wrote vory to "//trim(fname_uy)
  call save_field_hdf5 ( time,fname_uz,u(:,:,:,3),"vorz")
  if (mpirank==0) write(*,*) "Wrote vorz to "//trim(fname_uz)



  deallocate (u)
  deallocate (uk)
  call fft_free()

end subroutine convert_vorticity



!-------------------------------------------------------------------------------
! ./flusi --postprocess --vor2u vorx_00000.h5 vory_00000.h5 vorz_00000.h5
!-------------------------------------------------------------------------------
subroutine convert_velocity()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, order
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk, workc
  real(kind=pr),dimension(:,:,:,:),allocatable :: u,workr
  real(kind=pr) :: time

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if ((fname_ux(1:4).ne."vorx").or.(fname_uy(1:4).ne."vory").or.(fname_uz(1:4).ne."vorz")) then
     write (*,*) "Error in arguments, files do not start with vorx vory and vorz"
     write (*,*) "note files have to be in the right order"
     call abort()
  endif

  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi = 4.d0 *datan(1.d0)
  scalex = 2.d0*pi/xl
  scaley = 2.d0*pi/yl
  scalez = 2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  neq=3
  nd=3
  ncw=3
  nrw=3

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(workr(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))

  ! read vorticity to u
  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  workr = u ! workr is original vorticity

  call Vorticity2Velocity(uk,workc,u)
  call ifft3 (ink=uk, outx=u)   ! u is velocity in phys space

  ! now u contains the vorticity in physical space
  fname_ux='ux'//fname_ux(index(fname_ux,'_'):index(fname_ux,'.')-1)
  fname_uy='uy'//fname_uy(index(fname_uy,'_'):index(fname_uy,'.')-1)
  fname_uz='uz'//fname_uz(index(fname_uz,'_'):index(fname_uz,'.')-1)

  call save_field_hdf5 ( time,fname_ux,u(:,:,:,1),"ux")
  if (mpirank==0) write(*,*) "Wrote ux to "//trim(fname_ux)
  call save_field_hdf5 ( time,fname_uy,u(:,:,:,2),"uy")
  if (mpirank==0) write(*,*) "Wrote uy to "//trim(fname_uy)
  call save_field_hdf5 ( time,fname_uz,u(:,:,:,3),"uz")
  if (mpirank==0) write(*,*) "Wrote uz to "//trim(fname_uz)

  call divergence(uk,workc(:,:,:,1))
  call ifft(ink=workc(:,:,:,1),outx=u(:,:,:,1))
  !call save_field_hdf5 ( time,"divu_0000",u(:,:,:,1),"divu")
  write(*,*) "maximum divergence=", fieldmax(u(:,:,:,1)), fieldmin(u(:,:,:,1))

  call curl3_inplace(uk)
  call ifft3(ink=uk, outx=u)

  ! difference
  u(:,:,:,1)=u(:,:,:,1)-workr(:,:,:,1)
  u(:,:,:,2)=u(:,:,:,2)-workr(:,:,:,2)
  u(:,:,:,3)=u(:,:,:,3)-workr(:,:,:,3)

  write(*,*) "max diff vor_original-curl(result)", fieldmax(u(:,:,:,1))
  write(*,*) "max diff vor_original-curl(result)", fieldmax(u(:,:,:,2))
  write(*,*) "max diff vor_original-curl(result)", fieldmax(u(:,:,:,3))

  deallocate (u)
  deallocate (uk)
  deallocate (workc,workr)
  call fft_free()

end subroutine convert_velocity




!-------------------------------------------------------------------------------
! ./flusi --postprocess --keyvalues mask_00000.h5
!-------------------------------------------------------------------------------
! load the specified *.h5 file and creates a *.key file that contains
! min / max / mean / L2 norm of the field data. This is used for unit testing
! so that we don't need to store entire fields but rather the *.key only
subroutine keyvalues(filename)
  use vars
  use mpi
  implicit none
  character(len=*), intent(in) :: filename
  character(len=strlen) :: dsetname
  real(kind=pr) :: time
  real(kind=pr), dimension(:,:,:), allocatable :: field

  if (mpisize>1) then
    write (*,*) "--keyvalues is currently a serial version only, run it on 1CPU"
    call abort()
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
  call fetch_attributes( filename, dsetname, nx, ny, nz, xl, yl, zl, time )
  write(*,'("File is at time=",es12.4)') time
  allocate ( field(0:nx-1,0:ny-1,0:nz-1) )

  call read_single_file_serial (filename, field)

  open  (14, file = filename(1:index(filename,'.'))//'key', status = 'replace')
  write (14,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
  write (*,'(4(es17.10,1x))') maxval(field), minval(field), sum(field)/(nx*ny*nz), sum(field**2)/(nx*ny*nz)
  close (14)

  deallocate (field)
end subroutine keyvalues





!-------------------------------------------------------------------------------
! ./flusi --postprocess --compare-timeseries forces.t ref/forces.t
!-------------------------------------------------------------------------------
subroutine compare_timeseries()
  use fsi_vars
  implicit none
  character(len=strlen) :: file1,file2
  character(len=1024) :: header, line
  character(len=15) ::format
  real(kind=pr),dimension(:),allocatable :: values1, values2, error
  real(kind=pr)::diff
  integer :: i,columns,io_error,columns2

  call get_command_argument(3,file1)
  call get_command_argument(4,file2)

  call check_file_exists(file1)
  call check_file_exists(file2)

  !-----------------------------------------------------------------------------
  ! how many colums are in the *.t file?
  !-----------------------------------------------------------------------------
  open  (14, file = file1, status = 'unknown', action='read')
  read (14,'(A)') header
  read (14,'(A)') line
  columns=1
  do i=2,len_trim(line)
    if ((line(i:i)==" ").and.(line(i+1:i+1)/=" ")) columns=columns+1
  enddo
  close (14)

  !-----------------------------------------------------------------------------
  ! how many colums are in the second *.t file?
  !-----------------------------------------------------------------------------
  open  (14, file = file1, status = 'unknown', action='read')
  read (14,'(A)') header
  read (14,'(A)') line
  columns2=1
  do i=2,len_trim(line)
    if ((line(i:i)==" ").and.(line(i+1:i+1)/=" ")) columns2=columns2+1
  enddo
  close (14)

  if(columns/=columns2) then
    write(*,*) "trying to compare two t files with different #columns..."
    call exit(666)
  endif

  write(format,'("(",i2.2,"(es15.8,1x))")') columns
!   write(*,*) format
  !-----------------------------------------------------------------------------
  ! alloc arrays, then scan line by line for errors
  !-----------------------------------------------------------------------------
  allocate( values1(1:columns), values2(1:columns), error(1:columns) )

  ! read in files, line by line
  io_error=0
  open (20, file = file1, status = 'unknown', action='read')
  open (30, file = file2, status = 'unknown', action='read')
  read (20,'(A)') header
  read (30,'(A)') header
  do while (io_error==0)
    ! compare this line
    read (20,*,iostat=io_error) values1
    read (30,*,iostat=io_error) values2

    do i=1,columns
      diff = values1(i)-values2(i)
      ! ignore values smaller 1e-4 in ref file
      if ((dabs(values2(i))>1.d-4).and.(diff>1.d-7)) then
        error(i) = dabs(diff/values2(i))
      else
        error(i) = 0.0
      endif
    enddo

    if (maxval(error)>1.d-4) then
      write(*,*) "time series comparison failed..."
      write(*,format) values1
      write(*,format) values2
      write(*,format) error
      call exit(666)
    endif
  enddo
  close (20)
  close (30)
  deallocate (values1,values2,error)
end subroutine compare_timeseries

!-------------------------------------------------------------------------------
! ./flusi --postprocess --compare-keys mask_00000.key saved.key
!-------------------------------------------------------------------------------
! compares to *.key files if they're equal
subroutine compare_key(key1,key2)
  use vars
  use mpi
  implicit none
  character(len=*), intent(in) :: key1,key2
  real(kind=pr) :: a1,a2,b1,b2,c1,c2,d1,d2
  real(kind=pr) :: e1,e2,e3,e4

  call check_file_exists(key1)
  call check_file_exists(key2)

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

  ! errors:
  if (dabs(a2)>=1.0d-7) then
    e1 = dabs( (a2-a1) / a2 )
  else
    e1 = dabs( (a2-a1) )
  endif

  if (dabs(b2)>=1.0d-7) then
    e2 = dabs( (b2-b1) / b2 )
  else
    e2 = dabs( (b2-b1) )
  endif

  if (dabs(c2)>=1.0d-7) then
    e3 = dabs( (c2-c1) / c2 )
  else
    e3 = dabs( (c2-c1) )
  endif

  if (dabs(d2)>=1.0d-7) then
    e4 = dabs( (d2-d1) / d2 )
  else
    e4 = dabs( (d2-d1) )
  endif

  write (*,'("err(rel) : max=",es17.10," min=",es17.10," sum=",es17.10," sum**2=",es17.10)') &
        e1,e2,e3,e4

  if ((e1<1.d-4) .and. (e2<1.d-4) .and. (e3<1.d-4) .and. (e4<1.d-4)) then
    ! all cool
    write (*,*) "OKAY..."
    call exit(0)
  else
    ! very bad
    write (*,*) "ERROR"
    call exit(1)
  endif
end subroutine compare_key



!-------------------------------------------------------------------------------
! pressure_to_Qcriterion()
! converts a given pressure field and outputs the Q-criterion, computed with
! second order and periodic boundary conditions. Alternatively, one
! may use distint postprocessing tools, such as paraview, and compute the Q-crit
! there, but the accuracy may be different.
! note the precision is reduced to second order (by using the effective
! wavenumber), since spurious oscillations appear when computing it with
! spectral precision
!-------------------------------------------------------------------------------
! call:
! ./flusi --postprocess --p2Q p_00000.h5 Q_00000.h5
!-------------------------------------------------------------------------------
subroutine pressure_to_Qcriterion()
  use mpi
  use vars
  use basic_operators
  use p3dfft_wrapper
  implicit none
  character(len=strlen) :: fname_p, fname_Q
  complex(kind=pr),dimension(:,:,:),allocatable :: pk
  real(kind=pr),dimension(:,:,:),allocatable :: p
  real(kind=pr)::time,maxi,mini

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_p)
  call check_file_exists( fname_p )
  ! get filename to save Q criterion to
  call get_command_argument(4,fname_Q)

  ! read in information from the file
  call fetch_attributes( fname_p, "p", nx, ny, nz, xl, yl, zl, time )

  if (mpirank==0) then
    write(*,'("Computing Q criterion from  file ",A," saving to ",&
    & A," nx=",i4," ny=",i4," nz=",i4, &
    &"xl=",es12.4," yl=",es12.4," zl=",es12.4 )') &
    trim(fname_p), trim(fname_Q), nx,nx,nz,xl,yl,zl
  endif

  pi=4.d0 * datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl

  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate(p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  call read_single_file(fname_p, p)
  call fft(inx=p, outk=pk)
  ! Q-criterion is 0.5*laplace(P)
  ! see http://books.google.de/books?id=FWmsrNv3BYoC&pg=PA23&lpg=PA23&dq=q-criterion&source=bl&ots=CW-1NPn9p0&sig=Zi_Z2iw-ZuDqJqYctM9OUrb5WMA&hl=de&sa=X&ei=ZBCXU9nHGInVPPTRgagO&ved=0CCgQ6AEwADgK#v=onepage&q=q-criterion&f=false
  call laplacien_inplace_filtered(pk)
  call ifft(ink=pk, outx=p)
  p=0.5d0*p

  call save_field_hdf5(time, trim(fname_Q(1:index(fname_Q,'.h5')-1)), p, "Q")

  maxi = fieldmax(P)
  mini = fieldmin(P)

  if (mpirank==0) then
    write(*,'("Q-criterion 2nd order maxval=",es12.4," minval=",es12.4)') maxi,mini
  endif

  deallocate(p,pk)
  call fft_free()
end subroutine pressure_to_Qcriterion


!-------------------------------------------------------------------------------
! extract subset
! loads a file to memory and extracts a subset, writing to a different file. We
! assume here that you do this for visualization; in this case, one usually keeps
! the original, larger files. For simplicity, ensure all files follow FLUSI
! naming convention.
!---
! Note: through using helpers.f90::get_dsetname, it is finally possible to write
! to a subfolder *.h5 file directly from flusi
! (Thomas, 03/2015)
!---
! Using HDF5s hyperslab functions, we can read only a specific part into the
! memory - at no point we have to load the entire original file before we can
! downsample it. This is a good step forward. (Thomas 03/2015)
!-------------------------------------------------------------------------------
! call:
! ./flusi --postprocess --extract-subset ux_00000.h5 sux_00000.h5 128:1:256 128:2:1024 1:1:9999
!-------------------------------------------------------------------------------
subroutine extract_subset()
  use mpi
  use vars
  use hdf5
  use basic_operators
  use helpers
  implicit none

  character(len=strlen) :: fname_in, fname_out, dsetname_in, dsetname_out
  character(len=strlen) :: xset,yset,zset
  integer :: ix,iy,iz,i
  ! reduced domain size
  integer :: nx1,nx2, ny1,ny2, nz1,nz2, nxs,nys,nzs
  ! sizes of the new array
  integer :: nx_red, ny_red, nz_red, ix_red, iy_red, iz_red
  ! reduced domain extends
  real(kind=pr) :: xl1, yl1, zl1
  real(kind=pr), dimension(:,:,:), allocatable :: field

  integer, parameter            :: rank = 3 ! data dimensionality (2D or 3D)
  real (kind=pr)                :: time, xl_file, yl_file, zl_file
  character(len=80)             :: dsetname
  integer                       :: nx_file, ny_file, nz_file, mpierror

  integer(hid_t) :: file_id       ! file identifier
  integer(hid_t) :: dset_id       ! dataset identifier
  integer(hid_t) :: filespace     ! dataspace identifier in file
  integer(hid_t) :: memspace      ! dataspace identifier in memory
  integer(hid_t) :: plist_id      ! property list identifier

  ! dataset dimensions in the file.
  integer(hsize_t), dimension(rank) :: dimensions_file
  integer(hsize_t), dimension(rank) :: dimensions_local  ! chunks dimensions
  integer(hsize_t), dimension(rank) :: chunking_dims  ! chunks dimensions

  integer(hsize_t),  dimension(rank) :: count  = 1
  integer(hssize_t), dimension(rank) :: offset
  integer(hsize_t),  dimension(rank) :: stride = 1
  integer :: error  ! error flags

  ! what follows is for the attribute "time"
  integer, parameter :: arank = 1


  if (mpisize/=1) then
    write(*,*) "./flusi --postprocess --extract-subset is a SERIAL routine, use 1CPU only"
    call abort()
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  ! get filename to save subset to
  call get_command_argument(4,fname_out)

  dsetname_in  = get_dsetname( fname_in )
  dsetname_out = get_dsetname( fname_out )

  write(*,'("dsetname=",A,1x,A)') trim(adjustl(dsetname_in)),trim(adjustl(dsetname_out))

  call fetch_attributes( fname_in, dsetname_in, nx, ny, nz, xl, yl, zl, time )

  call get_command_argument(5,xset)
  call get_command_argument(6,yset)
  call get_command_argument(7,zset)

  ! red in subset from command line. it is given in the form
  ! ixmin:xspacing:ixmax as a string.
  read (xset(1:index(xset,':')-1) ,*) nx1
  read (xset(index(xset,':',.true.)+1:len_trim(xset)),*) nx2
  read (xset(index(xset,':')+1:index(xset,':',.true.)-1),*) nxs

  read (yset(1:index(yset,':')-1) ,*) ny1
  read (yset(index(yset,':',.true.)+1:len_trim(yset)),*) ny2
  read (yset(index(yset,':')+1:index(yset,':',.true.)-1),*) nys

  read (zset(1:index(zset,':')-1) ,*) nz1
  read (zset(index(zset,':',.true.)+1:len_trim(zset)),*) nz2
  read (zset(index(zset,':')+1:index(zset,':',.true.)-1),*) nzs


  ! stop if subset exceeds array bounds
  if ( nx1<0 .or. nx2>nx-1 .or. ny1<0 .or. ny2>ny-1 .or. nz1<0 .or. nz2>nz-1) then
    write (*,*) "subset indices exceed array bounds....proceed, but correct mistake"
    nx1 = max(nx1,0)
    ny1 = max(ny1,0)
    nz1 = max(nz1,0)
    nx2 = min(nx-1,nx2)
    ny2 = min(ny-1,ny2)
    nz2 = min(nz-1,nz2)
  endif


  write(*,'("Cropping field from " &
  &,"0:",i4," | 0:",i4," | 0:",i4,&
  &"   to subset   "&
  &,i4,":",i2,":",i4," | ",i4,":",i2,":",i4," | ",i4,":",i2,":",i4)')&
  nx-1,ny-1,nz-1,nx1,nxs,nx2,ny1,nys,ny2,nz1,nzs,nz2



  !-----------------------------------------------------------------------------
  ! compute dimensions of reduced subset:
  nx_red = nx1 + floor( dble(nx2-nx1)/dble(nxs) )  - nx1 + 1
  ny_red = ny1 + floor( dble(ny2-ny1)/dble(nys) )  - ny1 + 1
  nz_red = nz1 + floor( dble(nz2-nz1)/dble(nzs) )  - nz1 + 1
  write (*,'("Size of subset is ",3(i4,1x))') nx_red, ny_red, nz_red
  !-----------------------------------------------------------------------------

  ! we figured out how big the subset array is
  allocate ( field(0:nx_red-1,0:ny_red-1,0:nz_red-1) )


  call Fetch_attributes( fname_in, dsetname_in, nx_file,ny_file,nz_file,&
  xl_file,yl_file,zl_file,time )

  !-----------------------------------------------------------------------------
  ! load the file
  ! the basic idea is to just allocate the smaller field, and use hdf5
  ! to read just this field from the input file.
  !-----------------------------------------------------------------------------
  ! Initialize HDF5 library and Fortran interfaces.
  call h5open_f(error)

  ! Setup file access property list with parallel I/O access.  this
  ! sets up a property list ("plist_id") with standard values for
  ! FILE_ACCESS
  call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
  ! this modifies the property list and stores MPI IO
  ! comminucator information in the file access property list
  call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  ! open the file in parallel
  call h5fopen_f (fname_in, H5F_ACC_RDWR_F, file_id, error, plist_id)
  ! this closes the property list (we'll re-use it)
  call h5pclose_f(plist_id, error)

  ! Definition of memory distribution
  dimensions_file  = (/ nx_file, ny_file, nz_file/)
  dimensions_local = (/ nx_red , ny_red,  nz_red /)
  offset = (/ nx1, ny1, nz1 /)
  stride = (/ nxs, nys, nzs /)
  chunking_dims = 32 !min(nx_red,ny_red,nz_red)

  !----------------------------------------------------------------------------
  ! Read actual field from file (dataset)
  !----------------------------------------------------------------------------
  ! dataspace in the file: contains all data from all procs
  call h5screate_simple_f(rank, dimensions_file, filespace, error)
  ! dataspace in memory: contains only local data
  call h5screate_simple_f(rank, dimensions_local, memspace, error)

  ! Create chunked dataset
  call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
  call h5pset_chunk_f(plist_id, rank, chunking_dims, error)

  ! Open an existing dataset.
  call h5dopen_f(file_id, dsetname_in, dset_id, error)

  ! Select hyperslab in the file.
  call h5dget_space_f(dset_id, filespace, error)
  call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, count, &
  error , stride, dimensions_local)

  ! Create property list for collective dataset read
  call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
  call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
  ! actual read is the next command:
  call h5dread_f( dset_id, H5T_NATIVE_DOUBLE, field, dimensions_local, error, &
  mem_space_id = memspace, file_space_id = filespace, xfer_prp = plist_id )

  ! check if we loaded crap
  call checknan(field,"recently loaded field")

  call h5sclose_f(filespace, error)
  call h5sclose_f(memspace, error)
  call h5pclose_f(plist_id, error)
  call h5dclose_f(dset_id, error)
  call h5fclose_f(file_id,error)
  call H5close_f(error)

  ! set up dimensions in global variables, since save_field_hdf5 relies on this
  ra = 0
  rb(1) = nx_red-1
  rb(2) = ny_red-1
  rb(3) = nz_red-1
  nx = nx_red
  ny = ny_red
  nz = nz_red
  dx = xl_file / dble(nx_file)
  dy = yl_file / dble(ny_file)
  dz = zl_file / dble(nz_file)
  xl = dx + dble(nx1+(nx_red-1)*nxs)*dx - dble(nx1)*dx
  yl = dy + dble(ny1+(ny_red-1)*nys)*dy - dble(ny1)*dy
  zl = dz + dble(nz1+(nz_red-1)*nzs)*dz - dble(nz1)*dz

  ! Done! Write extracted subset to disk and be happy with the result
  call save_field_hdf5 ( time, fname_out(1:index(fname_out,'.h5')-1), field, dsetname_out )

end subroutine extract_subset


!-------------------------------------------------------------------------------
! copy one hdf5 file to another one, with different name. Not you cannot do this
! in terminal with "cp" since "cp" does not touch the dataset name in the file,
! which is then not conformal to flusi naming convention.
! call:
! ./flusi --postprocess --cp ux_00000.h5 new_00000.h5
!-------------------------------------------------------------------------------
subroutine copy_hdf_file()
  use mpi
  use vars
  implicit none
  character(len=strlen) :: fname_in, fname_out, dsetname_in, dsetname_out
  character(len=strlen) :: xset,yset,zset
  real(kind=pr)::time
  integer :: ix,iy,iz,i
  ! reduced domain size
  integer :: nx1,nx2, ny1,ny2, nz1,nz2, nxs,nys,nzs
  ! sizes of the new array
  integer :: nx_red, ny_red, nz_red, ix_red, iy_red, iz_red
  ! reduced domain extends
  real(kind=pr) :: xl1, yl1, zl1
  ! input field
  real(kind=pr), dimension(:,:,:), allocatable :: field_in


  if (mpisize/=1) then
    write(*,*) "./flusi --postprocess --cp is a SERIAL routine, use 1CPU only"
    call abort()
  endif



  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )

  ! get filename to save file to
  call get_command_argument(4,fname_out)

  dsetname_in = fname_in ( 1:index( fname_in, '_' )-1 )
  dsetname_out = fname_out ( 1:index( fname_out, '_' )-1 )

  call fetch_attributes( fname_in, dsetname_in, nx, ny, nz, xl, yl, zl, time )
  ra=0
  rb=(/nx-1,ny-1,nz-1/)
  write(*,*) "copying ",trim(adjustl(fname_in)), " to ", trim(adjustl(fname_out))

  allocate ( field_in(0:nx-1,0:ny-1,0:nz-1) )
  call read_single_file_serial(fname_in,field_in)

  call save_field_hdf5 ( time, fname_out(1:index(fname_out,'.h5')-1), field_in, dsetname_out )

  deallocate (field_in)
end subroutine copy_hdf_file




!-------------------------------------------------------------------------------
! Upsampling from a source resolution to a target resolution
! ./flusi -p --upsample source.h5 target.h5 256 256 526
!-------------------------------------------------------------------------------
! We first read in the original field from the source file, with it's resolution
! and domain size and timestamp.
!-------------------------------------------------------------------------------
subroutine upsample()
  use vars
  use p3dfft_wrapper
  implicit none
  character(len=strlen) :: fname_in, fname_out, dsetname_in, dsetname_out, tmp
  integer :: nx_new, ny_new, nz_new
  integer :: nx_org, ny_org, nz_org, ix_org,iy_org,iz_org, ix_new,iy_new,iz_new
  integer :: i,j,k
  real(kind=pr) :: time, kx_org,ky_org,kz_org, kx_new,ky_new,kz_new
  complex(kind=pr),dimension(:,:,:),allocatable :: uk_org, uk_new
  real(kind=pr),dimension(:,:,:),allocatable :: u_org, u_new
  integer, dimension(1:3) :: ra_org,rb_org,ca_org,cb_org,ra_new,rb_new,ca_new,cb_new
  if (mpisize/=1) then
    write(*,*) "./flusi --postprocess --upsample is a SERIAL routine, use 1CPU only"
    call abort()
  endif

  ! get file to read pressure from and check if this is present
  call get_command_argument(3,fname_in)
  call check_file_exists( fname_in )
  dsetname_in = fname_in ( 1:index( fname_in, '_' )-1 )

  call get_command_argument(4,fname_out)
  dsetname_out = fname_out ( 1:index( fname_out, '_' )-1 )

  ! read target resolution from command line
  call get_command_argument(5,tmp)
  read (tmp,*) nx_new
  call get_command_argument(6,tmp)
  read (tmp,*) ny_new
  call get_command_argument(7,tmp)
  read (tmp,*) nz_new

  write(*,'("Target resolution= ",3(i4,1x))') nx_new, ny_new, nz_new

  call fetch_attributes( fname_in, dsetname_in, nx_org, ny_org, nz_org, xl, yl, zl, time )
  write(*,'("Origin resolution= ",3(i4,1x))') nx_org,ny_org,nz_org

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl


  !-----------------------
  write(*,*) "Initializing small FFT and transforming source field to k-space"
  nx = nx_org
  ny = ny_org
  nz = nz_org
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  call fft_initialize
  ra_org = ra
  rb_org = rb
  ca_org = ca
  cb_org = cb
  !-----------------------

  allocate(u_org(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk_org(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))

  call fft_unit_test( u_org, uk_org )

  write(*,*) "Reading file "//trim(adjustl(fname_in))
  call read_single_file_serial(fname_in,u_org)

  call fft(inx=u_org, outk=uk_org)

  deallocate(u_org)

  call fft_free
  !-----------------------
  write(*,*) "Initializing big FFT and copying source Fourier coefficients to &
             & target field in k-space"

  nx=nx_new
  ny=ny_new
  nz=nz_new
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)
  call fft_initialize
  ra_new = ra
  rb_new = rb
  ca_new = ca
  cb_new = cb
  !-----------------------

  allocate(u_new(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  allocate(uk_new(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)))
  call fft_unit_test( u_new, uk_new )
  uk_new = dcmplx(0.d0,0.d0)

  !------------------------------------------------------------------
  do iz_org=ca_org(1),cb_org(1)
     nx=nx_org; ny=ny_org; nz=nz_org
     kz_org=wave_z(iz_org)
     ! find corresponding iz_new
     do iz_new=ca_new(1),cb_new(1)
       nx=nx_new; ny=ny_new; nz=nz_new
       kz_new = wave_z(iz_new)
       if (kz_new==kz_org) exit
     enddo
     ! we now have the pair (iz_org, iz_new)
     !------------------------------------------------------------------
     do iy_org=ca_org(2),cb_org(2)
          nx=nx_org; ny=ny_org; nz=nz_org
          ky_org=wave_y(iy_org)
          ! find corresponding iz_new
          do iy_new=ca_new(2),cb_new(2)
            nx=nx_new; ny=ny_new; nz=nz_new
            ky_new = wave_y(iy_new)
            if (ky_new==ky_org) exit
          enddo
          ! we now have the pair (iy_org, iy_new)
          !------------------------------------------------------------------
          do ix_org=ca_org(3),cb_org(3)
              nx=nx_org; ny=ny_org; nz=nz_org
              kx_org=wave_x(ix_org)
              ! find corresponding iz_new
              do ix_new=ca_new(3),cb_new(3)
                nx=nx_new; ny=ny_new; nz=nz_new
                kx_new = wave_x(ix_new)
                if (kx_new==kx_org) exit
              enddo
              ! we now have the pair (ix_org, ix_new)

              ! copy the old Fourier coefficients to the new field
              uk_new(iz_new,iy_new,ix_new) = uk_org(iz_org,iy_org,ix_org)
          enddo
     enddo
  enddo

  deallocate( uk_org )

  ! transform the zero-padded Fourier coefficients back to physical space. this
  ! is the upsampled (=interpolated) field.
  write(*,*) "transforming zero-padded Fourier coefficients back to x-space"
  call ifft(ink=uk_new,outx=u_new)

  deallocate( uk_new )

  ! save the final result to the specified file
  write(*,*) "Saving upsampled field to " // trim(adjustl(fname_out))
  call save_field_hdf5(time,fname_out,u_new,dsetname_out)

  deallocate( u_new )
end subroutine upsample


!-------------------------------------------------------------------------------
! ./flusi --postprocess --spectrum ux_00000.h5 uy_00000.h5 uz_00000.h5 spectrum.dat
!-------------------------------------------------------------------------------
! NOTE: I actually did not figure out what happens if xl=yl=zl/=2*pi
! which is a rare case in all isotropic turbulence situtations, and neither
! of the corresponding routines have been tested for that case.
subroutine post_spectrum()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, spectrum_file
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time, sum_u
  real(kind=pr), dimension(:), allocatable :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin
  integer :: mpicode, k

  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,spectrum_file)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
    write (*,*) "Error in arguments, files do not start with ux uy and uz"
    write (*,*) "note files have to be in the right order"
    call abort()
  endif

  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  if (mpirank==0) then
    write(*,'("Computing spectrum of ",A,1x,A,1x,A)') &
    trim(adjustl(fname_ux)), trim(adjustl(fname_uy)), trim(adjustl(fname_uz))
    write(*,'("Spectrum will be written to ",A)') trim(adjustl(spectrum_file))
    write(*,'("Resolution is ",i4,1x,i4,1x,i4)') nx, ny, nz
    write(*,'("Domain size is", es12.4,1x,es12.4,1x,es12.4)') xl, yl ,zl
  endif

  call fft_initialize() ! also initializes the domain decomp


  call MPI_barrier (MPI_COMM_world, mpicode)
  write (*,'("mpirank=",i5," x-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,&
  &") k-space=(",i4,":",i4," |",i4,":",i4," |",i4,":",i4,")")') &
  mpirank, ra(1),rb(1), ra(2),rb(2),ra(3),rb(3), ca(1),cb(1), ca(2),cb(2),ca(3),cb(3)
  call MPI_barrier (MPI_COMM_world, mpicode)


  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(S_Ekinx(0:nx-1),S_Ekiny(0:nx-1),S_Ekinz(0:nx-1),S_Ekin(0:nx-1))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft (uk(:,:,:,1),u(:,:,:,1))
  call fft (uk(:,:,:,2),u(:,:,:,2))
  call fft (uk(:,:,:,3),u(:,:,:,3))

  ! compute the actual spectrum
  call compute_spectrum( time,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

  ! on root, write it to disk
  if (mpirank == 0) then
    open(10,file=spectrum_file,status='replace')
    write(10,'(5(A15,1x))') '%   K ','E_u(K)','E_ux(K)','E_uy(K)','E_uz(K)'
    do k=0,nx-1
      write(10,'(5(1x,es15.8))') dble(k),S_Ekin(k),S_Ekinx(k),S_Ekiny(k),S_Ekinz(k)
    enddo

    sum_u=0.0d0
    do k=1,nx-1
      sum_u=sum_u +S_Ekin(k)
    enddo
    write(10,*) '% Etot = ',sum_u
    write(10,*) '% time = ',time
    close(10)
  endif

  deallocate (u)
  deallocate (uk)
  deallocate(S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
  call fft_free()

end subroutine post_spectrum


!-------------------------------------------------------------------------------
! ./flusi --postprocess --turbulence-analysis ux_00000.h5 uy_00000.h5 uz_00000.h5 nu outfile.dat
! NOTE: I actually did not figure out what happens if xl=yl=zl/=2*pi
! which is a rare case in all isotropic turbulence situtations, and neither
! of the corresponding routines have been tested for that case.
!-------------------------------------------------------------------------------
subroutine turbulence_analysis()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, viscosity, outfile
  complex(kind=pr),dimension(:,:,:,:),allocatable :: uk,vork
  real(kind=pr),dimension(:,:,:,:),allocatable :: u, vor
  real(kind=pr) :: time, epsilon_loc, epsilon, fact, E, u_rms,lambda_macro,lambda_micro
  real(kind=pr), dimension(:), allocatable :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin
  integer :: ix,iy,iz, mpicode
  real(kind=pr)::kx,ky,kz,kreal


  call get_command_argument(3,fname_ux)
  call get_command_argument(4,fname_uy)
  call get_command_argument(5,fname_uz)
  call get_command_argument(6,viscosity)
  call get_command_argument(7,outfile)
  read(viscosity,*) nu

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )

  if ((fname_ux(1:2).ne."ux").or.(fname_uy(1:2).ne."uy").or.(fname_uz(1:2).ne."uz")) then
    write (*,*) "Error in arguments, files do not start with ux uy and uz"
    write (*,*) "note files have to be in the right order"
    call abort()
  endif

  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    write(17,'(A)') "-----------------------------------"
    write(17,'(A)') "FLUSI turbulence analysis"
    write(17,'("call: ./flusi -p --turbulence-analysis ",5(A,1x))') trim(adjustl(fname_ux)),&
    trim(adjustl(fname_uy)),trim(adjustl(fname_uz)),trim(adjustl(viscosity)),trim(adjustl(outfile))
    write(17,'(A)') "-----------------------------------"
  endif



  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(vor(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(vork(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3))
  allocate(S_Ekinx(0:nx-1),S_Ekiny(0:nx-1),S_Ekinz(0:nx-1),S_Ekin(0:nx-1))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )

  call fft3 (inx=u,outk=uk)
  call curl (uk,vork)
  call ifft3 (ink=vork,outx=vor)

  ! compute spectrum
  call compute_spectrum( time,uk,S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

  !-----------------------------------------------------------------------------
  ! dissipation rate from velocity in Fourier space
  !-----------------------------------------------------------------------------
  do iz=ca(1),cb(1)
    kz=wave_z(iz)
    do iy=ca(2),cb(2)
      ky=wave_y(iy)
      do ix=ca(3),cb(3)
        kx=wave_x(ix)
        kreal = ( (kx*kx)+(ky*ky)+(kz*kz) )

        if ( ix==0 .or. ix==nx/2 ) then
          E=dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2)/2. &
          +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2)/2. &
          +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)/2.
        else
          E=dble(real(uk(iz,iy,ix,1))**2+aimag(uk(iz,iy,ix,1))**2) &
          +dble(real(uk(iz,iy,ix,2))**2+aimag(uk(iz,iy,ix,2))**2) &
          +dble(real(uk(iz,iy,ix,3))**2+aimag(uk(iz,iy,ix,3))**2)
        endif

        epsilon_loc = epsilon_loc + kreal * E
      enddo
    enddo
  enddo

  epsilon_loc = 2.d0 * nu * epsilon_loc

  call MPI_ALLREDUCE(epsilon_loc,epsilon,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)

  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') epsilon, "Dissipation rate from velocity in Fourier space"
  endif



  !-----------------------------------------------------------------------------
  ! dissipation rate from vorticty
  !-----------------------------------------------------------------------------
  epsilon_loc = nu * sum(vor(:,:,:,1)**2+vor(:,:,:,2)**2+vor(:,:,:,3)**2)
  call MPI_REDUCE(epsilon_loc,epsilon,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)

  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') epsilon/dble(nx*ny*nz), "Dissipation rate from vorticity"
  endif

  !-----------------------------------------------------------------------------
  ! dissipation rate from spectrum; see Ishihara, Kaneda "High
  ! resolution DNS of incompressible Homogeneous forced turbulence -time dependence
  ! of the statistics" or my thesis
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    epsilon=0.0
    do ix = 0,nx-1
      epsilon = epsilon + 2.d0 * nu * dble(ix**2) * S_Ekin(ix)
    enddo
    write(17,'(g15.8,5x,A)') epsilon, "Dissipation rate from spectrum"
  endif

  !-----------------------------------------------------------------------------
  ! energy from velocity
  !-----------------------------------------------------------------------------
  epsilon_loc = 0.5d0*sum(u(:,:,:,1)**2+u(:,:,:,2)**2+u(:,:,:,3)**2)
  call MPI_REDUCE(epsilon_loc,E,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
  MPI_COMM_WORLD,mpicode)

  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') E/dble(nx*ny*nz), "energy from velocity"
  endif

  !-----------------------------------------------------------------------------
  ! energy from spectrum
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    E=0.0
    do ix = 0,nx-1
      E = E + S_Ekin(ix)
    enddo
    write(17,'(g15.8,5x,A)') E, "energy from spectrum"
  endif

  u_rms=dsqrt(2.d0*E/3.d0)
  !-----------------------------------------------------------------------------
  ! kolmogrov scales
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    write(17,'(g15.8,5x,A)') (nu**3 / epsilon)**(0.25d0), "kolmogorov length scale"
    write(17,'(g15.8,5x,A)') (nu / epsilon)**(0.5d0), "kolmogorov time scale"
    write(17,'(g15.8,5x,A)') (nu*epsilon)**(0.25d0), "kolmogorov velocity scale"
    write(17,'(g15.8,5x,A)') u_rms, "RMS velocity"
  endif

  !-----------------------------------------------------------------------------
  ! taylor scales
  !-----------------------------------------------------------------------------
  if (mpirank==0) then
    lambda_micro = (15.d0*nu*u_rms**2 / epsilon)**(0.5d0)
    lambda_macro=0.0
    do ix = 1,nx-1
      lambda_macro = lambda_macro + pi/(2.d0*u_rms**2) * S_Ekin(ix) / dble(ix)
    enddo
    write(17,'(g15.8,5x,A)') lambda_micro, "taylor micro scale"
    write(17,'(g15.8,5x,A)') lambda_macro, "taylor macro scale"
    write(17,'(g15.8,5x,A)') u_rms*lambda_macro/nu, "Renolds taylor macro scale"
    write(17,'(g15.8,5x,A)') u_rms*lambda_micro/nu, "Renolds taylor micro scale"
    write(17,'(g15.8,5x,A)') lambda_macro/u_rms, "eddy turnover time"
    write(17,'(g15.8,5x,A)') (2./3.)*(dble(nx/2-1)), "kmax"
    write(17,'(g15.8,5x,A)') (2./3.)*(dble(nx/2-1))*(nu**3 / epsilon)**(0.25d0), "kmax*eta"
  endif

  if(mpirank==0) close(17)

  deallocate (u,vor,vork)
  deallocate (uk)
  deallocate (S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin)
  call fft_free()

end subroutine turbulence_analysis


!-------------------------------------------------------------------------------
! ./flusi -p --TKE-mean ekinavg_00.h5 uavgx_00.h5 uavgy_00.h5 uavgz_00.h5 tkeavg_000.h5
! From the time-avg kinetic energy field and the components of the time avg
! velocity field, compute the mean turbulent kinetic energy.
! See TKE note 18 feb 2015 (Thomas) and 13 feb 2015 (Dmitry)
! Can be done in parallel
!-------------------------------------------------------------------------------
subroutine TKE_mean()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: fname_ux, fname_uy, fname_uz, dsetname, fname_ekin, outfile
  real(kind=pr),dimension(:,:,: ),allocatable :: ekin
  real(kind=pr),dimension(:,:,:,:),allocatable :: u
  real(kind=pr) :: time

  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,fname_ux)
  call get_command_argument(5,fname_uy)
  call get_command_argument(6,fname_uz)
  call get_command_argument(7,outfile)

  call check_file_exists( fname_ux )
  call check_file_exists( fname_uy )
  call check_file_exists( fname_uz )
  call check_file_exists( fname_ekin )

  if (mpirank == 0) then
    write(*,*) "Processing "//trim(adjustl(fname_ux))//" "//trim(adjustl(fname_uy))//&
    &" "//trim(adjustl(fname_uz))//" and "//fname_ekin
  endif

  dsetname = fname_ux ( 1:index( fname_ux, '_' )-1 )
  call fetch_attributes( fname_ux, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate(u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3))
  allocate(ekin(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  call read_single_file ( fname_ux, u(:,:,:,1) )
  call read_single_file ( fname_uy, u(:,:,:,2) )
  call read_single_file ( fname_uz, u(:,:,:,3) )
  call read_single_file ( fname_ekin, ekin )

  ekin = ekin - 0.5d0*(u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2)

  dsetname = outfile ( 1:index( outfile, '_' )-1 )
  if (mpirank==0) write(*,*) "Wrote to "//trim(adjustl(outfile))//" "//trim(adjustl(dsetname))
  call save_field_hdf5 ( time,outfile,ekin,dsetname)


  deallocate (u,ekin)
  call fft_free()

end subroutine tke_mean


!-------------------------------------------------------------------------------
! ./flusi -p --max-over-x tkeavg_000.h5 outfile.dat
! This function reads in the specified *.h5 file and outputs the maximum value
! max_yz(x) into the specified ascii-outfile
! It may be used rarely, but we needed it for turbulent bumblebees.
! Can be done in parallel.
!-------------------------------------------------------------------------------
subroutine max_over_x()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: dsetname, fname_ekin, outfile
  real(kind=pr),dimension(:,:,:),allocatable :: u
  real(kind=pr), dimension(:), allocatable :: umaxx,umaxx_loc
  real(kind=pr) :: time
  integer :: ix, mpicode


  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,outfile)
  call check_file_exists( fname_ekin )

  dsetname = fname_ekin ( 1:index( fname_ekin, '_' )-1 )
  call fetch_attributes( fname_ekin, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( umaxx(0:nx-1) )
  allocate( umaxx_loc(0:nx-1) )

  call read_single_file ( fname_ekin, u )

  do ix=0,nx-1
    umaxx_loc(ix)=maxval(u(ix,:,:))
  enddo

  call MPI_ALLREDUCE(umaxx_loc,umaxx,nx,MPI_DOUBLE_PRECISION,MPI_MAX,&
       MPI_COMM_WORLD,mpicode)

  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    write(17,'(A)') "%-----------------------------------"
    write(17,'(A)') "%FLUSI max-over-x file="//trim(adjustl(fname_ekin))
    write(17,'(A)') "%-----------------------------------"
    do ix=0,nx-1
      write(17,'(es15.8)') umaxx(ix)
    enddo
    close(17)
  endif



  deallocate (u,umaxx,umaxx_loc)
  call fft_free()

end subroutine max_over_x


!-------------------------------------------------------------------------------
! ./flusi -p --max-over-x tkeavg_000.h5 outfile.dat
! This function reads in the specified *.h5 file and outputs the maximum value
! max_yz(x) into the specified ascii-outfile
! It may be used rarely, but we needed it for turbulent bumblebees.
! Can be done in parallel.
!-------------------------------------------------------------------------------
subroutine mean_over_x_subdomain()
  use vars
  use p3dfft_wrapper
  use basic_operators
  use mpi
  implicit none
  character(len=strlen) :: dsetname, fname_ekin, outfile
  real(kind=pr),dimension(:,:,:),allocatable :: u, mask
  real(kind=pr), dimension(:), allocatable :: umaxx,umaxx_loc
  real(kind=pr) :: time,x,y,z
  integer :: ix,iy,iz, mpicode


  call get_command_argument(3,fname_ekin)
  call get_command_argument(4,outfile)
  call check_file_exists( fname_ekin )

  dsetname = fname_ekin ( 1:index( fname_ekin, '_' )-1 )
  call fetch_attributes( fname_ekin, dsetname, nx, ny, nz, xl, yl, zl, time )

  pi=4.d0 *datan(1.d0)
  scalex=2.d0*pi/xl
  scaley=2.d0*pi/yl
  scalez=2.d0*pi/zl
  dx = xl/dble(nx)
  dy = yl/dble(ny)
  dz = zl/dble(nz)

  call fft_initialize() ! also initializes the domain decomp

  allocate( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  allocate( umaxx(0:nx-1) )
  allocate( umaxx_loc(0:nx-1) )

  call read_single_file ( fname_ekin, u )

  ! set a 1/0 mask to cancel values out of the bounds we want
  mask =0.d0
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x=dble(ix)*dx-0.5*xl
        y=dble(iy)*dy-0.5*yl
        z=dble(iz)*dz-0.5*zl

        if ( abs(y)<=1.3d0 .and. abs(z)<=1.3d0 ) then
          mask(ix, iy, iz) = 1.d0
        endif
      enddo
    enddo
  enddo


  u = u*mask

  do ix=0,nx-1
    umaxx_loc(ix)=sum(u(ix,:,:)) / sum(mask(ix,:,:))
  enddo

  call MPI_ALLREDUCE(umaxx_loc,umaxx,nx,MPI_DOUBLE_PRECISION,MPI_SUM,&
  MPI_COMM_WORLD,mpicode)



  if(mpirank==0) then
    write(*,*) " OUTPUT will be written to "//trim(adjustl(outfile))
    open(17,file=trim(adjustl(outfile)),status='replace')
    write(17,'(A)') "%-----------------------------------"
    write(17,'(A)') "%FLUSI mean over x (boxed) file="//trim(adjustl(fname_ekin))
    write(17,'(A)') "%-----------------------------------"
    do ix=0,nx-1
      write(17,'(es15.8)') umaxx(ix)
    enddo
    close(17)
  endif



  deallocate (u,umaxx,umaxx_loc,mask)
  call fft_free()

end subroutine mean_over_x_subdomain
