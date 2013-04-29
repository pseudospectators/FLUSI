!=========================================================
!      Collect the data from files
!      produced with FLUSI parallel code
!      and save them in a binary file
!=========================================================

program convert_mpiio
  use mpi
  implicit none

  integer, parameter :: pr_in = 8, pr_out = 4 ! precision for input and output (INPUT PRECISION MUST BE KNOWN!!!)
  integer, parameter :: mpireal = MPI_DOUBLE_PRECISION  ! double precision array for input
  integer, dimension(MPI_STATUS_SIZE) :: mpistatus
  integer :: nx, ny, nz
  integer :: mpicode, filedesc
  integer :: ix, iy, iz, iSaveAscii
  real (kind=pr_in ), dimension (:,:,:), allocatable :: field_in
  real (kind=pr_out), dimension (:,:,:), allocatable :: field_out
  character (len=128) :: fname, nx_str, ny_str, nz_str, ascii_str
  logical :: file_exists
  
  
  call MPI_INIT (mpicode)
  
  
  !--------------------------------------------------
  ! read in command line arguments
  !--------------------------------------------------  
  call get_command_argument(1, fname)
  if (len_trim(fname)==0) then
    write (*,*) "You forgot to tell me which file to process..."
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
  
  
  
  !--------------------------------------------------
  ! check if desired file exists
  !--------------------------------------------------
  inquire(file=trim(fname)//'.mpiio',exist=file_exists)   
  if (file_exists .eqv. .false.) then
    write (*,'("ERROR!! File ",A,".mpiio NOT found")') trim(fname)
    stop
  endif
  
  allocate ( field_in(0:(nx-1),0:(ny-1),0:(nz-1)), field_out(0:(nx-1),0:(ny-1),0:(nz-1)) )
 
  !--------------------------------------------------
  ! read MPI data
  !-------------------------------------------------- 
  
  call MPI_FILE_OPEN (MPI_COMM_WORLD,trim(fname)//'.mpiio',MPI_MODE_RDONLY,MPI_INFO_NULL,filedesc,mpicode)
  call MPI_FILE_READ_ORDERED (filedesc,field_in,nx*ny*nz,mpireal,mpistatus,mpicode)
  call MPI_FILE_CLOSE (filedesc,mpicode)

  write(*,'("Converting ",A,".mpiio  -> ",A,".binary Min:Max=",es12.4,":",es12.4,1x,"nx:ny:nz=",i3,":",i3,":",i3)') &
  trim(fname), trim(fname), minval (field_in), maxval (field_in),nx,ny,nz
  
  
  field_out = real(field_in, kind=pr_out)  ! convert field to output precision
  ! ------------------------------------
  ! Write *,binary file (for vapor et al)
  ! ------------------------------------
  open (12, file = trim(fname)//".binary", form='unformatted', status='replace')
  write (12) (((field_out (ix,iy,iz), ix=0, nx-1), iy=0, ny-1), iz=0, nz-1)
  close (12)


  deallocate ( field_in, field_out )

  
  call MPI_FINALIZE (mpicode)

end program convert_mpiio


subroutine help
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "		converter "
  write (*,*) " *.mpiio -> *.binary"
  write (*,*) "--------------------------------------------------------------"
  write (*,*) "  usage: ./convert_mpiio filename nx ny nz "
  write (*,*) ""
  write (*,*) " [filename]: the base file name WITHOUT the .mpiio suffix"
  write (*,*) " [nx, ny, nz] is the resolution"
  write (*,*) "--------------------------------------------------------------"

end subroutine

