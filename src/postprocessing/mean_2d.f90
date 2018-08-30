!-------------------------------------------------------------------------------
! ./flusi -p --mean-2D [x,y,z] infile_000.h5 outfile.dat
! This function reads in the specified *.h5 file and outputs the average over two
! directions as a function of the remaining direction.
! e.g., ./flusi -p --mean-2D z infile_000.h5 outfile.dat
! averages over the x and y direction
! e.g., ./flusi -p --mean-2D all infile_000.h5 outfile.dat
! will loop over x,y,z and output all three to different files
!-------------------------------------------------------------------------------
subroutine mean_2d(help)
  use vars
  use p3dfft_wrapper
  use basic_operators
  use module_helpers
  use mpi
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: infile, outfile, direction
  real(kind=pr),dimension(:,:,:),allocatable :: u
  real(kind=pr) :: time
  integer :: ix,iy,iz

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --mean-2D [x,y,z,all] infile_000.h5 outfile.dat"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "! This function reads in the specified *.h5 file and outputs the average over two"
    write(*,*) "! directions as a function of the remaining direction."
    write(*,*) "! e.g., ./flusi -p --mean-2D z infile_000.h5 outfile.dat"
    write(*,*) "! averages over the x and y direction"
    write(*,*) "! e.g., ./flusi -p --mean-2D all infile_000.h5 outfile.dat"
    write(*,*) "! will loop over x,y,z and output all three to different files"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: nope"
    return
  endif


  if (mpisize/=1) then
    ! the reason for this is simplicity. if the y-z direction is nonlocal in memory
    ! the avg is more complicated. only the x direction is always contiguous.
    ! Plus: this acts on one field only, and it usually fits in the memory.
    call abort(4343,"./flusi --postprocess --mean-2D is a SERIAL routine, use 1CPU only")
  endif

  call get_command_argument(3,direction)
  call get_command_argument(4,infile)
  call get_command_argument(5,outfile)
  call check_file_exists( infile )

  write(*,*) "computing average in a 2D plane"
  write(*,*) "infile="//trim(adjustl(infile))
  write(*,*) "outfile="//trim(adjustl(outfile))

  call fetch_attributes( infile, nx, ny, nz, xl, yl, zl, time, nu, origin )
  ! initialize code and scaling factors for derivatives, also domain decomposition
  call fft_initialize()

  ! allocate memory and read file
  allocate( u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
  call read_single_file ( infile, u )


  !-----------------------------------------------------------------------------
  ! the rest of the code depends on the direction
  !-----------------------------------------------------------------------------
  select case (direction)
  case ("x")
      open(17,file=trim(adjustl(outfile)),status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      write(*,*) "using X direction, thus averaging over y-z"
      ! compute actual mean value in the y-z plane, directly write to file
      do ix=0,nx-1
        write(17,'(es15.8)') sum( u(ix,:,:) ) / dble(ny*nz)
      enddo
      close(17)
  case ("y")
      open(17,file=trim(adjustl(outfile)),status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      write(*,*) "using Y direction, thus averaging over x-z"
      ! compute actual mean value in the x-z plane, directly write to file
      do iy=0,ny-1
        write(17,'(es15.8)') sum( u(:,iy,:) ) / dble(nx*nz)
      enddo
      close(17)
  case ("z")
      open(17,file=trim(adjustl(outfile)),status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      write(*,*) "using Z direction, thus averaging over x-y"
      ! compute actual mean value in the y-z plane, directly write to file
      do iz=0,nz-1
        write(17,'(es15.8)') sum( u(:,:,iz) ) / dble(ny*nx)
      enddo
      close(17)
  case ("all")
      write(*,*) "we compute all three possible averages"
      write(*,*) "--------------------------------------"
      write(*,*) "using X direction, thus averaging over y-z"
      ! compute actual mean value in the y-z planes, directly write to file
      open(17,file=trim(adjustl(outfile))//"_x",status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      do ix=0,nx-1
        write(17,'(es15.8)') sum( u(ix,:,:) ) / dble(ny*nz)
      enddo
      close(17)

      write(*,*) "using Y direction, thus averaging over x-z"
      ! compute actual mean value in the x-z plane, directly write to file
      open(17,file=trim(adjustl(outfile))//"_y",status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      do iy=0,ny-1
        write(17,'(es15.8)') sum( u(:,iy,:) ) / dble(nx*nz)
      enddo
      close(17)

      write(*,*) "using Z direction, thus averaging over x-y"
      ! compute actual mean value in the y-z plane, directly write to file
      open(17,file=trim(adjustl(outfile))//"_z",status='replace')
      call postprocessing_ascii_header(17)
      write(17,'(A)') "%-----------------------------------"
      write(17,'(A)') "% FLUSI --mean-2D file="//trim(adjustl(infile))
      write(17,'(A)') "% direction="//trim(adjustl(direction))
      write(17,'(A)') "%-----------------------------------"
      do iz=0,nz-1
        write(17,'(es15.8)') sum( u(:,:,iz) ) / dble(ny*nx)
      enddo
      close(17)
  case default
      call abort(131013,"Bad choice for direction "//trim(adjustl(direction)) )
  end select

  deallocate (u)
  call fft_free()

end subroutine mean_2d
