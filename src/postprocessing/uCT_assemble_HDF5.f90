!-------------------------------------------------------------------------------
! ./flusi --postprocess --time-avg file_list.txt avgx_0000.h5
! Reads in a list of files from a file, then loads one file after the other and
! computes the average field, which is then stored in the specified file.
!-------------------------------------------------------------------------------
subroutine uCT_assemble_HDF5(help)
  use vars
  use p3dfft_wrapper
  use module_helpers
  use ghosts
  implicit none
  logical, intent(in) :: help
  character(len=strlen) :: fname, fname_bin, fname_avg, fname_color
  real(kind=pr), dimension(:,:,:), allocatable :: field_avg, field, field_color
  integer :: ix, iy ,iz, io_error=0, i=0
  real(kind=pr) :: time

  if (help.and.root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "./flusi -p --uCT-assemble file_list.txt avgx_0000.h5 color_0000.h5"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) " Reads in a list of files from a file, computes their sum and saves the result to disk"
    write(*,*) " We also compute a color file with constant values for each non-zero region of a file."
    write(*,*) " This makes only sense for uCT data, where one typically has a bunch of files e.g. for legs"
    write(*,*) " head etc and one likes to assemble the general picture. Only one file per pixel is nonzero."
    write(*,*) " "
    write(*,*) " "
    write(*,*) " "
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Parallel: yes"
    return
  endif

  ! 1 ghost node
  ng = 1


  call get_command_argument(3,fname)
  call get_command_argument(4,fname_avg)
  call get_command_argument(5,fname_color)

  !-----------------------------------------------------------------------------
  ! check if input file exists, the file contains the list of h5 files to be avg
  !-----------------------------------------------------------------------------
  call check_file_exists ( fname )
  if (root) write(*,*) "Reading list of files from "//fname

  !-----------------------------------------------------------------------------
  ! read in the file, loop over lines
  !-----------------------------------------------------------------------------
  open( unit=14,file=fname, action='read', status='old')
  do while (io_error==0)
    ! fetch current filename
    read (14,'(A)', iostat=io_error) fname_bin

    if (io_error == 0) then
      if (root) write(*,*) "read "//trim(adjustl(fname_bin))

      call check_file_exists ( fname_bin )
      call fetch_attributes( fname_bin, nx, ny, nz, xl, yl, zl, time, nu, origin )

      !-------------------------
      ! initialization
      !-------------------------
      ! first time? allocate then.
      if ( .not. allocated(field_avg) ) then
        ! initialize code and domain decomposition, but do not use FFTs
        call decomposition_initialize()
        call setup_cart_groups()
        ! allocate memory
        allocate(field_avg(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
        allocate(field_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
        allocate(field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))

        field_avg = 0.d0
        field_color = 0.d0
      endif

      ! read the field from file
      call read_single_file( fname_bin, field(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )
      call synchronize_ghosts( field )

      field_avg = field_avg + field

      where ( field > 1.0d0 )
        field = 1.0d0
      end where

      where ( field /= 1.0d0 )
        field = 0.0d0
      end where

      do iz=ra(3),rb(3)
        do iy=ra(2),rb(2)
          do ix=ra(1),rb(1)
            if ( sum(field( ix-1:ix+1, iy-1:iy+1, iz-1:iz+1 )) >= 2.d0 ) then
              field_color(ix,iy,iz) = dble(i+1)*4.0d0
            endif
          enddo
        enddo
      enddo

      i = i+1
    endif
  enddo
  close (14)


  call save_field_hdf5(0.d0, fname_avg, field_avg(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))
  call save_field_hdf5(0.d0, fname_color, field_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)))

  deallocate(field_avg)
  deallocate(field)


end subroutine uCT_assemble_HDF5
