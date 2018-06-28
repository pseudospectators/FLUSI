subroutine read_mesh_data(wings)

  use vars
  implicit none
  type (wing), dimension(1:nWings), intent (inout) :: wings
  character(len=strlen) :: data_file
  character(len=1)  :: wingstr
  integer :: i, j
  real(kind=pr), allocatable :: tmp

  do i=1, nWings

    !-- for naming files..
    write (wingstr,'(i1)') i

    ! Read initial coordinates x,y,z of all points in 2nd,3rd,4th columms
    ! respectively in the points_coor.t data file
    data_file = 'points_coor'//wingstr//'.t'
    call  read_mesh_data_2D_array(data_file, tmp)
    do j=1, maxval(tmp(:,1))
      wings(i)%x(j) = tmp(j,2)
      wings(i)%y(j) = tmp(j,3)
      wings(i)%z(j) = tmp(j,4)
    enddo

    ! Read indices of three vertices (correnponding to 3rd, 4th and 5tn columms)
    ! of all triangle elements of the mesh
    data_file = 'mesh_triangle_elements'//wingstr//'.t'
    call  read_mesh_data_2D_array(data_file, tmp)
    do j=1, size(tmp,DIM=1)
      wings(i)%tri_elements(j,1) = j
      wings(i)%tri_elements(j,2) = int(tmp(j,3))
      wings(i)%tri_elements(j,3) = int(tmp(j,4))
      wings(i)%tri_elements(j,4) = int(tmp(j,5))
    enddo

  enddo

end subroutine read_mesh_data

subroutine read_mesh_data_2D_array(data_file, data_2D_array)

character(len=strlen), intent(in) :: data_file
real(kind=pr), intent(inout) :: data_2D_array
integer :: num_lines, num_cols, n_header=0

call count_lines_in_ascii_file_mpi(data_file, num_lines, n_header)
call count_cols_in_ascii_file_mpi(data_file, num_cols, n_header)
allocate(data_2D_array(1:nlines, 1:num_cols) )
call read_array_from_ascii_file_mpi(data_file, data_2D_array, n_header)

end subroutine read_mesh_data_2D_array
