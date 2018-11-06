subroutine SaveWingData( time, wings )
  use vars
  implicit none
  real (kind=pr), intent (in) :: time
  type (wing), dimension(1:nWings), intent (inout) :: wings
  character(len=16) :: format_ns1
  character(len=3)  :: ns1_string
  character(len=1)  :: wingstr
  integer :: np,step,i

  !-- loop over wings
  do i=1, nWings
    !-- for naming files..
    write (wingstr,'(i1)') i

    !Get total number of points
    np = wings(i)%np

    ! set up formats
      write(ns1_string, '(I3)') np+1
      format_ns1 = '('//ns1_string//'(es15.8,1x))'

    !-- save trailing edge data
    !open  (14, file = 'wing_data'//wingstr//'.t', status = 'unknown',position='append')
    !write (14, '(21(es15.8,1x))') &
    !  time,&
    !  wings(i)%x(np),&
    !  wings(i)%y(np),&
    !  wings(i)%z(np),&
    !  wings(i)%vx(np),&
    !  wings(i)%vy(np),&
    !  wings(i)%vz(np)
    !close (14)

    open (14, file = 'wing_x'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%x(1:np)
    close (14)

    open (14, file = 'wing_y'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%y(1:np)
    close (14)

    open (14, file = 'wing_z'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%z(1:np)
    close (14)

    open (14, file = 'wing_vx'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%vx(1:np)
    close (14)

    open (14, file = 'wing_vy'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%vy(1:np)
    close (14)

    open (14, file = 'wing_vz'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%vz(1:np)
    close (14)

  enddo


end subroutine SavewingData
