subroutine SaveWingData( time, wings )
  use vars
  implicit none
  real (kind=pr), intent (in) :: time
  type(flexible_wing), dimension(1:nWings), intent (inout) :: wings
  character(len=16) :: format_ns1
  character(len=3)  :: ns1_string, ctrlpoint
  character(len=1)  :: wingstr

  integer :: np,step,i

  !-- loop over wings
  do i=1, nWings
    !-- for naming files..
    write (wingstr,'(i1)') i
    write (ctrlpoint,'(i3)') wings(i)%ControlPoint

    !Get total number of points
    np = wings(i)%np

    ! set up formats
      write(ns1_string,'(i3)') np+1
      format_ns1 = '('//ns1_string//'(es15.8,1x))'

      !-- save control point data
      if (wings(i)%ControlPoint /= 0) then
      open  (14, file = 'wing_data'//wingstr//'_point'//ctrlpoint//'.t', status = 'unknown',position='append')
      write (14, '(21(es15.8,1x))') &
        time,&
        wings(i)%x(wings(i)%ControlPoint),&
        wings(i)%y(wings(i)%ControlPoint),&
        wings(i)%z(wings(i)%ControlPoint),&
        wings(i)%vx(wings(i)%ControlPoint),&
        wings(i)%vy(wings(i)%ControlPoint),&
        wings(i)%vz(wings(i)%ControlPoint)
      close (14)
    endif

    open (14, file = 'wing_x'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%x(1:np)
    close (14)

    open (14, file = 'wing_y'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%y(1:np)
    close (14)

    open (14, file = 'wing_z'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%z(1:np)
    close (14)

    open (14, file = 'wing_xrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(0*np+1:1*np)
    close (14)

    open (14, file = 'wing_yrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(1*np+1:2*np)
    close (14)

    open (14, file = 'wing_zrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(2*np+1:3*np)
    close (14)

    open (14, file = 'wing_vxrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(3*np+1:4*np)
    close (14)

    open (14, file = 'wing_vyrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(4*np+1:5*np)
    close (14)

    open (14, file = 'wing_vzrel'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%u_old(5*np+1:6*np)
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

    open (14, file = 'Fext_x'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%Fext(1:np)
    close (14)

    open (14, file = 'Fext_y'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%Fext(np+1:2*np)
    close (14)

    open (14, file = 'Fext_z'//wingstr//'.t', status = 'unknown',position='append')
    write (14, format_ns1) time, wings(i)%Fext(2*np+1:3*np)
    close (14)

  enddo


end subroutine SaveWingData

!-------------------------------------------------------------------------
! dump runtime backup for the flexible wing solver
!-------------------------------------------------------------------------
subroutine dump_flexible_wing_backup( time, wings, nbackup )
  use vars
  implicit none

  real(kind=pr), intent(in) :: time
  type(flexible_wing), dimension (1:nWings), intent (in) :: wings
  integer, intent (in):: nbackup
  integer :: i
  character(len=24) :: filename

  !-- only root rank dumps backup
  if (root) then
    write(*,'(A)',advance='no') "Backuping flexible wing solver..."
    write(filename,'("runtime_backup",i1,".fsi_bckp")') nbackup
    write(*,'(A)',advance='no') "file="//filename
    write(*,'(" time=",e11.4)',advance='no') time


    open(14,file=filename,status='replace',form='formatted')


    write(14,*) time
    write(14,*) nWings
    write(14,*) nVeins, nVeins_BC, nMembranes, nMembrane_edges

    do i=1,nWings
      write(14,*) wings(i)%np
      write(14,*) wings(i)%x, wings(i)%y, wings(i)%z
      write(14,*) wings(i)%vx, wings(i)%vy, wings(i)%vz
      write(14,*) wings(i)%u_old, wings(i)%u_oldold
      write(14,*) wings(i)%x0, wings(i)%y0, wings(i)%z0
      write(14,*) wings(i)%WingAngle_x, wings(i)%WingAngle_y, wings(i)%WingAngle_z
      write(14,*) wings(i)%vt0, wings(i)%at0
      write(14,*) wings(i)%vr0, wings(i)%ar0
      write(14,*) wings(i)%StartupStep
    enddo

    close(14)
    write(*,'(A)',advance='yes') "...DONE!"
  endif

end subroutine dump_flexible_wing_backup

!-------------------------------------------------------------------------
! read runtime backup for the flexible wing solver
!-------------------------------------------------------------------------
subroutine read_flexible_wing_backup( wings, filename )
  use vars
  implicit none

  type(flexible_wing), dimension (1:nWings), intent (inout) :: Wings
  integer :: nwings_file, nVeins_file, nVeins_BC_file, nMembranes_file, nMembrane_edges_file, i, mpicode
  character(len=*), intent(in) :: filename
  real(kind=pr) :: time

  ! only root shall read in the file, the results is then send to the other
  if (root) then
    write(*,'(A)') "Reading in backup of flexible wing solver: "//filename

    open(14,file=filename,status='old',form='formatted',action='read')
    read(14,*) time
    read(14,*) nwings_file
    read(14,*) nVeins_file, nVeins_BC_file, nMembranes_file, nMembrane_edges_file

    write(*,'("(time=",e11.4,")")') time

    if ((nwings_file/=nwings).or.(nVeins_file/=nVeins)) then
      call abort(95,"Cant resume flexible_wing backup: wing mesh doesnt match")
    endif

    do i =1, nwings
      read(14,*) wings(i)%np
      read(14,*) wings(i)%x, wings(i)%y, wings(i)%z
      read(14,*) wings(i)%vx, wings(i)%vy, wings(i)%vz
      read(14,*) wings(i)%u_old, wings(i)%u_oldold
      read(14,*) wings(i)%x0, wings(i)%y0, wings(i)%z0
      read(14,*) wings(i)%WingAngle_x, wings(i)%WingAngle_y, wings(i)%WingAngle_z
      read(14,*) wings(i)%vt0, wings(i)%at0
      read(14,*) wings(i)%vr0, wings(i)%ar0
      read(14,*) wings(i)%StartupStep
    enddo
    close(14)

    write(*,'(A)') "...done reading flexible wing backup."
  endif

  do i =1, nwings
    call MPI_Bcast( wings(i)%np, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%y, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%z, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vx, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vy, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vz, npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%u_old, 6*npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%u_oldold, 6*npmax, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%x0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%y0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%z0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%WingAngle_x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%WingAngle_y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%WingAngle_z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vt0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%at0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%vr0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%ar0, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
    call MPI_Bcast( wings(i)%StartupStep, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpicode )

  enddo
end subroutine read_flexible_wing_backup
