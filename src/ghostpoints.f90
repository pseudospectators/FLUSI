!-------------------------------------------------------------------------------
! FLUSI ghost point management
! allocates larger arrays with additional space for ghost points
! in real space and organizes their exchange between CPUs. We inherently rely
! on knowing every CPU the entire domain decomposition in a 2D array.
!
! 05 Oct 2014 (Thomas) Modify the ghost points to be added in all directions, also in
!       periodic ones. This makes applying finite differences easier and the
!       operation is not very expensive.
!
!-------------------------------------------------------------------------------
module ghosts

!-------------------------------------------------------------------------------
! interface for ghsot synchronization routines. you can
! call synchronize_ghosts(p)  (sync one field, or one component)
! call synchronize_ghosts(u,4)  (sync four components at the same time)
! note it is much more efficient to do the latter instead of 4 times the former
!-------------------------------------------------------------------------------
 interface synchronize_ghosts
   module procedure synchronize_ghosts, synchronize_ghosts_FD
 end interface


!!!!!!!!!!!!!!
 contains
!!!!!!!!!!!!!!


!-------------------------------------------------------------------------------
! Ghost point synchronization in all directions
! For only one 3d field
!-------------------------------------------------------------------------------
subroutine synchronize_ghosts ( fld )
  ! Routine assumes that heart of the matrix "fld" (thus the regular part of
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  implicit none
  ! Input/output
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))

  ! x direction is always local, i.e. it is not split among MPI procs. In this
  ! direction, which is assumed PERIODIC at present, we can locally copy data
  call synchronize_ghosts_FD_x_serial (fld,1) ! note nc=1


  ! y direction can be distributed or local, i.e. is is possible that it is split
  ! among processes.
  if (mpidims(2)>1) then
    call synchronize_ghosts_FD_y_mpi (fld,1)
  else
    call synchronize_ghosts_FD_y_serial (fld,1)
  endif


  ! z direction can be distributed or local
  ! p3dfft decomposes first in z, then in y
  if (mpidims(1)>1) then
    call synchronize_ghosts_FD_z_mpi (fld,1)
  else
    call synchronize_ghosts_FD_z_serial (fld,1)
  endif
end subroutine synchronize_ghosts



!-------------------------------------------------------------------------------
! Ghost point synchronization in all directions
! For nc 3d fields
!-------------------------------------------------------------------------------
subroutine synchronize_ghosts_FD ( fld, nc )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  integer, intent(in) :: nc
  real(kind=pr), intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)

  ! x direction is always local, i.e. it is not split among MPI procs. In this
  ! direction, which is assumed PERIODIC at present, we can locally copy data
  call synchronize_ghosts_FD_x_serial (fld,nc)


  ! y direction can be distributed or local
  if (mpidims(2)>1) then
    call synchronize_ghosts_FD_y_mpi (fld,nc)
  else
    call synchronize_ghosts_FD_y_serial (fld,nc)
  endif


  ! z direction can be distributed or local
  ! p3dfft decomposes first in z, then in y
  if (mpidims(1)>1) then
    call synchronize_ghosts_FD_z_mpi (fld,nc)
  else
    call synchronize_ghosts_FD_z_serial (fld,nc)
  endif
end subroutine synchronize_ghosts_FD

!-------------------------------------------------------------------------------

! Ghost point synchronization in z, distributed
subroutine synchronize_ghosts_FD_z_mpi ( fld,nc )
  use vars
  use mpi
  implicit none
  ! Input/output
  integer,intent(in) :: nc
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  ! Local variables
  integer :: statu(MPI_STATUS_SIZE),nnl,nnr,disp,dir,source,dest,mpicode
  real(kind=pr) :: fld_send_l(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1),1:nc),&
                   fld_send_r(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3),1:nc)
  real(kind=pr) :: fld_recv_l(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1),1:nc),&
                   fld_recv_r(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3),1:nc)
  ! Size of buffer arrays
  nnl = (ra(3)-ga(3))*(gb(2)-ga(2)+1)*(gb(1)-ga(1)+1)*nc
  nnr = (gb(3)-rb(3))*(gb(2)-ga(2)+1)*(gb(1)-ga(1)+1)*nc
  ! Copy data to buffer arrays
  fld_send_l(:,:,:,1:nc) = fld(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1),1:nc)
  fld_send_r(:,:,:,1:nc) = fld(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3),1:nc)
  ! Find rank of neighbors on the right
  disp=1                 !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommz,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the right
  call MPI_SENDRECV(fld_send_r,nnr,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_l,nnr,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommz,statu,mpicode)
  ! Find rank of neighbors on the left
  disp=-1                !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommz,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the left
  call MPI_SENDRECV(fld_send_l,nnl,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_r,nnl,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommz,statu,mpicode)
  ! Copy data to output array
  fld(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1),1:nc) = fld_recv_l(:,:,:,1:nc)
  fld(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3),1:nc) = fld_recv_r(:,:,:,1:nc)
end subroutine synchronize_ghosts_FD_z_mpi

!-------------------------------------------------------------------------------

! Ghost point synchronization in y
subroutine synchronize_ghosts_FD_y_mpi ( fld,nc )
  use vars
  use mpi
  implicit none
  ! Input/output
  integer,intent(in) :: nc
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  ! Local variables
  integer :: statu(MPI_STATUS_SIZE),nnl,nnr,disp,dir,source,dest,mpicode
  real(kind=pr) :: fld_send_l(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3),1:nc),&
                   fld_send_r(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3),1:nc)
  real(kind=pr) :: fld_recv_l(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3),1:nc),&
                   fld_recv_r(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3),1:nc)
  ! Size of buffer arrays
  nnl = (ra(2)-ga(2))*(gb(3)-ga(3)+1)*(gb(1)-ga(1)+1)*nc
  nnr = (gb(2)-rb(2))*(gb(3)-ga(3)+1)*(gb(1)-ga(1)+1)*nc
  ! Copy data to buffer arrays
  fld_send_l(:,:,:,1:nc) = fld(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3),1:nc)
  fld_send_r(:,:,:,1:nc) = fld(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3),1:nc)
  ! Find rank of neighbors on the right
  disp=1                 !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommy,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the right
  call MPI_SENDRECV(fld_send_r,nnr,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_l,nnr,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommy,statu,mpicode)
  ! Find rank of neighbors on the left
  disp=-1                !immediate neighbors
  dir=0
  call MPI_CART_SHIFT(mpicommy,dir,disp,source,dest,mpicode)
  ! Send data to the neighbors on the left
  call MPI_SENDRECV(fld_send_l,nnl,MPI_DOUBLE_PRECISION,dest,0,&
                    fld_recv_r,nnl,MPI_DOUBLE_PRECISION,source,0,&
                    mpicommy,statu,mpicode)
  ! Copy data to output array
  fld(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3),1:nc) = fld_recv_l(:,:,:,1:nc)
  fld(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3),1:nc) = fld_recv_r(:,:,:,1:nc)
end subroutine synchronize_ghosts_FD_y_mpi

!-------------------------------------------------------------------------------

! Ghost point synchronization in x. This is always local
subroutine synchronize_ghosts_FD_x_serial ( fld,nc )
  use vars
  implicit none
  ! Input/output
  integer,intent(in) :: nc
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  ! Copy data to ghost points. X direction is local
  ! Do nothing if there are no ghost points
  if ((ga(1)<ra(1)).and.(gb(1)>rb(1))) then
    fld((rb(1)+1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc) = fld(ra(1):(2*ra(1)-ga(1)-1),ga(2):gb(2),ga(3):gb(3),1:nc)
    fld(ga(1):(ra(1)-1),ga(2):gb(2),ga(3):gb(3),1:nc) = fld((2*rb(1)-gb(1)+1):rb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  endif
end subroutine synchronize_ghosts_FD_x_serial

!-------------------------------------------------------------------------------

! Ghost point synchronization in y on one process
subroutine synchronize_ghosts_FD_y_serial ( fld,nc )
  use vars
  implicit none
  ! Input/output
  integer,intent(in) :: nc
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  ! Copy data to ghost points. Y direction is local in this case
  ! Do nothing if there are no ghost points
  if ((ga(2)<ra(2)).and.(gb(2)>rb(2))) then
    fld(ga(1):gb(1),(rb(2)+1):gb(2),ga(3):gb(3),1:nc) = fld(ga(1):gb(1),ra(2):(2*ra(2)-ga(2)-1),ga(3):gb(3),1:nc)
    fld(ga(1):gb(1),ga(2):(ra(2)-1),ga(3):gb(3),1:nc) = fld(ga(1):gb(1),(2*rb(2)-gb(2)+1):rb(2),ga(3):gb(3),1:nc)
  endif
end subroutine synchronize_ghosts_FD_y_serial

!-------------------------------------------------------------------------------

! Ghost point synchronization in z on one process
subroutine synchronize_ghosts_FD_z_serial ( fld,nc )
  use vars
  implicit none
  ! Input/output
  integer,intent(in) :: nc
  real(kind=pr),intent(inout) :: fld(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nc)
  ! Copy data to ghost points. Z direction is local in this case
  ! Do nothing if there are no ghost points
  if ((ga(3)<ra(3)).and.(gb(3)>rb(3))) then
    fld(ga(1):gb(1),ga(2):gb(2),(rb(3)+1):gb(3),1:nc) = fld(ga(1):gb(1),ga(2):gb(2),ra(3):(2*ra(3)-ga(3)-1),1:nc)
    fld(ga(1):gb(1),ga(2):gb(2),ga(3):(ra(3)-1),1:nc) = fld(ga(1):gb(1),ga(2):gb(2),(2*rb(3)-gb(3)+1):rb(3),1:nc)
  endif
end subroutine synchronize_ghosts_FD_z_serial

end module ghosts
