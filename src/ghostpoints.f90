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
! TODO:
!       This must be generalized to arbitrary decompositions (e.g. 1D in 
!       x-direction, 3D and so on)
!-------------------------------------------------------------------------------
module ghosts
 
 interface synchronize_ghosts
   module procedure synchronize_ghosts, synchronize_ghosts_FD
 end interface

 
!!!!!!!!!!!!!! 
 contains 
!!!!!!!!!!!!!!
 
 
subroutine synchronize_ghosts ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  if ((decomposition=="1D").and.(ng>0)) then
    call sync_ghosts_1D_z ( field )
  elseif ((decomposition=="2D").and.(ng>0)) then
    call sync_ghosts_2D_yz ( field )
  else
    if(mpirank==0) write(*,*) "Error in ghost_sync (choice!)"
    call abort()
  endif
end subroutine synchronize_ghosts

!-------------------------------------------------------------------------------

subroutine sync_ghosts_1D_z ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer :: mpicode, destination, origin, status(MPI_status_size)
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)
  
  ! in the 1D decomposition case, we have two periodic borders. They are not 
  ! chunked so their dimensions are indeed 0:nx-1 and 0:ny-1 on all CPU
  ! x-direction:
  field(-ng:-1,:,:)=field(nx-1-(ng-1):nx-1,:,:)
  field(nx:nx+ng-1,:,:)=field(0:ng-1,:,:)
  
  ! y-direction:
  field(:,-ng:-1,:)=field(:,ny-1-(ng-1):ny-1,:)  
  field(:,ny:ny+ng-1,:)=field(:,0:ng-1,:)
  
  ! for scalar runs, we do not need mpi sync...
  if(mpisize==1) then
    field(:,:,nz:nz+ng-1)=field(:,:,0:ng-1)
    field(:,:,-ng:-1)=field(:,:,nz-1-(ng-1):nz-1)
    return
  endif
  
  ! in the third (z) direction, the data is spread accros processes
  destination = per( mpirank+1, mpisize )             ! send to your right
  origin      = per( mpirank-1, mpisize )             ! get from your left
  call MPI_sendrecv( &
        field(ra(1):rb(1),ra(2):rb(2),izmax-ng+1:izmax),&    ! send buffer
        nx*ny*ng,&                        ! send buffer size
        MPI_DOUBLE_PRECISION,&
        destination,&                     ! whom to send it to
        mpirank,&                         ! send tag
        field(ra(1):rb(1),ra(2):rb(2),izmin-1-ng+1:izmin-1),& ! receive buffer
        nx*ny*ng,&                        ! recvcount
        MPI_DOUBLE_PRECISION,&
        origin,&                          ! source
        origin,&                          ! recv tag
        MPI_COMM_WORLD,status,mpicode)
        
  destination = per( mpirank-1, mpisize ) ! send to your left
  origin      = per( mpirank+1, mpisize ) ! get from your right
  call MPI_sendrecv( &
        field(ra(1):rb(1),ra(2):rb(2),izmin:izmin+ng-1),&    ! send buffer
        nx*ny*ng,&                        ! send buffer size
        MPI_DOUBLE_PRECISION,&
        destination,&                     ! whom to send it to
        mpirank,&                         ! send tag
        field(ra(1):rb(1),ra(2):rb(2),izmax+1:izmax+1+ng-1),& ! receive buffer
        nx*ny*ng,&                        ! recvcount
        MPI_DOUBLE_PRECISION,&
        origin,&                          ! source
        origin,&                          ! recv tag
        MPI_COMM_WORLD,status,mpicode)     
end subroutine sync_ghosts_1D_z


!-------------------------------------------------------------------------------

subroutine sync_ghosts_2D_yz ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer :: mpicode, destination, origin, status(MPI_status_size)
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)

  ! in the 2D decomposition case, we have only one periodic direction, namely
  ! the x-one. Note that on all CPU, this had indeed dimensions 0:nx-1, since it
  ! is not chunked.
  field(nx:nx+ng-1,:,:)=field(0:ng-1,:,:)
  field(-ng:-1,:,:)=field(nx-1-(ng-1):nx-1,:,:)
  
  destination = yz_plane_ranks(iymin,per(izmax+1,nz)) ! send to your right
  origin      = yz_plane_ranks(iymin,per(izmin-1,nz)) ! get from your left
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymax,izmax-ng+1:izmax),&   ! send buffer
                    nx*(iymax-iymin+1)*ng,&                    ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                              ! whom to send it to
                    mpirank,&                                  ! send tag
                    field(ra(1):rb(1),iymin:iymax,izmin-1-ng+1:izmin-1),&! receive buffer
                    nx*(iymax-iymin+1)*ng,&                    ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                   ! source
                    origin,&                                   ! recv tag
                    MPI_COMM_WORLD,status,mpicode)
                    
  destination = yz_plane_ranks(iymin,per(izmin-1,nz)) ! send to your left
  origin      = yz_plane_ranks(iymin,per(izmax+1,nz)) ! get from your right
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymax,izmin:izmin+ng-1),&  ! send buffer
                    nx*(iymax-iymin+1)*ng,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymin:iymax,izmax+1:izmax+1+ng-1),&! receive buffer
                    nx*(iymax-iymin+1)*ng,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)     
                    
  destination = yz_plane_ranks(per(iymax+1,ny),izmin)         ! send to your top
  origin      = yz_plane_ranks(per(iymin-1,ny),izmin)         ! get from your bottom
  call MPI_sendrecv( field(ra(1):rb(1),iymax-ng+1:iymax,izmin:izmax),&  ! send buffer
                    nx*(izmax-izmin+1)*ng,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymin-1-ng+1:iymin-1,izmin:izmax),&! receive buffer
                    nx*(izmax-izmin+1)*ng,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)
                    
  destination = yz_plane_ranks(per(iymin-1,ny),izmin) ! send to your bottom
  origin      = yz_plane_ranks(per(iymax+1,ny),izmin) ! get from your top
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymin+ng-1,izmin:izmax),&  ! send buffer
                    nx*(izmax-izmin+1)*ng,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymax+1:iymax+1+ng-1,izmin:izmax),&! receive buffer
                    nx*(izmax-izmin+1)*ng,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)  
                      
    !-- corners   
       
!     destination = yz_plane_ranks(per(iymax+1,ny),per(izmax+1,nz)) ! send to top right
!     origin      = yz_plane_ranks(per(iymin-1,ny),per(izmin-1,nz)) ! get from bottom left
!     call MPI_sendrecv( field(ra(1):rb(1),iymax-ng+1:iymax,izmax-ng+1:izmax),& ! send buffer
!                       nx*ng*ng,& ! send buffer size
!                       MPI_DOUBLE_PRECISION,&
!                       destination,& ! whom to send it to
!                       mpirank,& ! send tag
!                       field(ra(1):rb(1),iymin-1-ng+1:iymin-1,izmin-1-ng+1:izmin-1),&! receive buffer
!                       nx*ng*ng,& ! recvcount
!                       MPI_DOUBLE_PRECISION,&
!                       origin,& ! source
!                       origin,& ! recv tag
!                       MPI_COMM_WORLD,status,mpicode)   
!                       
!     destination = yz_plane_ranks(per(iymin-1,ny),per(izmin-1,nz)) ! send to bottom left
!     origin      = yz_plane_ranks(per(iymax+1,ny),per(izmax+1,nz)) ! get from top right
!     call MPI_sendrecv( field(ra(1):rb(1),iymin:iymin+ng-1,izmin:izmin+ng-1),& ! send buffer
!                       nx*ng*ng,& ! send buffer size
!                       MPI_DOUBLE_PRECISION,&
!                       destination,& ! whom to send it to
!                       mpirank,& ! send tag
!                       field(ra(1):rb(1),iymax+1:iymax+1+ng-1,izmax+1:izmax+1+ng-1),&! receive buffer
!                       nx*ng*ng,& ! recvcount
!                       MPI_DOUBLE_PRECISION,&
!                       origin,& !source
!                       origin,& ! recv tag
!                       MPI_COMM_WORLD,status,mpicode)                  
!                       
!     destination = yz_plane_ranks(per(iymax+1,ny),per(izmin-1,nz))
!     origin      = yz_plane_ranks(per(iymin-1,ny),per(izmax+1,nz)) 
!     call MPI_sendrecv( field(ra(1):rb(1),iymax-ng+1:iymax,izmin:izmin+ng-1),& ! send buffer
!                       nx*ng*ng,& ! send buffer size
!                       MPI_DOUBLE_PRECISION,&
!                       destination,& ! whom to send it to
!                       mpirank,& ! send tag
!                       field(ra(1):rb(1),iymin-1-ng+1:iymin-1,izmax+1:izmax+1+ng-1),&! receive buffer
!                       nx*ng*ng,& ! recvcount
!                       MPI_DOUBLE_PRECISION,&
!                       origin,& !source
!                       origin,& ! recv tag
!                       MPI_COMM_WORLD,status,mpicode)       
!                       
!     destination = yz_plane_ranks(per(iymin-1,ny),per(izmax+1,nz)) 
!     origin      = yz_plane_ranks(per(iymax+1,ny),per(izmin-1,nz))
!     call MPI_sendrecv( field(ra(1):rb(1),iymin:iymin+ng-1,izmax-ng+1:izmax),& ! send buffer
!                       nx*ng*ng,& ! send buffer size
!                       MPI_DOUBLE_PRECISION,&
!                       destination,& ! whom to send it to
!                       mpirank,& ! send tag
!                       field(ra(1):rb(1),iymax+1:iymax+1+ng-1,izmin-1-ng+1:izmin-1),&! receive buffer
!                       nx*ng*ng,& ! recvcount
!                       MPI_DOUBLE_PRECISION,&
!                       origin,& !source
!                       origin,& ! recv tag
!                       MPI_COMM_WORLD,status,mpicode)                    
  
end subroutine sync_ghosts_2D_yz






subroutine synchronize_ghosts_FD ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  
  if ((decomposition=="1D").and.(ng>0)) then
    call sync_ghosts_1D_z_FD ( field )
  elseif ((decomposition=="2D").and.(ng>0)) then
    call sync_ghosts_2D_yz_FD ( field )
  else
    if(mpirank==0) write(*,*) "Error in ghost_sync (choice!)"
    call abort()
  endif
end subroutine synchronize_ghosts_FD

!-------------------------------------------------------------------------------

subroutine sync_ghosts_1D_z_FD ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer :: mpicode, destination, origin, status(MPI_status_size)
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)
  
  ! in the 1D decomposition case, we have two periodic borders. They are not 
  ! chunked so their dimensions are indeed 0:nx-1 and 0:ny-1 on all CPU
  ! x-direction:
  if (nx>1) then
    field(-ng:-1,:,:,1:neq)=field(nx-1-(ng-1):nx-1,:,:,1:neq)
    field(nx:nx+ng-1,:,:,1:neq)=field(0:ng-1,:,:,1:neq)
  endif
  
  ! y-direction:
  field(:,-ng:-1,:,1:neq)=field(:,ny-1-(ng-1):ny-1,:,1:neq)  
  field(:,ny:ny+ng-1,:,1:neq)=field(:,0:ng-1,:,1:neq)
  
  ! for scalar runs, we do not need mpi sync...
  if(mpisize==1) then
    field(:,:,nz:nz+ng-1,1:neq)=field(:,:,0:ng-1,1:neq)
    field(:,:,-ng:-1,1:neq)=field(:,:,nz-1-(ng-1):nz-1,1:neq)
    return
  endif
  
  ! in the third (z) direction, the data is spread accros processes
  destination = per( mpirank+1, mpisize )             ! send to your right
  origin      = per( mpirank-1, mpisize )             ! get from your left
  call MPI_sendrecv( &
        field(ra(1):rb(1),ra(2):rb(2),izmax-ng+1:izmax,1:neq),&    ! send buffer
        nx*ny*ng*neq,&                        ! send buffer size
        MPI_DOUBLE_PRECISION,&
        destination,&                     ! whom to send it to
        mpirank,&                         ! send tag
        field(ra(1):rb(1),ra(2):rb(2),izmin-1-ng+1:izmin-1,1:neq),& ! receive buffer
        nx*ny*ng*neq,&                        ! recvcount
        MPI_DOUBLE_PRECISION,&
        origin,&                          ! source
        origin,&                          ! recv tag
        MPI_COMM_WORLD,status,mpicode)
        
  destination = per( mpirank-1, mpisize ) ! send to your left
  origin      = per( mpirank+1, mpisize ) ! get from your right
  call MPI_sendrecv( &
        field(ra(1):rb(1),ra(2):rb(2),izmin:izmin+ng-1,1:neq),&    ! send buffer
        nx*ny*ng*neq,&                        ! send buffer size
        MPI_DOUBLE_PRECISION,&
        destination,&                     ! whom to send it to
        mpirank,&                         ! send tag
        field(ra(1):rb(1),ra(2):rb(2),izmax+1:izmax+1+ng-1,1:neq),& ! receive buffer
        nx*ny*ng*neq,&                        ! recvcount
        MPI_DOUBLE_PRECISION,&
        origin,&                          ! source
        origin,&                          ! recv tag
        MPI_COMM_WORLD,status,mpicode)     
end subroutine sync_ghosts_1D_z_FD


!-------------------------------------------------------------------------------

subroutine sync_ghosts_2D_yz_FD ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer :: mpicode, destination, origin, status(MPI_status_size)
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)

  if (nx>1) then
    ! in the 2D decomposition case, we have only one periodic direction, namely
    ! the x-one. Note that on all CPU, this had indeed dimensions 0:nx-1, since it
    ! is not chunked.
    field(nx:nx+ng-1,:,:,1:neq)=field(0:ng-1,:,:,1:neq)
    field(-ng:-1,:,:,1:neq)=field(nx-1-(ng-1):nx-1,:,:,1:neq)
  endif
  
  destination = yz_plane_ranks(iymin,per(izmax+1,nz)) ! send to your right
  origin      = yz_plane_ranks(iymin,per(izmin-1,nz)) ! get from your left
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymax,izmax-ng+1:izmax,1:neq),&   ! send buffer
                    nx*(iymax-iymin+1)*ng*neq,&                    ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                              ! whom to send it to
                    mpirank,&                                  ! send tag
                    field(ra(1):rb(1),iymin:iymax,izmin-1-ng+1:izmin-1,1:neq),&! receive buffer
                    nx*(iymax-iymin+1)*ng*neq,&                    ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                   ! source
                    origin,&                                   ! recv tag
                    MPI_COMM_WORLD,status,mpicode)
                    
  destination = yz_plane_ranks(iymin,per(izmin-1,nz)) ! send to your left
  origin      = yz_plane_ranks(iymin,per(izmax+1,nz)) ! get from your right
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymax,izmin:izmin+ng-1,1:neq),&  ! send buffer
                    nx*(iymax-iymin+1)*ng*neq,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymin:iymax,izmax+1:izmax+1+ng-1,1:neq),&! receive buffer
                    nx*(iymax-iymin+1)*ng*neq,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)     
                    
  destination = yz_plane_ranks(per(iymax+1,ny),izmin)         ! send to your top
  origin      = yz_plane_ranks(per(iymin-1,ny),izmin)         ! get from your bottom
  call MPI_sendrecv( field(ra(1):rb(1),iymax-ng+1:iymax,izmin:izmax,1:neq),&  ! send buffer
                    nx*(izmax-izmin+1)*ng*neq,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymin-1-ng+1:iymin-1,izmin:izmax,1:neq),&! receive buffer
                    nx*(izmax-izmin+1)*ng*neq,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)
                    
  destination = yz_plane_ranks(per(iymin-1,ny),izmin) ! send to your bottom
  origin      = yz_plane_ranks(per(iymax+1,ny),izmin) ! get from your top
  call MPI_sendrecv( field(ra(1):rb(1),iymin:iymin+ng-1,izmin:izmax,1:neq),&  ! send buffer
                    nx*(izmax-izmin+1)*ng*neq,&                   ! send buffer size
                    MPI_DOUBLE_PRECISION,&
                    destination,&                             ! whom to send it to
                    mpirank,&                                 ! send tag
                    field(ra(1):rb(1),iymax+1:iymax+1+ng-1,izmin:izmax,1:neq),&! receive buffer
                    nx*(izmax-izmin+1)*ng*neq,&                   ! recvcount
                    MPI_DOUBLE_PRECISION,&
                    origin,&                                  ! source
                    origin,&                                  ! recv tag
                    MPI_COMM_WORLD,status,mpicode)  
end subroutine sync_ghosts_2D_yz_FD

end module ghosts