!-------------------------------------------------------------------------------
! FLUSI ghost point management
! allocates larger arrays with additional space for ghost points
! in real space and organizes their exchange between CPUs. We inherently rely
! on knowing every CPU the entire domain decomposition in a 2D array. 
!-------------------------------------------------------------------------------


subroutine synchronize_ghosts ( field )
  ! Routine assumes that heart of the matrix "field" (thus the regular part of 
  ! it) have been filled previously. Here, we exchange only the ghosts
  use vars
  use mpi
  implicit none
  real(kind=pr), intent(inout) :: field(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  character(len=3)::str
  integer :: mpicode, destination, origin, status
  integer :: ixmin,ixmax,iymin,iymax,izmin,izmax
  integer :: iy
  ixmin=ra(1)
  ixmax=rb(1)
  iymin=ra(2)
  iymax=rb(2)
  izmin=ra(3)
  izmax=rb(3)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (decomposition=="1D") then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    destination = per( mpirank+1, mpisize )             ! send to your right
    origin      = per( mpirank-1, mpisize )             ! get from your left
    call MPI_sendrecv( field(:,:,izmax-ng+1:izmax),&    ! send buffer
                      nx*ny*ng,&                        ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                     ! whom to send it to
                      mpirank,&                         ! send tag
                      field(:,:,izmin-1-ng+1:izmin-1),& ! receive buffer
                      nx*ny*ng,&                        ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                          ! source
                      origin,&                          ! recv tag
                      MPI_COMM_WORLD,status,mpicode)
                      
    destination = per( mpirank-1, mpisize ) ! send to your left
    origin      = per( mpirank+1, mpisize ) ! get from your right
    call MPI_sendrecv( field(:,:,izmin:izmin+ng-1),&    ! send buffer
                      nx*ny*ng,&                        ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                     ! whom to send it to
                      mpirank,&                         ! send tag
                      field(:,:,izmax+1:izmax+1-ng+1),& ! receive buffer
                      nx*ny*ng,&                        ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                          ! source
                      origin,&                          ! recv tag
                      MPI_COMM_WORLD,status,mpicode)                 
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  elseif (decomposition=="2D") then
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    destination = yz_plane_ranks(iymin,per(izmax+1,nz)) ! send to your right
    origin      = yz_plane_ranks(iymin,per(izmin-1,nz)) ! get from your left
    call MPI_sendrecv( field(:,iymin:iymax,izmax-ng+1:izmax),&   ! send buffer
                      nx*(iymax-iymin+1)*ng,&                    ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                              ! whom to send it to
                      mpirank,&                                  ! send tag
                      field(:,iymin:iymax,izmin-1-ng+1:izmin-1),&! receive buffer
                      nx*(iymax-iymin+1)*ng,&                    ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                                   ! source
                      origin,&                                   ! recv tag
                      MPI_COMM_WORLD,status,mpicode)
                      
    destination = yz_plane_ranks(iymin,per(izmin-1,nz)) ! send to your left
    origin      = yz_plane_ranks(iymin,per(izmax+1,nz)) ! get from your right
    call MPI_sendrecv( field(:,iymin:iymax,izmin:izmin+ng-1),&  ! send buffer
                      nx*(iymax-iymin+1)*ng,&                   ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                             ! whom to send it to
                      mpirank,&                                 ! send tag
                      field(:,iymin:iymax,izmax+1:izmax+1+ng-1),&! receive buffer
                      nx*(iymax-iymin+1)*ng,&                   ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                                  ! source
                      origin,&                                  ! recv tag
                      MPI_COMM_WORLD,status,mpicode)     
                      
    destination = yz_plane_ranks(per(iymax+1,ny),izmin)         ! send to your top
    origin      = yz_plane_ranks(per(iymin-1,ny),izmin)         ! get from your bottom
    call MPI_sendrecv( field(:,iymax-ng+1:iymax,izmin:izmax),&  ! send buffer
                      nx*(izmax-izmin+1)*ng,&                   ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                             ! whom to send it to
                      mpirank,&                                 ! send tag
                      field(:,iymin-1-ng+1:iymin-1,izmin:izmax),&! receive buffer
                      nx*(izmax-izmin+1)*ng,&                   ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                                  ! source
                      origin,&                                  ! recv tag
                      MPI_COMM_WORLD,status,mpicode)
                      
    destination = yz_plane_ranks(per(iymin-1,ny),izmin) ! send to your bottom
    origin      = yz_plane_ranks(per(iymax+1,ny),izmin) ! get from your top
    call MPI_sendrecv( field(:,iymin:iymin+ng-1,izmin:izmax),&  ! send buffer
                      nx*(izmax-izmin+1)*ng,&                   ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,&                             ! whom to send it to
                      mpirank,&                                 ! send tag
                      field(:,iymax+1:iymax+1+ng-1,izmin:izmax),&! receive buffer
                      nx*(izmax-izmin+1)*ng,&                   ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,&                                  ! source
                      origin,&                                  ! recv tag
                      MPI_COMM_WORLD,status,mpicode)  
                      
    !-- corners   
       
    destination = yz_plane_ranks(per(iymax+1,ny),per(izmax+1,nz)) ! send to top right
    origin      = yz_plane_ranks(per(iymin-1,ny),per(izmin-1,nz)) ! get from bottom left
    call MPI_sendrecv( field(:,iymax-ng+1:iymax,izmax-ng+1:izmax),& ! send buffer
                      nx*ng*ng,& ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,& ! whom to send it to
                      mpirank,& ! send tag
                      field(:,iymin-1-ng+1:iymin-1,izmin-1-ng+1:izmin-1),&! receive buffer
                      nx*ng*ng,& ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,& ! source
                      origin,& ! recv tag
                      MPI_COMM_WORLD,status,mpicode)   
                      
    destination = yz_plane_ranks(per(iymin-1,ny),per(izmin-1,nz)) ! send to bottom left
    origin      = yz_plane_ranks(per(iymax+1,ny),per(izmax+1,nz)) ! get from top right
    call MPI_sendrecv( field(:,iymin:iymin+ng-1,izmin:izmin+ng-1),& ! send buffer
                      nx*ng*ng,& ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,& ! whom to send it to
                      mpirank,& ! send tag
                      field(:,iymax+1:iymax+1+ng-1,izmax+1:izmax+1+ng-1),&! receive buffer
                      nx*ng*ng,& ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,& !source
                      origin,& ! recv tag
                      MPI_COMM_WORLD,status,mpicode)                  
                      
    destination = yz_plane_ranks(per(iymax+1,ny),per(izmin-1,nz))
    origin      = yz_plane_ranks(per(iymin-1,ny),per(izmax+1,nz)) 
    call MPI_sendrecv( field(:,iymax-ng+1:iymax,izmin:izmin+ng-1),& ! send buffer
                      nx*ng*ng,& ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,& ! whom to send it to
                      mpirank,& ! send tag
                      field(:,iymin-1-ng+1:iymin-1,izmax+1:izmax+1+ng-1),&! receive buffer
                      nx*ng*ng,& ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,& !source
                      origin,& ! recv tag
                      MPI_COMM_WORLD,status,mpicode)       
                      
    destination = yz_plane_ranks(per(iymin-1,ny),per(izmax+1,nz)) 
    origin      = yz_plane_ranks(per(iymax+1,ny),per(izmin-1,nz))
    call MPI_sendrecv( field(:,iymin:iymin+ng-1,izmax-ng+1:izmax),& ! send buffer
                      nx*ng*ng,& ! send buffer size
                      MPI_DOUBLE_PRECISION,&
                      destination,& ! whom to send it to
                      mpirank,& ! send tag
                      field(:,iymax+1:iymax+1+ng-1,izmin-1-ng+1:izmin-1),&! receive buffer
                      nx*ng*ng,& ! recvcount
                      MPI_DOUBLE_PRECISION,&
                      origin,& !source
                      origin,& ! recv tag
                      MPI_COMM_WORLD,status,mpicode)                            
  endif
  
  write(str,'(i3.3)') mpirank
  do iy=ga(2),gb(2)
  open(14,file='fuck'//str,status='unknown',position='append')
  write(14,'(100(f6.2,1x))') field(0,iy,:)
  close(14)
  enddo
end subroutine synchronize_ghosts