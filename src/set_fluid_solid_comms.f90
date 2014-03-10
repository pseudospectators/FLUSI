subroutine set_fluid_solid_communicators
  !-----------------------------------------------------------------------------
  ! Sets ups the communicators for the solid and the fluid parts
  !-----------------------------------------------------------------------------
  use vars
  use mpi
  implicit none
  integer :: mpicode
  integer :: i
  integer, allocatable :: fluid_ranks(:), solid_ranks(:)
  integer :: original_group
  integer :: group_fluid, group_solid, mpirank_global
  
  !-- ncpu is the total number of processes
  call MPI_COMM_SIZE ( MPI_COMM_WORLD, ncpu, mpicode )
  call MPI_COMM_RANK ( MPI_COMM_WORLD, mpirank_global,mpicode ) 
  
  !-----------------------------------------------------------------------------
  !-- decide how many CPU we reserve for the solid solver
  !-----------------------------------------------------------------------------
  if ( method == "mhd" ) then
    !-- in the MHD case, we deal only with one group of CPU, no ncpu_solid
    ncpu_solid = 0
    ncpu_fluid = ncpu
  elseif ( method == "fsi" ) then
    !-- in the MHD case, we deal only with one group of CPU, no ncpu_solid
    ncpu_solid = 1
    ncpu_fluid = ncpu - ncpu_solid
  endif 
  
  
  !-----------------------------------------------------------------------------
  !-- define list of FLUID / SOLID ranks
  !-----------------------------------------------------------------------------
  allocate ( fluid_ranks(1:ncpu_fluid) )
  allocate ( solid_ranks(1:ncpu_solid) )
  
  do i = 1, ncpu_fluid
    fluid_ranks(i) = (i-1) ! note zero based indexing
    ! save globally that this is a fluid CPU, if my rank matches of these
    if ( mpirank_global == fluid_ranks(i) ) fluid_cpu = .true.
  enddo
    
  do i = 1, ncpu_solid
    solid_ranks(i) = (ncpu - (i-1)) - 1 ! note zero based indexing
    ! save globally that this is a solid CPU, if my rank matches of these
    if ( mpirank_global == solid_ranks(i) ) solid_cpu = .true.
  enddo
  
  !-----------------------------------------------------------------------------
  !-- create MPI communicators
  !-----------------------------------------------------------------------------  
  call MPI_COMM_GROUP ( MPI_COMM_WORLD, original_group, mpicode )    
  ! ----------------------------FLUID-------------------------------------------      
  call MPI_GROUP_INCL ( original_group, ncpu_fluid, fluid_ranks, group_fluid, mpicode )     
  call MPI_COMM_CREATE ( MPI_COMM_WORLD, group_fluid, MPI_COMM_FLUID, mpicode  )
  ! ----------------------------SOLID-------------------------------------------      
  call MPI_GROUP_INCL ( original_group, ncpu_solid, solid_ranks, group_solid, mpicode )     
  call MPI_COMM_CREATE ( MPI_COMM_WORLD, group_solid, MPI_COMM_SOLID, mpicode  )
  ! ----------------------------------------------------------------------------  
  
end subroutine