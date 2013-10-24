

subroutine cal_drag ( time, u )
  use mpi_header
  use fsi_vars
  implicit none
  
  real(kind=pr),intent(in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(in) :: time
  
  integer :: ix,iy,iz, mpicode, ixmin,ixmax,iymin,iymax,izmin,izmax
  real(kind=pr) :: penalx,penaly,penalz,forcex,forcey,forcez,xlev,ylev,zlev
  real(kind=pr) :: torquex,torquey,torquez,usx,usy,usz

  if ((iPenalization==1).and.(compute_forces==1)) then  
    forcex  = 0.d0
    forcey  = 0.d0
    forcez  = 0.d0  
    torquex = 0.d0
    torquey = 0.d0
    torquez = 0.d0  
    usx     = 0.d0
    usy     = 0.d0
    usz     = 0.d0  
    
    !-----------------------------------------
    ! define window for integration
    !-----------------------------------------      
    if (iCavity == "no") then
      ixmin = 0
      ixmax = nx-1
      iymin = 0
      iymax = ny-1
      izmin = 0
      izmax = nz-1      
    else
      ixmin = cavity_size+2
      ixmax = nx-1-cavity_size-2
      iymin = cavity_size+2
      iymax = ny-1-cavity_size-2
      izmin = cavity_size+2
      izmax = nz-1-cavity_size-2        
    endif
    
    !-----------------------------------------
    ! loop over penalization term (this saves a work array)
    !-----------------------------------------  
    do ix=ra(1),rb(1)
      do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)   
          if ((ix>=ixmin).and.(ix<=ixmax)) then
          if ((iy>=iymin).and.(iy<=iymax)) then
          if ((iz>=izmin).and.(iz<=izmax)) then
          
          if (iMoving==1) then
            usx = us(ix,iy,iz,1)
            usy = us(ix,iy,iz,2)
            usz = us(ix,iy,iz,3)
          endif
          ! actual penalization term
          penalx = -mask(ix,iy,iz)*(u(ix,iy,iz,1)-usx)
          penaly = -mask(ix,iy,iz)*(u(ix,iy,iz,2)-usy)
          penalz = -mask(ix,iy,iz)*(u(ix,iy,iz,3)-usz)
          
          ! for torque moment
          xlev = dble(ix)*dx - x0
          ylev = dble(iy)*dy - y0
          zlev = dble(iz)*dz - z0
          
          ! integrate forces + torques (note sign inversion!)
          forcex = forcex - penalx
          forcey = forcey - penaly
          forcez = forcez - penalz              
          torquex = torquex - (ylev*penalz - zlev*penaly)
          torquey = torquey - (zlev*penalx - xlev*penalz)
          torquez = torquez - (xlev*penaly - ylev*penalx)        
          
          endif
          endif
          endif
        enddo
      enddo
    enddo  
    
    ! save global forces
    forcex = forcex*dx*dy*dz
    forcey = forcey*dx*dy*dz
    forcez = forcez*dx*dy*dz  
    call MPI_REDUCE (forcex,GlobalIntegrals%Force(1),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode)  
    call MPI_REDUCE (forcey,GlobalIntegrals%Force(2),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode) 
    call MPI_REDUCE (forcez,GlobalIntegrals%Force(3),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode)    
                    
    torquex = torquex*dx*dy*dz
    torquey = torquey*dx*dy*dz
    torquez = torquez*dx*dy*dz  
    call MPI_REDUCE (torquex,GlobalIntegrals%Torque(1),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode)  
    call MPI_REDUCE (torquey,GlobalIntegrals%Torque(2),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode) 
    call MPI_REDUCE (torquez,GlobalIntegrals%Torque(3),1,mpireal,MPI_SUM,0,&
                    MPI_COMM_WORLD,mpicode)    
                    
    !-------------------------------------------------------
    ! write time series to disk
    ! note: we also dump the unst corrections and thus suppose that they
    ! have been computed
    !-------------------------------------------------------     
    if(mpirank == 0) then
      open(14,file='forces.t',status='unknown',position='append')
      write (14,'(13(e12.5,1x))') time, GlobalIntegrals%Force, &
      GlobalIntegrals%Force_unst, GlobalIntegrals%Torque, &
      GlobalIntegrals%Torque_unst
      close(14)
    endif
    
  endif
end subroutine cal_drag



!-------------------------------------------------------------------------------
! computes the unsteady corrections, if this is set in the params file
! INPUT: 
!       time, dt: current time level and OLD time step
!       NOTE: the mask is at time t(n), so dt=t(n)-t(n-1)
! OUTPUT:
!       saves the unsteady corrections in the global force datatype
! NOTES:
!       this code uses persistent variables to determine the time derivative.
!       it will fail in the very first time step, also after retaking a backup. 
!       if it is not possible to compute the unst corrections, we return 0 here.
!-------------------------------------------------------------------------------
subroutine cal_unst_corrections ( time, dt )
  use mpi_header
  use fsi_vars
  implicit none
  real(kind=pr),intent(in) :: time, dt
  ! is it possible to compute unsteady forces?
  logical, save :: is_possible = .false.
  ! the old value of the integral
  real(kind=pr),dimension(1:3),save :: force_old = 0.d0, torque_old=0.d0
  real(kind=pr),dimension(1:3) :: force_new, force_new_loc, torque_new, torque_new_loc
  real(kind=pr) :: xlev,ylev,zlev,usx,usy,usz
  integer :: mpicode, ix, iy, iz
  
  
  ! we currently do not need to take care of the cavity since us is zero there
  
  !------------------------------------
  ! force
  !------------------------------------
  force_new_loc(1) = sum( mask*us(:,:,:,1) )*dx*dy*dz
  force_new_loc(2) = sum( mask*us(:,:,:,2) )*dx*dy*dz
  force_new_loc(3) = sum( mask*us(:,:,:,3) )*dx*dy*dz  
  
  call MPI_REDUCE(force_new_loc(1),force_new(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)  
  call MPI_REDUCE(force_new_loc(2),force_new(2),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)   
  call MPI_REDUCE(force_new_loc(3),force_new(3),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)   
  
  if (is_possible) then
    GlobalIntegrals%Force_unst(1) = force_new(1)-force_old(1)
    GlobalIntegrals%Force_unst(2) = force_new(2)-force_old(2)
    GlobalIntegrals%Force_unst(3) = force_new(3)-force_old(3)    
    GlobalIntegrals%Force_unst = GlobalIntegrals%Force_unst / dt
  else
    GlobalIntegrals%Force_unst = 0.d0    
  endif  
  force_old = force_new
  
  
  !------------------------------------
  ! torque
  !------------------------------------
  torque_new_loc = 0.d0
  
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)  
        ! for torque moment
        xlev = dble(ix)*dx - x0
        ylev = dble(iy)*dy - y0
        zlev = dble(iz)*dz - z0    
        
        usx = us(ix,iy,iz,1)
        usy = us(ix,iy,iz,2)
        usz = us(ix,iy,iz,3)
        
        torque_new_loc(1) = torque_new_loc(1) + (ylev*usz - zlev*usy)
        torque_new_loc(2) = torque_new_loc(2) + (zlev*usx - xlev*usz)
        torque_new_loc(3) = torque_new_loc(3) + (xlev*usy - ylev*usx)   
      enddo
    enddo
  enddo
  
  torque_new_loc = torque_new_loc*dx*dy*dz
  
  call MPI_REDUCE(torque_new_loc(1),torque_new(1),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)  
  call MPI_REDUCE(torque_new_loc(2),torque_new(2),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)   
  call MPI_REDUCE(torque_new_loc(3),torque_new(3),1,mpireal,MPI_SUM,0,MPI_COMM_WORLD,mpicode)     
  
  
  if (is_possible) then
    GlobalIntegrals%Torque_unst(1) = torque_new(1)-torque_old(1)
    GlobalIntegrals%Torque_unst(2) = torque_new(2)-torque_old(2)
    GlobalIntegrals%Torque_unst(3) = torque_new(3)-torque_old(3)    
    GlobalIntegrals%Torque_unst = GlobalIntegrals%Torque_unst / dt
  else
    GlobalIntegrals%Torque_unst = 0.d0    
  endif  
  torque_old = torque_new  
  
  
  ! now we sure have the old value in the next step
  is_possible = .true.
end subroutine cal_unst_corrections