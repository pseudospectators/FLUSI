subroutine get_surface_pressure_jump (time, beam, p)
  use mpi
  use fsi_vars
  use interpolation
  implicit none
  real(kind=pr),intent (in) :: time
  type(solid),intent (inout) :: beam
  real(kind=pr), intent (in)   :: p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: xf,yf,zf,dh, soft_startup, t0
  real(kind=pr) :: psi,gamma,tmp,tmp2,psi_dt,beta_dt,gamma_dt,beta
  real(kind=pr),dimension(:,:,:,:), allocatable :: surfaces
  real(kind=pr),dimension(:,:,:), allocatable :: p_surface, p_surface_local
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate
  real(kind=pr),dimension(1:3) :: u_tmp,rot_body,v_tmp,v0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate
  real(kind=pr),dimension(:,:),allocatable::ghostsy,ghostsz
  integer :: nh,is,ih,isurf,mpicode
  t0 = MPI_wtime()
  !-- get relative coordinate system
  call plate_coordinate_system( time,x0_plate,v0_plate,psi,beta,gamma,psi_dt,beta_dt,gamma_dt,M_plate)
  
  !-- number of interpolation points in rigid direction (span)
  nh = nint( L_span/min(dx,dy,dz)  )
  !-- spacing in span direction
  dh = L_span/dble(nh)
  
  ! this array holds the interpolation points (2 2D arrays of 3D vectors = 4 indices)
  allocate(surfaces(0:ns-1,0:nh,1:2,1:3))
  allocate(p_surface(0:ns-1,0:nh,1:2))
  allocate(p_surface_local(0:ns-1,0:nh,1:2))

  surfaces = 17.d0
  p_surface = 17.d0
  p_surface_local = 17.d0
  
  !-----------------------------------------------------------------------------
  ! Compute the two 2D surfaces in 3D space (top and bottom)
  !-----------------------------------------------------------------------------
  do is=0,ns-1
    do ih=0,nh
      !-- top surface points
      xf = beam%x(is)-t_beam*dsin(beam%theta(is))
      yf = beam%y(is)+t_beam*dcos(beam%theta(is))    
      x_plate = (/ xf,yf,dble(ih)*dh-0.5d0*L_span /)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,1,1:3) = x
      
      !-- bottom surface points
      xf = beam%x(is)+t_beam*dsin(beam%theta(is))
      yf = beam%y(is)-t_beam*dcos(beam%theta(is))      
      x_plate = (/ xf,yf,dble(ih)*dh-0.5d0*L_span /)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,2,1:3) = x    
      
    enddo
  enddo
  
  !-- this just copies ra(1:3) and rb(1:3) in more readble names
  call init_interpolation()
  
  !-----------------------------------------------------------------------------
  ! Tri-linear interpolation of all points on the surface. 
  ! If the point does not lie in the locally stores memory, zero is set
  ! Then, we can sum over all MPI processes and have the complete pressure on
  ! the surface on all ranks
  !-----------------------------------------------------------------------------
  if ((mpidims(2)==1).and.(mpisize>1)) then
    !------------------------
    !-- 1D data decompositon
    !------------------------
    allocate (ghostsz(ixmin:ixmax,iymin:iymax) )
    !--fetch a sheet from neighbor
    call extend_array_1D( p, ghostsz )
    
    do is=0,ns-1
      do ih=0,nh
        do isurf=1,2
          x = surfaces(is,ih,isurf,1:3)
          call trilinear_interp_1Ddecomp( x, p, ghostsz, p_surface_local(is,ih,isurf))
        enddo
      enddo
    enddo
    
    deallocate (ghostsz)    

  elseif (mpisize==1) then
    !------------------------
    !-- serial run one CPU
    !------------------------
    do is=0,ns-1
      do ih=0,nh
        do isurf=1,2
          x = surfaces(is,ih,isurf,1:3)
          call trilinear_interp( x, p, p_surface_local(is,ih,isurf))
        enddo
      enddo
    enddo
    
  else
    !------------------------
    !-- 2D data decomposition
    !------------------------
    allocate(ghostsz(ixmin:ixmax,iymin:iymax+1))
    allocate(ghostsy(ixmin:ixmax,izmin:izmax+1))
    !-- fetch two ghost surfaces from other CPU
    call extend_array_2D( p, ghostsz, ghostsy)
    
    do is=0,ns-1
      do ih=0,nh
        do isurf=1,2
          x = surfaces(is,ih,isurf,1:3)
          call trilinear_interp_2Ddecomp( x, p, ghostsz,ghostsy, p_surface_local(is,ih,isurf))
        enddo
      enddo
    enddo
    
    deallocate(ghostsz)
    deallocate(ghostsy)
  endif
  
  
  !-----------------------------------------------------------------------------
  ! All CPUs now interpolated some values of p on the surfaces, if they lie in
  ! their local memory (or on the lines between two CPUs). If not, they saved 
  ! zero. The sum of all pressure distributions is now what we wanted to have
  ! p_surface is available on all CPU (and it has to since all evolve the solid)
  !-----------------------------------------------------------------------------
  call MPI_ALLREDUCE ( p_surface_local,p_surface,size(p_surface_local),&
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode) 
  
  
  
  !-----------------------------------------------------------------------------
  ! To avoid startup problems, we can smoothly "turn on" the coupling by 
  ! multiplying the pressure with a startup conditioner.
  ! Note this does not apply for gravity or imposed motion.
  !-----------------------------------------------------------------------------  
  if (time <= T_release) then
      soft_startup = 0.0
  elseif ( ( time >T_release ).and.(time<(T_release + tau)) ) then
      soft_startup =  ((time-T_release)**3)/(-0.5*tau**3)   + 3.*((time-T_release)**2)/tau**2
  else
      soft_startup = 1.0
  endif
  
  
  !-----------------------------------------------------------------------------
  ! Average the pressure in the spanwise direction and multiply with startup 
  ! conditioner.
  !-----------------------------------------------------------------------------    
  do is=0,ns-1 
    beam%pressure_old(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)
    beam%pressure_new(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)
    
    beam%pressure_old(is) = beam%pressure_old(is)*soft_startup
    beam%pressure_new(is) = beam%pressure_new(is)*soft_startup
  enddo
      
!   if (mpirank==0) then
!   open (17, file='surfaces.m', status='replace')
!   
!   write (17,'(A)') 'x1=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,1,1)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'y1=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,1,2)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'z1=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,1,3)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'p1=['
!   do is=0,ns-1
!   write(17,'(256(f10.2,1x))') p_surface(is,:,1)
!   enddo  
!   write (17,'(A)') '];'
!   
!   
!   
!   write (17,'(A)') 'x2=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,2,1)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'y2=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,2,2)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'z2=['
!   do is=0,ns-1
!   write(17,'(256(f5.2,1x))') surfaces(is,:,2,3)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') 'p2=['
!   do is=0,ns-1
!   write(17,'(256(f10.2,1x))') p_surface(is,:,2)
!   enddo  
!   write (17,'(A)') '];'
!   
!   write (17,'(A)') "surface(x1,y1,z1,p1); hold on; surface(x2,y2,z2,p2)"
!   write (17,'(A)') "axis([0 4 0 4 0 4]);xlabel('x');ylabel('y');zlabel('z');grid"
!   close(17)
!     endif
    
    
  deallocate(surfaces)
  deallocate(p_surface)
  deallocate(p_surface_local)
  
  time_surf = time_surf + MPI_wtime() - t0
end subroutine get_surface_pressure_jump