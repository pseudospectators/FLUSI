subroutine get_surface_pressure_jump (time, beam, p)
  use mpi
  use fsi_vars
  use interpolation
  implicit none
  real(kind=pr),intent (in) :: time
  type(solid),intent (inout) :: beam
  real(kind=pr), intent (in)   :: p(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: xf,yf,zf,dh
  real(kind=pr) :: psi,gamma,tmp,tmp2,psi_dt,beta_dt,gamma_dt,beta
  real(kind=pr),dimension(:,:,:,:), allocatable :: surfaces
  real(kind=pr),dimension(:,:,:), allocatable :: p_surface, p_surface_local
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate
  real(kind=pr),dimension(1:3) :: u_tmp,rot_body,v_tmp,v0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate
  real(kind=pr),dimension(:,:),allocatable::ghostsy,ghostsz
  integer :: nh,is,ih,isurf,mpicode
  
  !-- get relative coordinate system
  call plate_coordinate_system( time,x0_plate,v0_plate,psi,beta,gamma,psi_dt,beta_dt,gamma_dt,M_plate)
  
  ! number of interpolation points in rigid direction (span)
  nh = nint( L_span/min(dx,dy,dz)  )
  ! spacing in span direction
  dh = L_span/dble(nh)
  
  ! this array holds the interpolation points (2 2D arrays of 3D vectors = 4 indices)
  allocate(surfaces(0:ns-1,0:nh,1:2,1:3))
  allocate(p_surface(0:ns-1,0:nh,1:2))
  allocate(p_surface_local(0:ns-1,0:nh,1:2))

  !-----------------------------------------------------------------------------
  ! Compute the two 2D surfaces in 3D space (top and bottom)
  !-----------------------------------------------------------------------------
  do is=0,ns-1
    do ih=0,nh
      !-- top points
      xf = beam%x(is)-t_beam*dsin(beam%theta(is))
      yf = beam%y(is)+t_beam*dcos(beam%theta(is))
    
      x_plate = (/ xf,yf,dble(ih)*dh-0.5d0*L_span /)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,1,1:3) = x
      
      !-- bottom points
      xf = beam%x(is)+t_beam*dsin(beam%theta(is))
      yf = beam%y(is)-t_beam*dcos(beam%theta(is))
      
      x_plate = (/ xf,yf,dble(ih)*dh-0.5d0*L_span /)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,2,1:3) = x    
      
    enddo
  enddo
  
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
    do is=0,ns-1
    write(*,*) '!-- serial run one CPU'
    enddo
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
  
  
  call MPI_ALLREDUCE ( p_surface_local,p_surface,ns*nh*2,MPI_DOUBLE_PRECISION,&
                       MPI_SUM,MPI_COMM_WORLD,mpicode) 
  
  if (mpirank==0) then
  open (14, file='surfaces.dat', status='replace')
  
  write (14,'(A)') 'x1=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,1,1)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'y1=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,1,2)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'z1=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,1,3)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'p1=['
  do is=0,ns-1
  write(14,'(256(f10.2,1x))') p_surface(is,:,1)
  enddo  
  write (14,'(A)') '];'
  
  
  
  write (14,'(A)') 'x2=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,2,1)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'y2=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,2,2)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'z2=['
  do is=0,ns-1
  write(14,'(256(f5.2,1x))') surfaces(is,:,2,3)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') 'p2=['
  do is=0,ns-1
  write(14,'(256(f10.2,1x))') p_surface(is,:,2)
  enddo  
  write (14,'(A)') '];'
  
  write (14,'(A)') "surface(x1,y1,z1,p1); hold on; surface(x2,y2,z2,p2)"
  write (14,'(A)') "axis([0 4 0 4 0 4]);xlabel('x');ylabel('y');zlabel('z');grid"
    endif
  deallocate(surfaces, p_surface, p_surface_local)
!   stop
end subroutine get_surface_pressure_jump