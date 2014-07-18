subroutine get_surface_pressure_jump (time, beam, p, testing, timelevel)
  use mpi
  use fsi_vars
  use interpolation
  implicit none
  real(kind=pr),intent (in) :: time
  type(solid),intent (inout) :: beam
  character(len=*), intent(in), optional :: testing, timelevel
  real(kind=pr), intent (in)   :: p(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr) :: xf,yf,zf,dh, soft_startup, t0
  real(kind=pr) :: psi,gamma,tmp,tmp2,psi_dt,beta_dt,gamma_dt,beta
  real(kind=pr),dimension(:,:,:,:), allocatable :: surfaces
  real(kind=pr),dimension(:,:,:), allocatable :: p_surface, p_surface_local
  real(kind=pr),dimension(1:3) :: x, x_plate, x0_plate
  real(kind=pr),dimension(1:3) :: u_tmp,rot_body,v_tmp,v0_plate
  real(kind=pr),dimension(1:3,1:3) :: M_plate
  real(kind=pr),dimension(:,:),allocatable::ghostsy,ghostsz
  real(kind=pr),dimension(:),allocatable::heights
  integer :: nh,is,ih,isurf,mpicode
  t0 = MPI_wtime()
  
  !-- get relative coordinate system
  call plate_coordinate_system( time,x0_plate,v0_plate,psi,beta,gamma,psi_dt,beta_dt,gamma_dt,M_plate)
  
  
  if (nx>1) then
    !-- "true" 3D case
    !-- number of interpolation points in rigid direction (span)
    nh = nint( L_span/min(dx,dy,dz)  )
    !-- spacing in span direction
    dh = L_span/dble(nh)
    allocate(heights(0:nh))
    do ih=0,nh
      heights(ih) = dble(ih)*dh -0.5d0*L_span
    enddo
    
  elseif (nx==1) then
  
    nh = 0
    dh = 1.d0
    allocate(heights(0:nh))
    !-- in 2D case, just interpolate always the same position (subsequent avg
    !-- leaves value untouched)
    heights = 0.d0
  endif
  
  
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
      x_plate = (/ xf,yf,heights(ih)/)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,1,1:3) = x

      
      !-- bottom surface points
      xf = beam%x(is)+t_beam*dsin(beam%theta(is))
      yf = beam%y(is)-t_beam*dcos(beam%theta(is))      
      x_plate = (/ xf,yf,heights(ih)/)
      x = matmul( transpose(M_plate) , x_plate )
      x = x + x0_plate
      surfaces(is,ih,2,1:3) = x    
      
      if (root) then
      if((x(1)<-dx).or.(x(1)>xl+dx).or.(x(2)<-dy).or.(x(2)>yl+dy).or.(x(3)<-dz).or.(x(3)>zl+dz)) then
        write(*,*) "ERROR: Surface constuction yielded values outside valid bound."
        call abort()
      endif
      endif
    enddo
  enddo
  
  deallocate( heights )
    
  !-----------------------------------------------------------------------------
  ! Tri-linear interpolation of all points on the surface. 
  ! If the point does not lie in the locally stores memory, zero is set
  ! Then, we can sum over all MPI processes and have the complete pressure on
  ! the surface on all ranks
  !-----------------------------------------------------------------------------
  call synchronize_ghosts ( p )
  
  do is=0,ns-1
    do ih=0,nh
      do isurf=1,2
        x = surfaces(is,ih,isurf,1:3)
        if (interp=='linear') then
          call trilinear_interp_ghosts( x, p, p_surface_local(is,ih,isurf))
        elseif (interp=='delta') then
          call delta_interpolation( x, p, p_surface_local(is,ih,isurf))
        endif
      enddo
    enddo
  enddo
  
  !-----------------------------------------------------------------------------
  ! All CPUs now interpolated some values of p on the surfaces, if they lie in
  ! their local memory (or on the lines between two CPUs). If not, they saved 
  ! -9.9e10. The max of all pressure distributions is now what we wanted to have
  ! p_surface is available on all CPU (and it has to since all evolve the solid)
  !-----------------------------------------------------------------------------
  call MPI_ALLREDUCE ( p_surface_local,p_surface,size(p_surface_local),&
                       MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpicode) 
  
  !-----------------------------------------------------------------------------
  ! Testing. Interpolate the mask function (given by p from external) at the beam
  ! surfaces. the value should be in the range 0.05 and 0.60 (roughly). ideally,
  ! the value should be 0.5
  !-----------------------------------------------------------------------------
  if (present(testing)) then
  if (((testing=="mask").or.(testing=="linear")).and.(root)) then
    write(*,'(80("-"))')
    write(*,*) testing//" Interpolation test, top surface:"
    do ih=0,nh
      do is=0,ns-1      
        if (testing=="mask") then
          write(*,'((f5.3,1x))',advance='no') p_surface(is,ih,1)
        elseif (testing=="linear") then
          write(*,'((f5.3,1x))',advance='no') p_surface(is,ih,1)&
          -sin(surfaces(is,ih,1,2))*cos(surfaces(is,ih,1,3))
        endif
      enddo
      write(*,*) " "
    enddo
    write(*,*) testing//" Interpolation test, bottom surface:"
    do ih=0,nh
      do is=0,ns-1
        if (testing=="mask") then
          write(*,'((f5.3,1x))',advance='no') p_surface(is,ih,2)
        elseif (testing=="linear") then
          write(*,'((f5.3,1x))',advance='no') p_surface(is,ih,2)&
          -sin(surfaces(is,ih,2,2))*cos(surfaces(is,ih,2,3))
        endif
      enddo
      write(*,*) " "
    enddo
    write(*,'(80("-"))  ')
  endif
  endif
  
  !-----------------------------------------------------------------------------
  ! To avoid startup problems, we can smoothly "turn on" the coupling by 
  ! multiplying the pressure with a startup conditioner.
  ! Note this does not apply for gravity or imposed motion.
  !-----------------------------------------------------------------------------  
  if (time <= T_release) then
      soft_startup = 0.0
  elseif ( ( time >T_release ).and.(time<(T_release + tau)) ) then
      soft_startup =  ((time-T_release)**3)/(-0.5*tau**3) + 3.*((time-T_release)**2)/tau**2
  else
      soft_startup = 1.0
  endif
  
  
  !-----------------------------------------------------------------------------
  ! Average the pressure in the spanwise direction and multiply with startup 
  ! conditioner. 
  !-----------------------------------------------------------------------------    
  if (.not.(present(testing))) then
    do is=0,ns-1 
      if (present(timelevel)) then
        if (timelevel=="old") then
          !-- store result of interpolation at OLD timelevel
          beam%pressure_old(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)
          beam%pressure_old(is) = beam%pressure_old(is)*soft_startup
        elseif (timelevel=="new") then
          !-- store result of interpolation at NEW timelevel
          beam%pressure_new(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)        
          beam%pressure_new(is) = beam%pressure_new(is)*soft_startup
        endif
      else
        ! write into both
        beam%pressure_old(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)
        beam%pressure_new(is) = sum(p_surface(is,:,1)-p_surface(is,:,2)) / dble(nh+1)
        
        beam%pressure_old(is) = beam%pressure_old(is)*soft_startup
        beam%pressure_new(is) = beam%pressure_new(is)*soft_startup
      endif      
    enddo
  endif
      
  deallocate(surfaces)
  deallocate(p_surface)
  deallocate(p_surface_local)
  
  time_surf = time_surf + MPI_wtime() - t0
end subroutine get_surface_pressure_jump


!-------------------------------------------------------------------------------
! Testing routine for the surface interpolation, called in every run. First, we
! interpolate the mask function at the beams's surface - the value must be close 
! to 0.5. If not, then mask function and solid model are not consistently defned.
! Second, we interpolate and analytic function f=(y*z) to compare with the exact
! value. the error should be zero
!
! subroutine GET_SURFACE_PRESSURE_JUMP has an optional argument for this testing
!
!-------------------------------------------------------------------------------
subroutine surface_interpolation_testing( time, beams, work )
  use mpi
  use fsi_vars
  use interpolation
  use insect_module
  implicit none
  real(kind=pr),intent (in) :: time
  type(solid),dimension(1),intent (inout) :: beams
  real(kind=pr),intent(inout):: work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer :: ix,iy,iz
  real(kind=pr) :: x,y,z
  type(diptera)::insectdummy
  
  !-- create mask
!   call Draw_flexible_plate(time, beam)
  call create_mask(time,insectdummy,beams)
  !-- copy mask in extended array (the one with ghost points)
  work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) = mask*eps
  !-- interpolate the mask at the beam surfaces
  call get_surface_pressure_jump (time, beams(1), work, testing="mask")

  work=0.d0
  !-- the array "mask" with some values (analytic)
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)
        x=dble(ix)*dx
        y=dble(iy)*dy
        z=dble(iz)*dz
        work(ix,iy,iz)=sin(y)*cos(z)
      enddo
    enddo
  enddo 
  !-- interpolate these values on the surface and compare to exact solution
  call get_surface_pressure_jump (time, beams(1), work, testing="linear")
  
  
end subroutine
