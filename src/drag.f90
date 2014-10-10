!-------------------------------------------------------------------------------
! Compute hydrodynamic forces and torque
!-------------------------------------------------------------------------------
! The forces are the volume integral of the penalization term. Different parts of
! the mask are colored differently.
! Color         Description
!   0           Boring parts (channel walls, cavity)
!   1           Interesting parts (e.g. a cylinder), for the insects this is BODY
!   2           Other parts, for the insects, this is LEFT WING
!   3           For the insects, this is RIGHT WING
! Currently, we store the torque / forces over all colors greater than 0 in the
! global structure "GlobalIntegrals". If we're running in "insects" mode, the 
! colors 1 and 2 are the forces on wings and body, respectively. These are stored
! in the "Insect" global struct.
!-------------------------------------------------------------------------------
subroutine cal_drag ( t, u, mask, mask_color, us, Insect )
  use vars
  use insect_module
  implicit none
  
  type(timetype),intent(in) :: t
  real(kind=pr),intent(in) :: u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(in) :: mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(in) :: us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(in)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(diptera), intent(inout) :: Insect
  
  integer :: ix,iy,iz,mpicode
  integer(kind=2) :: color
  real(kind=pr) :: penalx,penaly,penalz,xlev,ylev,zlev,apowtotal,ipowtotal
  ! we can choose up to 6 different colors
  real(kind=pr),dimension(0:5) :: torquex,torquey,torquez,forcex,forcey,forcez
  real(kind=pr) :: torquex0,torquey0,torquez0
  ! power (flux of energy) 
  real(kind=pr) :: power,powerx,powery,powerz, time
  character(len=1024) :: forcepartfilename
  
  time = t%time
  
  forcex  = 0.d0
  forcey  = 0.d0
  forcez  = 0.d0  
  torquex = 0.d0
  torquey = 0.d0
  torquez = 0.d0
  torquex0 = 0.d0
  torquey0 = 0.d0
  torquez0 = 0.d0
  powerx = 0.d0
  powery = 0.d0
  powerz = 0.d0
  
  !---------------------------------------------------------------------------
  ! loop over penalization term (this saves a work array)
  !---------------------------------------------------------------------------
  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        ! actual penalization term
        penalx = -mask(ix,iy,iz)*(u(ix,iy,iz,1)-us(ix,iy,iz,1))
        penaly = -mask(ix,iy,iz)*(u(ix,iy,iz,2)-us(ix,iy,iz,2))
        penalz = -mask(ix,iy,iz)*(u(ix,iy,iz,3)-us(ix,iy,iz,3))
        
        ! what color does the point have?
        color = mask_color(ix,iy,iz)
        
        ! integrate forces + torques (note sign inversion!)
        forcex(color) = forcex(color) - penalx
        forcey(color) = forcey(color) - penaly
        forcez(color) = forcez(color) - penalz

        ! for torque moment with respect to (x0,y0,z0)
        xlev = dble(ix)*dx - x0
        ylev = dble(iy)*dy - y0
        zlev = dble(iz)*dz - z0
 
        ! compute moment with respect to (x0,y0,z0)
        torquex0 = torquex0 - (ylev*penalz - zlev*penaly)
        torquey0 = torquey0 - (zlev*penalx - xlev*penalz)
        torquez0 = torquez0 - (xlev*penaly - ylev*penalx)
        
        ! input power due to penalization term
        powerx = powerx + u(ix,iy,iz,1)*penalx
        powery = powery + u(ix,iy,iz,2)*penaly
        powerz = powerz + u(ix,iy,iz,3)*penalz        

        ! for insects, moment of the body is computed with respect to (x0,y0,z0)
        ! but for the wings it is computed with respect tot the pivot points
        if (iMask=="Insect") then
          if (color==2) then
            xlev = xlev - Insect%x_pivot_l_glob(1)
            ylev = ylev - Insect%x_pivot_l_glob(2)
            zlev = zlev - Insect%x_pivot_l_glob(3)
          elseif (color==3) then
            xlev = xlev - Insect%x_pivot_r_glob(1)
            ylev = ylev - Insect%x_pivot_r_glob(2)
            zlev = zlev - Insect%x_pivot_r_glob(3)
          endif
        endif

        ! Compute moments relative to each part
        torquex(color) = torquex(color) - (ylev*penalz - zlev*penaly)
        torquey(color) = torquey(color) - (zlev*penalx - xlev*penalz)
        torquez(color) = torquez(color) - (xlev*penaly - ylev*penalx)
      enddo
    enddo
  enddo  
  
  !---------------------------------------------------------------------------
  ! save global forces
  !---------------------------------------------------------------------------
  forcex = forcex*dx*dy*dz
  forcey = forcey*dx*dy*dz
  forcez = forcez*dx*dy*dz  
  
  powerx = powerx*dx*dy*dz
  powery = powery*dx*dy*dz
  powerz = powerz*dx*dy*dz
  power  = powerx+powery+powerz
  
  ! in the global structure, we store all contributions with color > 0, so we
  ! only EXCLUDE channel / cavity walls (the boring stuff)
  call MPI_ALLREDUCE ( sum(forcex(1:5)),GlobalIntegrals%Force(1),1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE ( sum(forcey(1:5)),GlobalIntegrals%Force(2),1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE ( sum(forcez(1:5)),GlobalIntegrals%Force(3),1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode) 
                  
  ! penalization power ( -u.[chi*(u-us)] )
  call MPI_ALLREDUCE ( power,GlobalIntegrals%penalization_power,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( powerx,GlobalIntegrals%penalization_power_x,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( powery,GlobalIntegrals%penalization_power_y,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
  call MPI_ALLREDUCE ( powerz,GlobalIntegrals%penalization_power_z,1,&
                  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)                  
                  
                  
  ! the insects have forces on the wing and body separate
  if (iMask=="Insect") then
    do color = 1,3
      call MPI_ALLREDUCE (forcex(color),Insect%PartIntegrals(color)%Force(1),1,&
                      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
      call MPI_ALLREDUCE (forcey(color),Insect%PartIntegrals(color)%Force(2),1,&
                      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode) 
      call MPI_ALLREDUCE (forcez(color),Insect%PartIntegrals(color)%Force(3),1,&
                      MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
    enddo
  endif
  
  ! in the global structure, we store all contributions with color > 0, so we
  ! only EXCLUDE channel / cavity walls (the boring stuff)                  
  torquex0 = torquex0*dx*dy*dz
  torquey0 = torquey0*dx*dy*dz
  torquez0 = torquez0*dx*dy*dz  
  torquex = torquex*dx*dy*dz
  torquey = torquey*dx*dy*dz
  torquez = torquez*dx*dy*dz  
 
  ! Integral moment with respect to (x0,y0,z0) 
  call MPI_ALLREDUCE ( torquex0,GlobalIntegrals%Torque(1),1,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE ( torquey0,GlobalIntegrals%Torque(2),1,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode) 
  call MPI_ALLREDUCE ( torquez0,GlobalIntegrals%Torque(3),1,&
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
          
  ! the insects have torques on the wing and body separate
  if (iMask=="Insect") then
    do color = 1,3
      call MPI_ALLREDUCE (torquex(color),Insect%PartIntegrals(color)%Torque(1),1,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
      call MPI_ALLREDUCE (torquey(color),Insect%PartIntegrals(color)%Torque(2),1,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode) 
      call MPI_ALLREDUCE (torquez(color),Insect%PartIntegrals(color)%Torque(3),1,&
              MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)
    enddo
  endif
 
  ! compute aerodynamic power
  if (iMask=="Insect") then
    call aero_power (Insect,apowtotal)
    call inert_power(Insect,ipowtotal)
  endif
  
  !---------------------------------------------------------------------------
  ! write time series to disk
  ! note: we also dump the unst corrections and thus suppose that they
  ! have been computed ( call cal_unst_corrections first! )
  !---------------------------------------------------------------------------
  if(mpirank == 0) then
    if (iMask=="Insect") then
      ! Aerodynamic power is only computed for insects
      open(14,file='forces.t',status='unknown',position='append')
      write (14,'(15(es15.8,1x))') time, GlobalIntegrals%Force, &
        GlobalIntegrals%Force_unst, GlobalIntegrals%Torque, &
        GlobalIntegrals%Torque_unst, apowtotal, ipowtotal
      close(14)
      ! currently, only insects have different colors
      do color=1,3
        write (forcepartfilename, "(A11,I1,A2)") "forces_part", color, ".t"
        open(14,file=trim(forcepartfilename),status='unknown',position='append')
        write (14,'(15(es15.8,1x))') time, Insect%PartIntegrals(color)%Force, &
          Insect%PartIntegrals(color)%Force_unst, Insect%PartIntegrals(color)%Torque, &
          Insect%PartIntegrals(color)%Torque_unst, Insect%PartIntegrals(color)%APow,&
          Insect%PartIntegrals(color)%IPow
        close(14)
      enddo
    else    
      ! Not an insect, we have only one part of the force
      open(14,file='forces.t',status='unknown',position='append')
      write (14,'(13(es15.8,1x))') time, GlobalIntegrals%Force, &
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
subroutine cal_unst_corrections ( t, mask, mask_color, us, Insect )
  use vars
  use insect_module
  implicit none
  
  type(timetype),intent(in) :: t  
  type(diptera), intent(inout) :: Insect
  real(kind=pr),intent(in)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(in)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(in)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  
  ! is it possible to compute unsteady forces?
  logical, save :: is_possible = .false.
  ! the old value of the integral, component by component, for each color
  real(kind=pr),dimension(0:5),save :: force_oldx = 0.d0
  real(kind=pr),dimension(0:5),save :: force_oldy = 0.d0
  real(kind=pr),dimension(0:5),save :: force_oldz = 0.d0
  real(kind=pr),dimension(0:5),save :: torque_oldx = 0.d0
  real(kind=pr),dimension(0:5),save :: torque_oldy = 0.d0
  real(kind=pr),dimension(0:5),save :: torque_oldz = 0.d0
  real(kind=pr),save :: torque_oldx0 = 0.d0
  real(kind=pr),save :: torque_oldy0 = 0.d0
  real(kind=pr),save :: torque_oldz0 = 0.d0
  real(kind=pr) :: xlev,ylev,zlev,norm,usx,usy,usz
  ! we have up to 6 different colors
  real(kind=pr),dimension(0:5) :: force_new_locx, force_new_locy, force_new_locz
  real(kind=pr),dimension(0:5) :: torque_new_locx, torque_new_locy, torque_new_locz
  real(kind=pr) :: torque_new_locx0, torque_new_locy0, torque_new_locz0
  real(kind=pr),dimension(0:5) :: force_newx, force_newy, force_newz
  real(kind=pr),dimension(0:5) :: torque_newx, torque_newy, torque_newz
  real(kind=pr) :: torque_newx0, torque_newy0, torque_newz0
  
  integer :: mpicode, ix, iy, iz
  integer(kind=2) :: color
  real(kind=pr) :: dt, time
  
  time = t%time
  dt = t%dt_new
  
  !-----------------------------------------------------------------------------
  ! force
  !-----------------------------------------------------------------------------
  force_new_locx = 0.0
  force_new_locy = 0.0
  force_new_locz = 0.0
  
  norm = dx*dy*dz*eps
  
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)
        color = mask_color(ix,iy,iz)
        ! sum up new integral as a function of color
        force_new_locx(color) = force_new_locx(color)+mask(ix,iy,iz)*us(ix,iy,iz,1)*norm
        force_new_locy(color) = force_new_locy(color)+mask(ix,iy,iz)*us(ix,iy,iz,2)*norm
        force_new_locz(color) = force_new_locz(color)+mask(ix,iy,iz)*us(ix,iy,iz,3)*norm
      enddo
    enddo
  enddo  
    
  call MPI_ALLREDUCE(force_new_locx,force_newx,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE(force_new_locy,force_newy,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)   
  call MPI_ALLREDUCE(force_new_locz,force_newz,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)   
  
  if (is_possible) then
    ! here, we store the unsteady corrections for all nonzero colors in the global struct
    GlobalIntegrals%Force_unst(1) = sum(force_newx(1:5))-sum(force_oldx(1:5))
    GlobalIntegrals%Force_unst(2) = sum(force_newy(1:5))-sum(force_oldy(1:5))
    GlobalIntegrals%Force_unst(3) = sum(force_newz(1:5))-sum(force_oldz(1:5))    
    GlobalIntegrals%Force_unst = GlobalIntegrals%Force_unst / dt
    
    if (iMask=="Insect") then
      ! for the insects, we save separately the WINGs and the BODY
      do color = 1,3
        Insect%PartIntegrals(color)%Force_unst(1) = force_newx(color)-force_oldx(color)
        Insect%PartIntegrals(color)%Force_unst(2) = force_newy(color)-force_oldy(color)
        Insect%PartIntegrals(color)%Force_unst(3) = force_newz(color)-force_oldz(color)
        Insect%PartIntegrals(color)%Force_unst = Insect%PartIntegrals(color)%Force_unst / dt
      enddo
    endif
  else
    ! we cannot compute the time derivative, because we lack the old value of the
    ! integral. As a hack, return zero.
    GlobalIntegrals%Force_unst = 0.d0    
    if (iMask=="Insect") then
      do color = 1,3
        Insect%PartIntegrals(color)%Force_unst = 0.d0
      enddo
    endif
  endif  
  
  ! iterate
  force_oldx = force_newx
  force_oldy = force_newy
  force_oldz = force_newz
    
  !-----------------------------------------------------------------------------
  ! torque
  !-----------------------------------------------------------------------------
  torque_new_locx = 0.d0
  torque_new_locy = 0.d0
  torque_new_locz = 0.d0
  torque_new_locx0 = 0.d0
  torque_new_locy0 = 0.d0
  torque_new_locz0 = 0.d0
  
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)  
        ! for torque moment
        xlev = dble(ix)*dx - x0
        ylev = dble(iy)*dy - y0
        zlev = dble(iz)*dz - z0    
        
        usx = us(ix,iy,iz,1)*mask(ix,iy,iz)*eps
        usy = us(ix,iy,iz,2)*mask(ix,iy,iz)*eps
        usz = us(ix,iy,iz,3)*mask(ix,iy,iz)*eps

        ! moment with respect to (x0,y0,z0)
        torque_new_locx0 = torque_new_locx0 + (ylev*usz - zlev*usy)
        torque_new_locy0 = torque_new_locy0 + (zlev*usx - xlev*usz)
        torque_new_locz0 = torque_new_locz0 + (xlev*usy - ylev*usx)   

        color = mask_color(ix,iy,iz)

        ! for insects, moment of the body is computed with respect to (x0,y0,z0)
        ! but for the wings it is computed with respect tot the pivot points
        if (iMask=="Insect") then
          if (color==2) then
            xlev = xlev - Insect%x_pivot_l_glob(1)
            ylev = ylev - Insect%x_pivot_l_glob(2)
            zlev = zlev - Insect%x_pivot_l_glob(3)
          elseif (color==3) then
            xlev = xlev - Insect%x_pivot_r_glob(1)
            ylev = ylev - Insect%x_pivot_r_glob(2)
            zlev = zlev - Insect%x_pivot_r_glob(3)
          endif
        endif
        
        torque_new_locx(color) = torque_new_locx(color) + (ylev*usz - zlev*usy)
        torque_new_locy(color) = torque_new_locy(color) + (zlev*usx - xlev*usz)
        torque_new_locz(color) = torque_new_locz(color) + (xlev*usy - ylev*usx)   
      enddo
    enddo
  enddo
  
  torque_new_locx = torque_new_locx*dx*dy*dz
  torque_new_locy = torque_new_locy*dx*dy*dz
  torque_new_locz = torque_new_locz*dx*dy*dz
  torque_new_locx0 = torque_new_locx0*dx*dy*dz
  torque_new_locy0 = torque_new_locy0*dx*dy*dz
  torque_new_locz0 = torque_new_locz0*dx*dy*dz
 
  call MPI_ALLREDUCE(torque_new_locx,torque_newx,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE(torque_new_locy,torque_newy,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)   
  call MPI_ALLREDUCE(torque_new_locz,torque_newz,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)     
  call MPI_ALLREDUCE(torque_new_locx0,torque_newx0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)  
  call MPI_ALLREDUCE(torque_new_locy0,torque_newy0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)   
  call MPI_ALLREDUCE(torque_new_locz0,torque_newz0,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpicode)     
  
  
  if (is_possible) then
    ! here, we store the unsteady corrections for all nonzero colors in the global struct
    GlobalIntegrals%Torque_unst(1) = torque_newx0-torque_oldx0
    GlobalIntegrals%Torque_unst(2) = torque_newy0-torque_oldy0
    GlobalIntegrals%Torque_unst(3) = torque_newz0-torque_oldz0
    GlobalIntegrals%Torque_unst = GlobalIntegrals%Torque_unst / dt
    if (iMask=="Insect") then
      ! for the insects, we save separately the WINGs and the BODY
      do color = 1,3
        Insect%PartIntegrals(color)%Torque_unst(1) = torque_newx(color)-torque_oldx(color)
        Insect%PartIntegrals(color)%Torque_unst(2) = torque_newy(color)-torque_oldy(color)
        Insect%PartIntegrals(color)%Torque_unst(3) = torque_newz(color)-torque_oldz(color)
        Insect%PartIntegrals(color)%Torque_unst = Insect%PartIntegrals(color)%Torque_unst / dt
      enddo
    endif    
  else
    ! we cannot compute the time derivative, because we lack the old value of the
    ! integral. As a hack, return zero.
    GlobalIntegrals%Torque_unst = 0.d0   
    if (iMask=="Insect") then
      do color = 1,3
        Insect%PartIntegrals(color)%Torque_unst = 0.d0
      enddo
    endif
  endif  
  
  ! iterate
  torque_oldx = torque_newx
  torque_oldy = torque_newy
  torque_oldz = torque_newz  
  torque_oldx0 = torque_newx0
  torque_oldy0 = torque_newy0
  torque_oldz0 = torque_newz0
  
  ! now we sure have the old value in the next step
  is_possible = .true.
end subroutine cal_unst_corrections
