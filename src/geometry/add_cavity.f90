!-------------------------------------------------------------------------------
! add a cavity mask around the domain
! note: if the obstacle extends in this zone (e.g. the flying insect touching
! the lower end of the domain), it will be overwritten.
! the cavity walls have the color "0" which helps excluding them 
! for example when computing the integral forces
!-------------------------------------------------------------------------------
subroutine add_cavity()
  use mpi
  use fsi_vars
  use penalization ! mask array etc
  implicit none
  integer :: ix,iy,iz
  real(kind=pr) :: ux,uy,uz
  
  if (iCavity=="yes") then
    ! this is the normal cavity forcing homogeneous dirichlet BC
    ux = 0.d0
    uy = 0.d0
    uz = 0.d0
  elseif (iCavity=="freestream") then
    ! cavity to force free stream inflow
    ux = uxmean
    uy = uymean
    uz = uzmean
  endif
  
  do ix = ra(1), rb(1)
    do iy = ra(2), rb(2)
       do iz = ra(3), rb(3)
       
          if ((ix<=cavity_size-1).or.(ix>=nx-1-cavity_size+1)) then            
            mask(ix,iy,iz) = 1.d0
            us(ix,iy,iz,:) = (/ux,uy,uz/)
            ! external boxes have color 0 (important for forces)
            mask_color(ix,iy,iz) = 0
          endif
          
          if ((iy<=cavity_size-1).or.(iy>=ny-1-cavity_size+1)) then       
            mask(ix,iy,iz) = 1.d0
            us(ix,iy,iz,:) = (/ux,uy,uz/)
            ! external boxes have color 0 (important for forces)
            mask_color(ix,iy,iz) = 0
          endif     
          
          if ((iz<=cavity_size-1).or.(iz>=nz-1-cavity_size+1)) then            
            mask(ix,iy,iz) = 1.d0
            us(ix,iy,iz,:) = (/ux,uy,uz/)
            ! external boxes have color 0 (important for forces)
            mask_color(ix,iy,iz) = 0
          endif
       enddo
    enddo
  enddo
 
end subroutine Add_Cavity
