! Spherical obstacle
subroutine draw_sphere
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, tmp, R, N_smooth, smoothing, safety

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  N_smooth = 1.5d0
  ! thickness of smoothing layer and safety distance
  if (nx==1) then
    smoothing = N_smooth*max(dy,dz)
    safety = 2.d0*N_smooth*max(dy,dz)
  else
    smoothing = N_smooth*max(dx,dy,dz)
    safety = 2.d0*N_smooth*max(dx,dy,dz)
  endif

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        x=dble(ix)*dx
        y=dble(iy)*dy
        z=dble(iz)*dz
        R = dsqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
        if ( R <= 0.5d0*length+safety ) then
          call SmoothStep (tmp, R, 0.5d0*length , smoothing)
          mask (ix, iy, iz) = tmp

          ! assign color "1" where >0 indicates something "useful"
          if (tmp > 1.0e-12) mask_color(ix,iy,iz) = 1
        endif
      enddo
    enddo
  enddo
end subroutine draw_sphere


! draws a circular cylinder, the axis is the x-axis (then it can be used for 2D
! runs as well)
subroutine draw_cylinder_x
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, tmp, R, N_smooth

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  N_smooth = 1.5d0

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        y=dble(iy)*dy
        z=dble(iz)*dz
        R = dsqrt( (y-y0)**2 + (z-z0)**2 )
        if ( R <= 0.5d0*length+2.d0*N_smooth*max(dy,dz) ) then
          call SmoothStep (tmp, R, 0.5d0*length , N_smooth*max(dy,dz))
          mask (ix, iy, iz) = tmp
          us(ix,iy,iz,1) = 0.d0
          us(ix,iy,iz,2) = 0.d0
          us(ix,iy,iz,3) = 0.d0

          ! assign color "1" where >0 indicates something "useful"
          if (tmp > 1.0e-12) mask_color(ix,iy,iz) = 1
        endif
      enddo
    enddo
  enddo
end subroutine draw_cylinder_x

! draws a moving circular cylinder, the axis is the x-axis (then it can be used for 2D
! runs as well)
subroutine draw_moving_cylinder_x (time)
  use vars
  use vars
  use penalization ! mask array etc
  implicit none

  real(kind=pr),intent(in) :: time
  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, tmp, R, N_smooth, Vc, yc

  ! Velocity and positions of the cylinder
  Vc = -1.0d0
  yc = y0 + Vc * time

  N_smooth = 1.5d0

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)
        y=dble(iy)*dy
        z=dble(iz)*dz
        R = dsqrt( (y-yc)**2 + (z-z0)**2 )
        if ( R <= 0.5d0*length+2.d0*N_smooth*max(dy,dz) ) then
          call SmoothStep (tmp, R, 0.5d0*length , N_smooth*max(dy,dz))
          mask (ix, iy, iz) = tmp
          us(ix,iy,iz,1) = 0.d0
          us(ix,iy,iz,2) = Vc
          us(ix,iy,iz,3) = 0.d0

          ! assign color "1" where >0 indicates something "useful"
          if (tmp > 1.0e-12) mask_color(ix,iy,iz) = 1
        endif
      enddo
    enddo
  enddo
end subroutine draw_moving_cylinder_x


! draws a vertical circular cylinder (z direction) with a turbulizer.
! No wall smoothing.
subroutine draw_cylinder_asym_z
  use mpi
  use vars
  use penalization ! mask array etc
  implicit none

  integer :: ix, iy, iz
  real (kind=pr) :: x, y, z, Rsqr

  ! Defaults for the similar mast as a channel are
  ! length = 2*thick_wall = 2*0.947d0
  ! x0 = pos_wall = 1.5d0
  ! y0 = 0.5d0*yl

  ! reset everything
  mask = 0.d0
  mask_color = 0
  us = 0.d0

  do iz=ra(3),rb(3)
    do iy=ra(2),rb(2)
      do ix=ra(1),rb(1)

        x = dble(ix)*dx
        y = dble(iy)*dy
        z = dble(iz)*dz

        ! Main cylinder
        Rsqr = (x-x0)**2 + (y-y0)**2
        if ( Rsqr <= 0.5d0**2*length**2 ) then
          mask (ix, iy, iz) = 1.0d0
          us(ix,iy,iz,1) = 0.d0
          us(ix,iy,iz,2) = 0.d0
          us(ix,iy,iz,3) = 0.d0

          ! assign color "1" where >0 indicates something "useful"
          mask_color(ix,iy,iz) = 1
        endif

        ! Small cylinder to break symmetry
        Rsqr = (x-(x0+0.5d0*length*dcos(0.25d0*pi)))**2+ &
           (y-(y0+0.5d0*length*dsin(0.25d0*pi)))**2
        if ( ( Rsqr < 0.05d0**2 ) .and. ( dabs(z) > 0.05d0 ) ) then
          mask (ix, iy, iz) = 1.0d0
          us(ix,iy,iz,1) = 0.d0
          us(ix,iy,iz,2) = 0.d0
          us(ix,iy,iz,3) = 0.d0

          ! assign color "1" where >0 indicates something "useful"
          mask_color(ix,iy,iz) = 1
        endif

      enddo
    enddo
  enddo
end subroutine draw_cylinder_asym_z

