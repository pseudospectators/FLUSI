!-------------------------------------------------------------------------------
! Passive scalar initial conditions, called from init_fields_fsi
! when resuming a backup file, this is automatically skipped (the scalar is then
! loaded from the runtime backup in init_fields_fsi-->Read_Runtime_Backup)
!
! INPUT/OUTPUT
!   scalars: the actual passive scalar (odor) fields
!   scalars_rhs: their right hand sides (set to zero here)
!-------------------------------------------------------------------------------
subroutine init_passive_scalar(scalars,scalars_rhs,Insect,beams,Wings)
  use mpi
  use vars
  use penalization ! mask array etc
  use solid_model
  use flexible_model
  use module_insects
  use passive_scalar_module
  use basic_operators
  implicit none

  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)

  integer :: ix,iy,iz,j
  real (kind=pr) :: x,y,z
  type(flexible_wing),dimension(1:nWings), intent(inout) :: Wings
  type(solid),dimension(1:nBeams), intent(inout) :: beams
  type(diptera),intent(inout)::Insect

  scalars = 0.d0
  scalars_rhs = 0.d0

  ! loop over passive scalars (we have up to nine different ones)
  do j=1, n_scalars

    if (mpirank==0) write(*,'("Initializing scalar #",i1)') j

    select case(scalar_props(j)%inicond)
    case("cosine")
      !-- set half the domain to one
      do iz=ra(3), rb(3)
        z = dble(iz)*dz
        do iy=ra(2), rb(2)
          y = dble(iy)*dy
          do ix=ra(1), rb(1)
            x = dble(ix)*dx - 0.5d0*xl
            scalars(ix,iy,iz,j) = dcos(pi*x) + dcos(4.d0*pi*x)
          enddo
        enddo
      enddo

  case("Kadoch2012")
      !-- set half the domain to one
      do iz = ra(3), rb(3)
          z = dble(iz)*dz - 1.0_pr - length
          do iy = ra(2), rb(2)
              y = dble(iy)*dy - 1.0_pr - length

              scalars(:,iy,iz,j) = cos(pi*z)*( cos(4.0_pr*pi*y) + cos(pi*y) )
          enddo
      enddo

    case("right_left_discontinuous")
      !---------------------------------------------------------------------------
      ! one half of the domain is 1, the other is 0, divided along y-axis
      ! no scalar inside mask (we build this here since FLUSI first loads inicond,
      ! then creates the mask!)
      !---------------------------------------------------------------------------
      call create_mask( 0.d0, Insect,beams,Wings )

      !-- set half the domain to one
      do iz=ra(3), rb(3)
        z = dble(iz)*dz
        do iy=ra(2), rb(2)
          y = dble(iy)*dy
          do ix=ra(1), rb(1)
            x = dble(ix)*dx
            if (y>0.5d0*yl) then
              scalars(ix,iy,iz,j) = 1.d0
            endif
          enddo
        enddo
      enddo

      !-- kill scalar inside obstacle
      scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) = &
      scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) * (1.0-mask)

    case ("right_left_smooth")
      !---------------------------------------------------------------------------
      ! smoothed heaviside function, covering approx. half the domain
      !---------------------------------------------------------------------------
      call create_mask( 0.d0, Insect,beams,Wings )

      !-- set half the domain to one
      do iz=ra(3), rb(3)
        z = dble(iz)*dz
        do iy=ra(2), rb(2)
          y = dble(iy)*dy
          do ix=ra(1), rb(1)
            x = dble(ix)*dx
            call smoothstep( scalars(ix,iy,iz,j), abs(y-0.5*yl), 0.25*yl, 4.0*dy)
          enddo
        enddo
      enddo

      !-- kill scalar inside obstacle
      scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) = scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) * (1.0-mask)

    case ("right_left_smooth2")
      !---------------------------------------------------------------------------
      ! smoothed heaviside function, covering approx. half the domain
      !---------------------------------------------------------------------------
      call create_mask( 0.d0, Insect,beams,Wings )

      !-- set half the domain to one
      do iz=ra(3), rb(3)
        z = dble(iz)*dz
        do iy=ra(2), rb(2)
          y = dble(iy)*dy
          do ix=ra(1), rb(1)
            x = dble(ix)*dx
            call smoothstep( scalars(ix,iy,iz,j), abs(z-0.5*zl), 0.25*zl, 2.0*dz)
          enddo
        enddo
      enddo

      !-- kill scalar inside obstacle
      scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) = scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) * (1.0-mask)

    case ("empty")
      !---------------------------------------------------------------------------
      !-- initialize scalar zero everywhere
      scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) = 0.d0


    case default
      if (mpirank==0) then
        write(*,*) "init_passive_scalar:: unknown inicond_scalar="//scalar_props(j)%inicond
        call abort(777811,"init_passive_scalar:: unknown inicond_scalar")
      endif
    end select
  end do

end subroutine init_passive_scalar
