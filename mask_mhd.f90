! MHD wrapper for different mask functions 
subroutine create_mask_mhd()
  use mpi_header
  use mhd_vars
  implicit none

  if(iPenalization == 1) then
     select case (iMask)
     case("taylorcouette")
        call  tc_mask_mhd()
     case("smc")
        call smc_mask_mhd()
     case default
        if(mpirank == 0) then
           write (*,*) &
                "iMask not properly set for mhd in create_mask_mhd; stopping."
           stop
        endif
     end select
  endif
end subroutine create_mask_mhd


! MHD wrapper for setting (possibly velocity-dependent) imposed field.
subroutine update_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none

  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  if(iPenalization == 1) then
     select case (iMask)
     case("taylorcouette")
        call tc_us_mhd()
     case("smc")
        call smc_us_mhd(ub)
     case default
        if(mpirank == 0) then
           write (*,*) &
                "iMask not properly set for mhd in update_us_mhd; stopping."
           stop
        endif
     end select
  endif
end subroutine update_us_mhd


! Set the solid velocity for MHD Taylor-Couette flow.
subroutine tc_us_mhd()
  use mpi_header
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz
  
  us=0.d0

  do ix = ra(1), rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy = ra(2), rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)

        if(r <= R1) then
           do iz = ra(3), rb(3)
              us(ix,iy,iz,1)=0.d0 
           enddo
        endif
        if(r >= R2) then
           do iz = ra(3), rb(3)
              us(ix,iy,iz,1)=0.d0
           enddo
        endif

     enddo
  enddo
end subroutine tc_us_mhd


! Set the mask function for MHD Taylor-Couette flow.
subroutine tc_mask_mhd()
  use mpi_header
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz

  mask=0.d0

  do ix = ra(1), rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy = ra(2), rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)

        if(r <= R1 .or. r >= R2) then
           do iz = ra(3), rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine tc_mask_mhd


! Set the solid velocity for Sean-Montgomery-Chen flow.
subroutine smc_us_mhd(ub)
  use mpi_header
  use mhd_vars
  implicit none
  
  real(kind=pr),intent(in)::ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz
  
  us=0.d0

  do ix = ra(1), rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy = ra(2), rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x + y*y)

        if(r >= R1) then
           do iz = ra(3), rb(3)

              ! Velocity is no-slip:
              us(ix,iy,iz,1)=0.d0
              us(ix,iy,iz,2)=0.d0
              us(ix,iy,iz,3)=0.d0

              ! Magnetic fields is penalized to a particular value:
              ! FIXME: add smoothing, stuff here.
              us(ix,iy,iz,4)=ub(ix,iy,iz,4) + y*Bc/R1
              us(ix,iy,iz,5)=ub(ix,iy,iz,5) - x*Bc/R1
              us(ix,iy,iz,6)=B0
           enddo
        endif
     enddo
  enddo
end subroutine smc_us_mhd


! Set the mask function for Sean-Montgomery-Chen flow.
subroutine smc_mask_mhd()
  use mpi_header
  use mhd_vars
  implicit none
  
  real (kind=pr) :: r, x, y
  integer :: ix, iy, iz

  mask=0.d0

  do ix = ra(1), rb(1)  
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy = ra(2), rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)

        if(r >= R1) then
           do iz = ra(3), rb(3)
              mask(ix,iy,iz)=1.d0
           enddo
        endif

     enddo
  enddo
end subroutine smc_mask_mhd
