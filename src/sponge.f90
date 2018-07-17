! ------------------------------------------------------------------------------
! Vorticity sponge technology
!
! adds a damping term to the vorticity and thus removes partly the periodicity
!
! input is penalized vorticity in fourier space, output is the penalization
! term to be added to the Navier--Stokes eqn.
!
! the vorticity sponge term is chi_sp * (vort) / eta_sponge
! but this is in the vorticity formulation
! so we compute the streamfunctions vort = - LAPLACE(psi)
! and then take the curl sponge = nabla \cross psi
!
! INPUT:
!       vort: the vorticity vector in x-space
!       work1: real valued work array (that will hold penalized vorticity)
! OUTPUT:
!       workc: cmplx work array for sponge term (vector)
! TO DO:
!       merge with vorticity2velocity in init_fields_fsi
! ------------------------------------------------------------------------------
subroutine vorticity_sponge( vort, work1, workc, Insect )
  use mpi
  use p3dfft_wrapper
  use module_insects
  use vars
  implicit none
  complex (kind=pr) :: im, spx,spy,spz
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  integer :: i, ix, iy, iz
  real(kind=pr),intent(in):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout):: work1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  type(diptera),intent(inout) :: Insect

  if (iVorticitySponge == "yes") then
    ! loop over components
    do i=1,3
      ! penalize this vorticity component in the work array
      call penalize_vort ( work1, vort(:,:,:,i), Insect )
      ! transform it to fourier space and store it in the sponge array workc
      call fft ( inx=work1, outk=workc(:,:,:,i)  )
    enddo


    ! imaginary unit
    im=dcmplx(0.d0,1.d0)

    do ix=ca(3), cb(3)
      kx=wave_x(ix)
      kx2=kx*kx

      do iy=ca(2), cb(2)
        ky=wave_y(iy)
        ky2=ky*ky

        do iz=ca(1), cb(1)
          kz=wave_z(iz)
          kz2=kz*kz

          k_abs_2=kx2+ky2+kz2
          if (abs(k_abs_2) .ne. 0.0) then
            ! we first "solve" the poisson eqn
            ! which gives us the streamfunction components
            spx = workc(iz,iy,ix,1) / k_abs_2
            spy = workc(iz,iy,ix,2) / k_abs_2
            spz = workc(iz,iy,ix,3) / k_abs_2

            ! we then take the curl of the streamfunction
            workc(iz,iy,ix,1)=im*(ky*spz - kz*spy)
            workc(iz,iy,ix,2)=im*(kz*spx - kx*spz)
            workc(iz,iy,ix,3)=im*(kx*spy - ky*spx)

          else
            workc(iz,iy,ix,1)=dcmplx(0.d0,0.d0)
            workc(iz,iy,ix,2)=dcmplx(0.d0,0.d0)
            workc(iz,iy,ix,3)=dcmplx(0.d0,0.d0)
          endif
        enddo
      enddo
    enddo
  endif
end subroutine vorticity_sponge


!-------------------------------------------------------------------------------
! vorticity penalization
!
! computes chi_sponge * (vort-vort0) / eta_sponge in physical space
! we currently do not allocate a mask_sponge array
!
! for one component only
!-------------------------------------------------------------------------------
subroutine penalize_vort ( vort_penalized, vort, Insect )
  use mpi
  use vars
  use penalization ! mask array etc
  use module_insects
  implicit none
  type(diptera),intent(inout) :: Insect
  ! input: vorticity in phys space
  real(kind=pr),intent(in):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  ! output: penalized vorticity in phys space
  real(kind=pr),intent(out):: vort_penalized(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: eps_inv,x(1:3)
  integer :: ix,iy,iz


  eps_inv= 1.d0 / eps_sponge
  vort_penalized = 0.d0

  select case (iSpongeType)
  case ("moving_insect")
    ! this sponge moves with the center of gravity Insect%xc of an insect or particle
    ! it is used since e.g a falling sphere that periodically re-enters the domain
    ! hits its own wake if no sponge is applied.
    ! on the other hand, traditional sponges were fixed in space (e.g. at the faces
    ! of the domain), thus the insect or particle would move through the sponge,
    ! which is very unphysical. thus the sponge has to move with the center of gravity
    ! MAY CAUSE TROUBLE IF MORE THAN ONE IS PRESENT!
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          x = (/dble(ix)*dx,dble(iy)*dy,dble(iz)*dz/) - Insect%xc_body_g
          x = periodize_coordinate(x, (/xl,yl,zl/))

          if ( (x(1)>xl/2.d0-dble(sponge_thickness)*dx).or.(x(1)<-xl/2.d0+dble(sponge_thickness)*dx)  ) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          if ( (x(2)>yl/2.d0-dble(sponge_thickness)*dy).or.(x(2)<-yl/2.d0+dble(sponge_thickness)*dy)  ) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          if ( (x(3)>zl/2.d0-dble(sponge_thickness)*dz).or.(x(3)<-zl/2.d0+dble(sponge_thickness)*dz)  ) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
        enddo
      enddo
    enddo

  case ("moving_insect_vertical")
    ! this sponge moves with the center of gravity Insect%xc of an insect or particle
    ! it is used since e.g a falling sphere that periodically re-enters the domain
    ! hits its own wake if no sponge is applied.
    ! on the other hand, traditional sponges were fixed in space (e.g. at the faces
    ! of the domain), thus the insect or particle would move through the sponge,
    ! which is very unphysical. thus the sponge has to move with the center of gravity
    ! MAY CAUSE TROUBLE IF MORE THAN ONE IS PRESENT!
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          x = (/dble(ix)*dx,dble(iy)*dy,dble(iz)*dz/) - Insect%xc_body_g
          x = periodize_coordinate(x, (/xl,yl,zl/))
          if ( (x(3)>zl/2.d0-dble(sponge_thickness)*dz).or.(x(3)<-zl/2.d0+dble(sponge_thickness)*dz)  ) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
        enddo
      enddo
    enddo



  case ("cavity")
    !--------------------------------------------
    ! sponge for the cavity type (ie we set a solid wall
    ! around the domain and kill the vorticity in front
    ! of if). But you can also use it without the solid wall
    ! ie iCavity=no and iSponge=yes, iSpongeType=cavity
    !--------------------------------------------
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1.0d-12) then
            !------------------------
            if (nx>4) then ! skip this direction for 2D runs...
            if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then
              vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
            endif
            endif
            !------------------------
            if ((iy<=sponge_thickness-1).or.(iy>=ny-1-sponge_thickness+1)) then
              vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
            endif
            !------------------------
            if ((iz<=sponge_thickness-1).or.(iz>=nz-1-sponge_thickness+1)) then
              vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
            endif
            !------------------------
          endif
        enddo
      enddo
    enddo

  case ("cavity_open_z")
    !--------------------------------------------
    ! sponge for the cavity type (ie we set a solid wall
    ! around the domain and kill the vorticity in front
    ! of if). But you can also use it without the solid wall
    ! ie iCavity=no and iSponge=yes, iSpongeType=cavity
    !--------------------------------------------
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1.0d-12) then
            !------------------------
            if (nx>4) then ! skip this direction for 2D runs...
            if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then
              vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
            endif
            endif
            !------------------------
            if ((iy<=sponge_thickness-1).or.(iy>=ny-1-sponge_thickness+1)) then
              vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
            endif
            !------------------------
          endif
        enddo
      enddo
    enddo

  case ("xmin_xmax_ymin_ymax")
    !--------------------------------------------
    ! Dmitry, 25 Oct 2013
    ! At xmin, xmax, ymin, ymax walls only. No sponge at zmin, zmax walls.
    !--------------------------------------------
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1.0d-12) then
          if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif

          if ((iy<=sponge_thickness-1).or.(iy>=ny-1-sponge_thickness+1)) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          endif
        enddo
      enddo
    enddo

  case ("xmin_xmax")
    !--------------------------------------------
    ! Dmitry, 25 Oct 2013
    ! At xmin, xmax walls only.
    !--------------------------------------------
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1.0d-12) then
          if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          endif
        enddo
      enddo
    enddo

  case ("outlet_x")
    !--------------------------------------------
    ! sponge is one at right outflow, ie for ix>nx-1-sponge_thickness+1
    !--------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1e-12) then
          if (ix>=nx-1-sponge_thickness+1) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          endif
        enddo
      enddo
    enddo
  case ("xmin_xmax_zmin_zmax")
    !--------------------------------------------
    ! Dmitry, 25 Oct 2013
    ! At xmin, xmax, zmin, zmax walls only.
    !--------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        do ix=ra(1),rb(1)
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1.0d-12) then
          if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif

          if ((iz<=sponge_thickness-1).or.(iz>=nz-1-sponge_thickness+1)) then
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          endif
        enddo
      enddo
    enddo

  case ("2D_frame")
    !--------------------------------------------
    ! Thomas, 3 Jul 2014
    ! frame around domain, "cavity" for 2D. note nx=1.
    !--------------------------------------------
    do iz = ra(3), rb(3)
      do iy = ra(2), rb(2)
        if ((iy<=sponge_thickness-1).or.(iy>=ny-1-sponge_thickness+1)) then
          vort_penalized(:,iy,iz) = -vort(:,iy,iz)*eps_inv
        endif

        if ((iz<=sponge_thickness-1).or.(iz>=nz-1-sponge_thickness+1)) then
          vort_penalized(:,iy,iz) = -vort(:,iy,iz)*eps_inv
        endif
      enddo
    enddo

  case ("top_cover")
    !--------------------------------------------
    ! sponge as a cover on top of the domain
    ! (top=positive z)
    !--------------------------------------------
    do iz=ra(3),rb(3)
      ! note we currenly do not allocate a mask for this
      if ( iz>nz-sponge_thickness ) then
        vort_penalized(:,:,iz) = - vort(:,:,iz)*eps_inv
      endif
    enddo
  end select
end subroutine penalize_vort
