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
! and then take the curl sponge = nabla \crossproduct psi
!
! INPUT: 
!       work: real valued work array that will hold penalized vorticity
!       vort: the vorticity vector
! OUTPUT:
!       sponge (global): the vorticity sponge term in Fourier space
!
! TO DO:
!       merge with vorticity2velocity in init_fields_fsi
! ------------------------------------------------------------------------------
subroutine vorticity_sponge( work, vort )
  use mpi_header
  use fsi_vars  
  implicit none  
  complex (kind=pr) :: im, spx,spy,spz
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  integer :: i, ix, iy, iz
  
  real(kind=pr),intent(inout):: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(in):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  
  if (iVorticitySponge == "yes") then    
    ! loop over components
    do i=1,3  
      ! penalize this vorticity component in the work array
      call penalize_vort ( work, vort(:,:,:,i) )
      ! then transform it to fourier space and store it in the global sponge array
      call fft ( sponge(:,:,:,i), work )  
    enddo  
    
    
    ! imaginary unit
    im=dcmplx(0.d0,1.d0)    
    do iy=ca(3), cb(3)    ! ky : 0..ny/2-1 ,then, -ny/2..-1     
      ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
      ky2=ky*ky
      do ix=ca(2), cb(2)  ! kx : 0..nx/2
        kx=scalex*dble(ix)                
        kx2=kx*kx
        do iz=ca(1),cb(1)  ! kz : 0..nz/2-1 ,then, -nz/2..-1           
          kz     =scalez*dble(modulo(iz+nz/2,nz)-nz/2)
          kz2    =kz*kz
          k_abs_2=kx2+ky2+kz2
          if (abs(k_abs_2) .ne. 0.0) then  
            ! we first "solve" the poisson eqn 
            ! which gives us the streamfunction components
            spx = sponge(iz,ix,iy,1) / k_abs_2
            spy = sponge(iz,ix,iy,2) / k_abs_2
            spz = sponge(iz,ix,iy,3) / k_abs_2
            
            ! we then take the curl of the streamfunction
            sponge(iz,ix,iy,1)=im*(ky*spz - kz*spy)
            sponge(iz,ix,iy,2)=im*(kz*spx - kx*spz)
            sponge(iz,ix,iy,3)=im*(kx*spy - ky*spx)          
            
          else
            sponge(iz,ix,iy,1)=dcmplx(0.d0,0.d0)
            sponge(iz,ix,iy,2)=dcmplx(0.d0,0.d0)
            sponge(iz,ix,iy,3)=dcmplx(0.d0,0.d0)
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
!
! in this version, the mask is applied in a layer on top (in Z-direction)
! of the domain.
!-------------------------------------------------------------------------------
subroutine penalize_vort ( vort_penalized, vort )
  use mpi_header
  use fsi_vars
  
  ! input: vorticity in phys space
  real(kind=pr),intent(in):: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  ! output: penalized vorticity in phys space
  real(kind=pr),intent(out):: vort_penalized(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: eps_inv
  integer :: iz
  eps_inv= 1.d0 / eps_sponge
  vort_penalized = 0.d0
  
  select case (iSpongeType)
  case ("cavity")
    !--------------------------------------------
    ! sponge for the cavity type (ie we set a solid wall
    ! around the domain and kill the vorticity in front
    ! of if). But you can also use it without the solid wall
    ! ie iCavity=no and iSponge=yes, iSpongeType=cavity
    !--------------------------------------------
    do ix = ra(1), rb(1)
      do iy = ra(2), rb(2)
        do iz = ra(3), rb(3) 
          ! do not use vorticity sponge and solid wall simulateously
          if (mask(ix,iy,iz) < 1e-12) then
          if ((ix<=sponge_thickness-1).or.(ix>=nx-1-sponge_thickness+1)) then            
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          
          if ((iy<=sponge_thickness-1).or.(iy>=ny-1-sponge_thickness+1)) then            
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif     
          
          if ((iz<=sponge_thickness-1).or.(iz>=nz-1-sponge_thickness+1)) then            
            vort_penalized(ix,iy,iz) = -vort(ix,iy,iz)*eps_inv
          endif
          endif
        enddo
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