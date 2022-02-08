module passive_scalar_module
  use mpi
  use vars
  use penalization
  use basic_operators
  use module_helpers
  use ghosts
  implicit none

  type scalar_params
    real(kind=pr) :: kappa
    real(kind=pr) :: eps, source_x0(1:4)
    character(len=strlen) :: inicond, sourceterm
    integer :: imin(1:3), imax(1:3)
  end type


  ! module global variables
  type(scalar_params),allocatable,dimension(:) :: scalar_props


contains



  !-------------------------------------------------------------------------------
  ! right hand side for passive scalars. we compute it using finite differences
  ! (AND NOT SPECTRAL) since this is cheaper. we treat up to 9 different scalars
  ! The discretization in space is 4th order central, first derivatives treated
  ! with optimized scheme.
  !-------------------------------------------------------------------------------
  subroutine cal_nlk_scalar( time, it, u, scalars, scalars_rhs )
    implicit none

    real(kind=pr),intent(in) :: time
    integer, intent(in) :: it
    real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
    real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
    real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
    real(kind=pr),allocatable,dimension(:,:,:) :: mask2
    real(kind=pr) :: kx,ky,kz,maxi,t1
    real(kind=pr) :: k(1:3) ! wavenumber vector
    complex(kind=pr) :: imag
    integer :: ix,iy,iz,id,mpicode,j

    real(kind=pr)::a(-3:+3)
    real(kind=pr)::b1,b2,b3,b4,b5,dxinv,dyinv,dzinv,dx2inv,dy2inv,dz2inv,ux,uy,uz,&
    usx,usy,usz,wx,wy,wz,gx,gy,gz,D,chi,chidx,chidz,chidy,D_dx,D_dy,D_dz,gxx,gyy,gzz

    t1 = MPI_wtime()

    !-- initialization
    scalars_rhs = 0.0d0

    ! this is a fnite differences code with explicit discretization; therefore
    ! a redundant layer of ghost nodes has to be synchronzied between processes
    call synchronize_ghosts(scalars,n_scalars)

    ! we need to derive the mask function, but it has no ghost nodes. therefore,
    ! extend it to mask2 with ghosts and sync them
    allocate (mask2(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3)))
    mask2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) = mask
    call synchronize_ghosts(mask2)

    ! finite differences coefficients
    ! Tam & Webb, 4th order optimized (first derivative)
    a = (/-0.02651995d0, +0.18941314d0, -0.79926643d0, 0.0d0, &
    0.79926643d0, -0.18941314d0, 0.02651995d0/)

    ! for second derivatives:
    b1 =-1.d0/12.d0
    b2 = 4.d0/3.d0
    b3 =-5.d0/2.d0
    b4 = 4.d0/3.d0
    b5 =-1.d0/12.d0

    dxinv = 1.d0/dx
    dyinv = 1.d0/dy
    dzinv = 1.d0/dz

    dx2inv = 1.d0/(dx**2)
    dy2inv = 1.d0/(dy**2)
    dz2inv = 1.d0/(dz**2)

    ! loop over scalars
    do j = 1,n_scalars
      !-----------------------------------------------------------------------------
      ! Sanity test: if the scalar values are out of range, we skip the entire
      ! right hand side, but we do not abort (since the fluid is fine, it's just
      ! the scalar that fails.
      !-----------------------------------------------------------------------------
      if (compute_scalar) then
        !-- global maximum
        maxi = fieldmax( scalars(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),j) )
        !-- skip passive scalar from now on, if out of range
        if ( maxi > 1.0d4 ) then
          if (mpirank==0) then
            write(*,'("WARNING! PASSIVE SCALAR FAILED! time=",es12.4)') time
            write(*,'("THIS MEANS WE WILL NO LONGER COMPUTE IT FROM NOW ON!!!")')
          endif
          !-- kill run if scalar is considered crutial for results
          if (stop_on_fail=="yes") call abort(8181, "PASSIVE SCALAR FAILED")
          !-- otherwise, warn and continue without scalar
          compute_scalar = .false. ! callees will skip
        endif
      endif

      !------------------------------------------------------------------------------
      ! compute right hand side of passive scalar
      !------------------------------------------------------------------------------
      do iz=ra(3),rb(3)
        do iy=ra(2),rb(2)
          do ix=ra(1),rb(1)
            ux = u(ix,iy,iz,1)
            uy = u(ix,iy,iz,2)
            uz = u(ix,iy,iz,3)

            usx = us(ix,iy,iz,1)
            usy = us(ix,iy,iz,2)
            usz = us(ix,iy,iz,3)

            chi = mask2(ix,iy,iz)

            ! penalized diffusion coefficient
            D = scalar_props(j)%kappa*(1.d0-chi) + scalar_props(j)%eps*chi

            ! this is the vector in front of the gradient
            wx = -((1.d0-chi)*ux + chi*usx)
            wy = -((1.d0-chi)*uy + chi*usy)
            wz = -((1.d0-chi)*uz + chi*usz)

            ! gradient of passive scalar
            gx = (a(-3)*scalars(ix-3,iy,iz,j)&
            +a(-2)*scalars(ix-2,iy,iz,j)&
            +a(-1)*scalars(ix-1,iy,iz,j)&
            +a( 0)*scalars(ix  ,iy,iz,j)&
            +a(+3)*scalars(ix+3,iy,iz,j)&
            +a(+2)*scalars(ix+2,iy,iz,j)&
            +a(+1)*scalars(ix+1,iy,iz,j))*dxinv

            gy = ( a(-3)*scalars(ix,iy-3,iz,j)&
            +a(-2)*scalars(ix,iy-2,iz,j)&
            +a(-1)*scalars(ix,iy-1,iz,j)&
            +a( 0)*scalars(ix,iy  ,iz,j)&
            +a(+3)*scalars(ix,iy+3,iz,j)&
            +a(+2)*scalars(ix,iy+2,iz,j)&
            +a(+1)*scalars(ix,iy+1,iz,j))*dyinv

            gz = ( a(-3)*scalars(ix,iy,iz-3,j)&
            +a(-2)*scalars(ix,iy,iz-2,j)&
            +a(-1)*scalars(ix,iy,iz-1,j)&
            +a( 0)*scalars(ix,iy,iz  ,j)&
            +a(+3)*scalars(ix,iy,iz+3,j)&
            +a(+2)*scalars(ix,iy,iz+2,j)&
            +a(+1)*scalars(ix,iy,iz+1,j))*dzinv

            ! gradient of mask function ( we need that for the diffusive term)
            ! since this guy reads div( (kappa(1-mask) + eps*mask) * grad(phi) )
            ! so this boils down to d/dx (D*gx) = D_dx*gx + D*gxx
            ! so we need D_dx and this is kappa*(1-mask_dx)+ eps*mask_dx
            chidx = (a(-3)*mask2(ix-3,iy,iz)&
            +a(-2)*mask2(ix-2,iy,iz)&
            +a(-1)*mask2(ix-1,iy,iz)&
            +a( 0)*mask2(ix  ,iy,iz)&
            +a(+3)*mask2(ix+3,iy,iz)&
            +a(+2)*mask2(ix+2,iy,iz)&
            +a(+1)*mask2(ix+1,iy,iz))*dxinv

            chidy = ( a(-3)*mask2(ix,iy-3,iz)&
            +a(-2)*mask2(ix,iy-2,iz)&
            +a(-1)*mask2(ix,iy-1,iz)&
            +a( 0)*mask2(ix,iy  ,iz)&
            +a(+3)*mask2(ix,iy+3,iz)&
            +a(+2)*mask2(ix,iy+2,iz)&
            +a(+1)*mask2(ix,iy+1,iz))*dyinv

            chidz = ( a(-3)*mask2(ix,iy,iz-3)&
            +a(-2)*mask2(ix,iy,iz-2)&
            +a(-1)*mask2(ix,iy,iz-1)&
            +a( 0)*mask2(ix,iy,iz  )&
            +a(+3)*mask2(ix,iy,iz+3)&
            +a(+2)*mask2(ix,iy,iz+2)&
            +a(+1)*mask2(ix,iy,iz+1))*dzinv

            D_dx = scalar_props(j)%kappa*(1.d0-chidx) + scalar_props(j)%eps*chidx
            D_dy = scalar_props(j)%kappa*(1.d0-chidy) + scalar_props(j)%eps*chidy
            D_dz = scalar_props(j)%kappa*(1.d0-chidz) + scalar_props(j)%eps*chidz

            ! second derivatives of passive scalar
            gxx = (b1*scalars(ix-2,iy,iz,j)+b2*scalars(ix-1,iy,iz,j)+b3*scalars(ix,iy,iz,j)&
            +b4*scalars(ix+1,iy,iz,j)+b5*scalars(ix+2,iy,iz,j))*dx2inv
            gyy = (b1*scalars(ix,iy-2,iz,j)+b2*scalars(ix,iy-1,iz,j)+b3*scalars(ix,iy,iz,j)&
            +b4*scalars(ix,iy+1,iz,j)+b5*scalars(ix,iy+2,iz,j))*dy2inv
            gzz = (b1*scalars(ix,iy,iz-2,j)+b2*scalars(ix,iy,iz-1,j)+b3*scalars(ix,iy,iz,j)&
            +b4*scalars(ix,iy,iz+1,j)+b5*scalars(ix,iy,iz+2,j))*dz2inv

            ! assemble everything
            scalars_rhs(ix,iy,iz,j) = wx*gx + wy*gy + wz*gz & ! penalized convection term
            + D_dx*gx + D*gxx & ! penalized laplacian
            + D_dy*gy + D*gyy &
            + D_dz*gz + D*gzz
          enddo
        enddo
      enddo
    enddo

    call add_scalar_source( time, scalars, scalars_rhs )
    deallocate(mask2)
    
    ! save timing
    call toc( "Pasisve scalar (rhs)", MPI_wtime() - t1)
  end subroutine cal_nlk_scalar


  !-------------------------------------------------------------------------------
  ! source term for passive scalar, dirichlet condition
  ! note here we have to use eps instead of eps_scalar, since
  ! this term imposed a dt<eps stability condition
  !-------------------------------------------------------------------------------
  subroutine add_scalar_source ( time, scalars, scalars_rhs )
    use mpi
    use p3dfft_wrapper
    use vars
    implicit none

    real(kind=pr),intent(in) :: time
    real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
    real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
    real(kind=pr) :: x,y,z,source,masksource,R
    integer :: ix,iy,iz,j


    do j=1, n_scalars
      select case (scalar_props(j)%sourceterm)
      case ("gauss_blob")
        ! gaussian blob source
        do iz=ra(3),rb(3)
          do iy=ra(2),rb(2)
            do ix=ra(1),rb(1)
              x = dble(ix)*dx
              y = dble(iy)*dy
              z = dble(iz)*dz

              R = ( (x-scalar_props(j)%source_x0(1))**2 &
                   +(y-scalar_props(j)%source_x0(2))**2 &
                   +(z-scalar_props(j)%source_x0(3))**2 )

              masksource = dexp( -R / (scalar_props(j)%source_x0(4))  )

              if (masksource > 1.0d-4) then
                source = -masksource*(scalars(ix,iy,iz,j)-1.d0) / eps
                scalars_rhs(ix,iy,iz,j) = scalars_rhs(ix,iy,iz,j) + source
              endif
            enddo
          enddo
        enddo
      case ("upstream_random")
        ! in this case, we inject in a layer upstream x<x0 a randomized scalar
        ! with mean 0.5


      case ("none","empty","no","without")
        ! simply do not do anything
      case default
        call abort(0909090,'Passive scalar, source term undefined.')
      end select
    enddo

  end subroutine add_scalar_source


end module
