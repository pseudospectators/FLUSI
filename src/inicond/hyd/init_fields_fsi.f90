! Set initial conditions for fsi code.
subroutine init_fields_fsi(time,u,nlk,work,mask,mask_color,us,Insect,beams)
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none

  type(timetype),intent(inout) :: time
  real(kind=pr),intent(inout)::u(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  real(kind=pr),intent(inout)::nlk(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq,1:nrhs)
  real(kind=pr),intent(inout)::work(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:nrw)
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  type(solid), dimension(1:nBeams),intent(inout) :: beams
  type(diptera), intent(inout) :: Insect
    
  integer :: ix,iy,iz
  real (kind=pr) :: x,y,z,r,gamma0,x00,r00,omega,a,b
  real (kind=pr) :: R1,R2,uu

  ! Assign zero values
  time%time = 0.0d0
  time%dt_new  = tsave
  time%it = 0
  
  u = 0.d0
  us = 0.d0
  nlk = 0.d0
  work = 0.d0
  mask = 0.d0
  mask_color = 0
  
  
  select case(inicond)
  case ("couette")
    !--------------------------------------------------
    ! couette flow
    !--------------------------------------------------  
    R1=0.5d0
    R2=1.0d0
    omega=1.25d0
    
    a = omega*(-R1**2 / (R2**2 - R1**2))
    b = omega*(R1**2 * R2**2) / (R2**2 - R1**2)
    
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy - 0.5d0*yl
        z = dble(iz)*dz - 0.5d0*zl
        R = dsqrt(y**2 + z**2)
        
        if ((R>R1).and.(R<R2)) then
          ! fluid domain
          uu = a*R + b/R
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) =+uu*z/R 
          u(:,iy,iz,3) =-uu*y/R
        elseif (R>=R2) then
          ! outer cylinder
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) = 0.d0
          u(:,iy,iz,3) = 0.d0
        elseif (R<=R1) then
          ! inner cylinder
          u(:,iy,iz,1) = 0.d0
          u(:,iy,iz,2) = +omega*z
          u(:,iy,iz,3) = -omega*y
        endif
      enddo
    enddo
    
  case ("taylor_green_2d")
    !--------------------------------------------------
    ! taylor green vortices
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy
        z = dble(iz)*dz
        u(:,iy,iz,2) = dsin( y ) * dcos( z )
        u(:,iy,iz,3) =-dcos( y ) * dsin( z )
        ! remember that p is the total pressure, i.e. p+0.5u*u
        u(:,iy,iz,4) =0.25d0*(dcos(2.d0*y)+dcos(2.d0*z)) + &
                      0.5d0*(u(:,iy,iz,2)**2+u(:,iy,iz,3)**2)
      enddo
    enddo
    
  case("infile")
     !--------------------------------------------------
     ! read HDF5 files
     !--------------------------------------------------  
     if (mpirank==0) write (*,*) "*** inicond: reading infiles"
     call Read_Single_File ( file_ux, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1) )
     call Read_Single_File ( file_uy, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),2) )
     call Read_Single_File ( file_uz, u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),3) )
     call Read_Single_File ( file_p , u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),4) )
     if (mpirank==0) write (*,*) "*** done reading infiles"

  case("MeanFlow")
     !--------------------------------------------------
     ! mean flow only
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: mean flow"
     
     u(:,:,:,1) = Uxmean
     u(:,:,:,2) = Uymean
     u(:,:,:,3) = Uzmean
       
  case("quiescent")
     !--------------------------------------------------
     ! fluid at rest
     !--------------------------------------------------  
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     u = 0.d0

  case default
     if(inicond(1:8) == "backup::") then
        !--------------------------------------------------
        ! read from backup
        !--------------------------------------------------   
        if (mpirank==0) write (*,*) "*** inicond: retaking backup " // &
             inicond(9:len(inicond))
        call read_runtime_backup(inicond(9:len(inicond)),time,u,Insect,beams)
     else
        !--------------------------------------------------
        ! unknown inicond : error
        !--------------------------------------------------
        if (mpirank==0) write (*,*) inicond
        if (mpirank==0) write (*,*) '??? ERROR: Invalid initial condition'
        call abort()
     endif
  end select

  
  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the solid solver
  !-----------------------------------------------------------------------------
  if (use_solid_model=="yes") then
    if(mpirank==0) write(*,*) "Initializing solid solver and testing..."
    call init_beams( beams )
    call show_beam(beams(1))
    call surface_interpolation_testing( time%time, beams(1),work(:,:,:,1),mask,mask_color,us )
    call init_beams( beams )
  endif
  
  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the passive scalar(s)
  !-----------------------------------------------------------------------------
!   if ((use_passive_scalar==1).and.(index(inicond,"backup::")==0)) then
!     ! only if not resuming a backup
!     call init_passive_scalar(uk(:,:,:,4),vort,workc(:,:,:,1),Insect,beams)
!   endif
  
end subroutine init_fields_fsi


! Computes the divergence-free velocity in Fourier space u given vort
! in physical space.  work is a work array
subroutine Vorticity2Velocity(uk,work,vort)
  use vars
  use p3dfft_wrapper
  implicit none

  complex (kind=pr),intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout)::work(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  real (kind=pr), intent (in) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  integer :: ix, iy, iz, i
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  complex (kind=pr) :: im
  ! imaginary unit
  im=dcmplx(0.d0,1.d0)
  !-------------------------------------------------
  ! Compute vorticity in Fourier space
  !-------------------------------------------------
  do i=1,3
     call fft(work(:,:,:,i),vort(:,:,:,i))
  enddo
  
  !-------------------------------------------------
  ! Compute streamfunction in Fourier space
  ! work(:,:,:,1:3, 1) will contain the three components of
  ! streamfunction
  !------------------------------------------------- 
  do iz=ca(1), cb(1)
    kz=wave_z(iz)
    kz2=kz*kz
    do iy=ca(2), cb(2)
      ky=wave_y(iy)
      ky2=ky*ky
      do ix=ca(3),cb(3)
        kx     =wave_x(ix)
        kx2    =kx*kx
        k_abs_2=kx2+ky2+kz2
        if (abs(k_abs_2) .ne. 0.0) then  
          work(iz,iy,ix,1)=-work(iz,iy,ix,1) / k_abs_2
          work(iz,iy,ix,2)=-work(iz,iy,ix,2) / k_abs_2
          work(iz,iy,ix,3)=-work(iz,iy,ix,3) / k_abs_2
        else
          work(iz,iy,ix,1)=dcmplx(0.d0,0.d0)
          work(iz,iy,ix,2)=dcmplx(0.d0,0.d0)
          work(iz,iy,ix,3)=dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo

  !-----------------------------------------------
  !-- compute velocity as curl of streamfunction
  !-----------------------------------------------
  do iz=ca(1),cb(1)          
    !-- wavenumber in z-direction
    kz = wave_z(iz)      
    do iy=ca(2), cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)      
      do ix=ca(3), cb(3)
        !-- wavenumber in x-direction
        kx = wave_x(ix)
        uk(iz,iy,ix,1)=im*(ky*work(iz,iy,ix,3) - kz*work(iz,iy,ix,2))
        uk(iz,iy,ix,2)=im*(kz*work(iz,iy,ix,1) - kx*work(iz,iy,ix,3))
        uk(iz,iy,ix,3)=im*(kx*work(iz,iy,ix,2) - ky*work(iz,iy,ix,1))
      enddo
    enddo
  enddo
end subroutine Vorticity2Velocity
