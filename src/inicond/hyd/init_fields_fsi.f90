! Set initial conditions for fsi code.
subroutine init_fields_fsi(time,it,dt0,dt1,n0,n1,uk,nlk,vort,explin,workc,press,Insect,beams)
  use mpi
  use fsi_vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none

  integer,intent (inout) :: n1,it,n0
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  ! the workc array is not always allocated, ensure allocation before using
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:1)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))  
  type(solid),dimension(1:nBeams), intent(inout) :: beams
  type(diptera),intent(inout)::Insect 
  integer :: ix,iy,iz
  real (kind=pr) :: x,y,z,r,a,b,gamma0,x00,r00,omega
  real (kind=pr) :: uu

  ! Assign zero values
  time = 0.0d0
  dt1  = tsave
  it = 0
  
  uk = dcmplx(0.0d0,0.0d0)
  nlk = dcmplx(0.0d0,0.0d0)
  explin = 0.0
  vort = 0.0d0
  
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
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) =+uu*z/R 
          vort(:,iy,iz,3) =-uu*y/R
        elseif (R>=R2) then
          ! outer cylinder
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) = 0.d0
          vort(:,iy,iz,3) = 0.d0
        elseif (R<=R1) then
          ! inner cylinder
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) = +omega*z
          vort(:,iy,iz,3) = -omega*y
        endif
      enddo
    enddo
    call fft3 ( uk,vort )
    
  case ("taylor_green_2d")
    !--------------------------------------------------
    ! taylor green vortices
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy
        z = dble(iz)*dz
        vort(:,iy,iz,2) = dsin( y ) * dcos( z )
        vort(:,iy,iz,3) =-dcos( y ) * dsin( z )
      enddo
    enddo
    call fft3 ( uk,vort )
    
  case("infile")
     !--------------------------------------------------
     ! read HDF5 files
     !--------------------------------------------------  
     if (mpirank==0) write (*,*) "*** inicond: reading infiles"
     call Read_Single_File ( file_ux, vort(:,:,:,1) )
     call Read_Single_File ( file_uy, vort(:,:,:,2) )
     call Read_Single_File ( file_uz, vort(:,:,:,3) )
     call fft3 ( uk,vort )
     if (mpirank==0) write (*,*) "*** done reading infiles"
  case("VortexRing")
     !--------------------------------------------------
     ! Vortex ring
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: vortex ring initial condition"
     r00=yl/8.d0
     a  =0.4131d0 * r00 
     a  =0.82d0 * r00 
     gamma0=12.0d0
     x00=0.5d0 * xl

     ! define vorticity in phy space
     do iz=ra(3), rb(3)
        z=dble(iz) * zl / dble(nz)
        do iy=ra(2), rb(2)
           y=dble(iy) * yl / dble(ny)
           do ix=ra(1), rb(1)
              x=dble(ix) * xl / dble(nx)
              r=dsqrt( (y-0.5d0*yl)**2 + (z-0.5d0*zl)**2 )
              omega=(gamma0 / (pi*a**2))&
                   *dexp(-( (x-x00)**2 + (r-r00)**2 )/(a**2) )

              if ( dabs(r)> 1.0d-12) then
                 vort (ix,iy,iz,2)=-omega * ( (z-0.5d0*zl)/r) ! sin(theta)
                 vort (ix,iy,iz,3)= omega * ( (y-0.5d0*yl)/r) ! cos(theta)
              else
                 vort (ix,iy,iz,2)=0.d0
                 vort (ix,iy,iz,3)=0.d0
              endif

           end do
        end do
     end do
     call Vorticity2Velocity(uk, nlk(:,:,:,:,0), vort)
     
     call set_mean_flow(uk,time)
     
  case ("vortex")
     if (mpirank==0) write (*,*) "*** inicond: vortex ring initial condition"
     r00=yl/8.d0
     a  =0.4131d0 * r00 
     a  =0.82d0 * r00 
     gamma0=12.0d0
     x00=0.5d0 * xl

     ! define vorticity in phy space
     do iz=ra(3), rb(3)
        z=dble(iz) * zl / dble(nz)
        do iy=ra(2), rb(2)
           y=dble(iy) * yl / dble(ny)
           do ix=ra(1), rb(1)
              x=dble(ix) * xl / dble(nx)
              
              r=dsqrt( (x-0.5d0*xl)**2+ (y-0.5d0*yl)**2 + (z-0.5d0*zl)**2 )
              
              omega = 10.d0*exp( -(r / (0.1*zl))**2 )

              vort (ix,iy,iz,3)=omega
           enddo
        enddo
     enddo
     
     call Vorticity2Velocity(uk, nlk(:,:,:,:,0), vort)
     
     call set_mean_flow(uk,time)
  
  case("turbulence")
     !--------------------------------------------------
     ! random vorticity
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: turbulence (random vorticity) initial condition"
     call random_seed()
     call create_mask( 0.d0, Insect, beams )
     do iz=ra(3), rb(3)
        do iy=ra(2), rb(2)
           do ix=ra(1), rb(1)
              vort (ix,iy,iz,1)=50.d0*(2.0d0*rand_nbr() - 1.d0)
              vort (ix,iy,iz,2)=50.d0*(2.0d0*rand_nbr() - 1.d0)
              vort (ix,iy,iz,3)=50.d0*(2.0d0*rand_nbr() - 1.d0)
           end do
        end do
     end do
     
!      call Read_Single_File ('vorx_inicond.h5', vort(:,:,:,1))
!      call Read_Single_File ('vory_inicond.h5', vort(:,:,:,2))
!      call Read_Single_File ('vorz_inicond.h5', vort(:,:,:,3))
     
     call cal_vis( 1.0e-2/nu, explin(:,:,:,1))
     call fft3( inx=vort, outk=nlk(:,:,:,:,0) )
     nlk(:,:,:,1,0)=nlk(:,:,:,1,0)*explin(:,:,:,1)
     nlk(:,:,:,2,0)=nlk(:,:,:,2,0)*explin(:,:,:,1)
     nlk(:,:,:,3,0)=nlk(:,:,:,3,0)*explin(:,:,:,1)
     call ifft3( ink=nlk(:,:,:,:,0), outx=vort )
     call Vorticity2Velocity (uk, nlk(:,:,:,:,0), vort)

  case("MeanFlow")
     !--------------------------------------------------
     ! mean flow only
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: mean flow"
     uk=dcmplx(0.0d0,0.0d0)
     ! note this inicond also works without meanflow forcing, it is then
     ! really just an inicond
     
     ! forcing = zeroth Fourier mode only
     if ( (ca(1) == 0) .and. (ca(2) == 0) .and. (ca(3) == 0) ) then
       uk(0, 0, 0,1) = Uxmean
       uk(0, 0, 0,2) = Uymean
       uk(0, 0, 0,3) = Uzmean
     endif

  case("quiescent")
     !--------------------------------------------------
     ! fluid at rest
     !--------------------------------------------------  
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     uk=dcmplx(0.0d0,0.0d0)

  case default
     if(inicond(1:8) == "backup::") then
        !--------------------------------------------------
        ! read from backup
        !--------------------------------------------------  
        if (mpirank==0) write (*,*) "*** inicond: retaking backup " // &
             inicond(9:len(inicond))
        call Read_Runtime_Backup(inicond(9:len(inicond)),time,dt0,dt1,n1,it,uk,&
             nlk,explin,vort(:,:,:,1))
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
    call surface_interpolation_testing( time, beams(1), press )
    call init_beams( beams )
  endif
  
  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the passive scalar(s)
  !-----------------------------------------------------------------------------
  if ((use_passive_scalar==1).and.(index(inicond,"backup::")==0)) then
    ! only if not resuming a backup
    call init_passive_scalar(uk(:,:,:,4),vort,workc(:,:,:,1),Insect,beams)
  endif
  
  !-----------------------------------------------------------------------------
  ! when computing running time avg, initialize
  !-----------------------------------------------------------------------------
  if (time_avg=="yes") then
    uk_avg(:,:,:,1:3) = uk(:,:,:,1:3)
    if (inicond(1:8) == "backup::") then
      if (mpirank==0) write(*,*) "Resuming backup and we are computing time avg"
      if (mpirank==0) write(*,*) "trying to load old avg  uavgx_0000.h5"
      call check_file_exists( "uavgx_0000.h5" )
      call check_file_exists( "uavgy_0000.h5" )
      call check_file_exists( "uavgz_0000.h5" )
      
      call Read_Single_File ( "uavgx_0000.h5", vort(:,:,:,1) )
      call fft ( inx=vort(:,:,:,1) , outk=uk_avg(:,:,:,1) )
      
      call Read_Single_File ( "uavgy_0000.h5", vort(:,:,:,2) )
      call fft ( inx=vort(:,:,:,2) , outk=uk_avg(:,:,:,2) )
      
      call Read_Single_File ( "uavgz_0000.h5", vort(:,:,:,3) )
      call fft ( inx=vort(:,:,:,3) , outk=uk_avg(:,:,:,3) )
    endif
  endif
  
end subroutine init_fields_fsi


! Computes the divergence-free velocity in Fourier space u given vort
! in physical space.  work is a work array
subroutine Vorticity2Velocity(uk,work,vort)
  use mpi
  use fsi_vars
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
  call fft3(inx=vort,outk=work)
  
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
          work(iz,iy,ix,1)=+work(iz,iy,ix,1) / k_abs_2
          work(iz,iy,ix,2)=+work(iz,iy,ix,2) / k_abs_2
          work(iz,iy,ix,3)=+work(iz,iy,ix,3) / k_abs_2
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
