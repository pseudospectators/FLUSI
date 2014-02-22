! Set initial conditions for fsi code.
subroutine init_fields_fsi(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  use mpi_header
  use fsi_vars
  implicit none

  integer,intent (inout) :: n1,it
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex (kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd),&
       intent(inout) :: uk
  complex (kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1),&
       intent(inout) :: work_nlk
  real (kind=pr),dimension(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd),&
       intent(inout) :: vort
  real (kind=pr),dimension(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)),&
       intent(inout) :: explin
  integer :: ix,iy,iz
  real (kind=pr) :: x,y,z,r,a,gamma0,x00,r00,omega

  ! Assign zero values
  time=0.0d0
  dt1=0.0d0
  uk=dcmplx(0.0d0,0.0d0)
  work_nlk=dcmplx(0.0d0,0.0d0)
  explin=0.0
  it=0
  vort=0.0d0

  select case(inicond)
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
     call Vorticity2Velocity(uk, work_nlk(:,:,:,:,0), vort)

  case("turbulence")
     !--------------------------------------------------
     ! random vorticity
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: turbulence (random vorticity) initial condition"
     call random_seed()
     do iz=ra(3), rb(3)
        do iy=ra(2), rb(2)
           do ix=ra(1), rb(1)
              call RANDOM_NUMBER(r)
              vort (ix,iy,iz,1)=50.d0*(2.0d0*r - 1.d0) 
              !* (1.d0-eps*mask(ix,iy,iz))
              call RANDOM_NUMBER(r)
              vort (ix,iy,iz,2)=50.d0*(2.0d0*r - 1.d0)
              !* (1.d0-eps*mask(ix,iy,iz))
              call RANDOM_NUMBER(r)
              vort (ix,iy,iz,3)=50.d0*(2.0d0*r - 1.d0)
              !* (1.d0-eps*mask(ix,iy,iz))
           end do
        end do
     end do
     call Vorticity2Velocity (uk, work_nlk(:,:,:,:,0), vort)

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
             work_nlk,explin,vort(:,:,:,1))
     else
        if (mpirank==0) write (*,*) inicond
        if (mpirank==0) write (*,*) '??? ERROR: Invalid initial condition'
        call abort
     endif
  end select
end subroutine init_fields_fsi


! Computes the divergence-free velocity in Fourier space u given vort
! in physical space.  work is a work array
subroutine Vorticity2Velocity(uk,work,vort)
  use mpi_header
  use fsi_vars
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
          work(iz,ix,iy,1)=-work(iz,ix,iy,1) / k_abs_2
          work(iz,ix,iy,2)=-work(iz,ix,iy,2) / k_abs_2
          work(iz,ix,iy,3)=-work(iz,ix,iy,3) / k_abs_2
        else
          work(iz,ix,iy,1)=dcmplx(0.d0,0.d0)
          work(iz,ix,iy,2)=dcmplx(0.d0,0.d0)
          work(iz,ix,iy,3)=dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo

  !-----------------------------------------------
  !-- compute velocity as curl of streamfunction
  !-----------------------------------------------
  do iy=ca(3), cb(3) ! ky : 0..ny/2-1 ,then, -ny/2..-1     
    ky=scaley*dble(modulo(iy+ny/2,ny)-ny/2)     
    do ix=ca(2), cb(2) ! kx : 0..nx/2
      kx=scalex*dble(ix)
      do iz=ca(1),cb(1) ! kz : 0..nz/2-1 ,then, -nz/2..-1
        kz=scalez*dble(modulo(iz+nz/2,nz)-nz/2)
        uk(iz,ix,iy,1)=im*(ky*work(iz,ix,iy,3) - kz*work(iz,ix,iy,2))
        uk(iz,ix,iy,2)=im*(kz*work(iz,ix,iy,1) - kx*work(iz,ix,iy,3))
        uk(iz,ix,iy,3)=im*(kx*work(iz,ix,iy,2) - ky*work(iz,ix,iy,1))
      enddo
    enddo
  enddo
end subroutine Vorticity2Velocity
