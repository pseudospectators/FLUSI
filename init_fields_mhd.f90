! Initialize fields for mhd simulations
subroutine init_fields_mhd(n1,time,it,dt0,dt1,ubk,nlk,wj,explin)
  use mpi_header
  use fsi_vars
  implicit none

  integer,intent (inout) :: n1,it
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex (kind=pr),intent (out):: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent (out)::&
       nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real(kind=pr),intent (inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  
  ! Assign zero values
  time=0.0d0
  dt1=0.0d0
  ubk=dcmplx(0.0d0,0.0d0)
  nlk=dcmplx(0.0d0,0.0d0)
  explin=0.0
  it=0
  wj=0.0d0

  select case(inicond)
  case("quiescent")
     ubk=dcmplx(0.0d0,0.0d0)
  case("constant")
     call init_const(ubk,wj)
  case("orszagtang")
     call init_orszagtang(ubk,wj)
  case("smc")
     call init_smc(ubk,wj)
  case default
     if(inicond(1:8) == "backup::") then
        call Read_Runtime_Backup(inicond(9:len(inicond)),&
             time,dt0,dt1,n1,it,ubk,nlk,explin,wj(:,:,:,1))
     else
        if (mpirank==0) then
           write (*,*) inicond
           write (*,*) '??? ERROR: Invalid initial condition'
        endif
        call abort
     endif
  end select

  ! Ensure that initial conditions are divergence-free by performing a
  ! Helmholtz decomposition:
  call div_field_nul(ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3))
  call div_field_nul(ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6))
end subroutine init_fields_mhd


! The Orszag-Tang initial conditions for mhd.
subroutine init_orszagtang(ubk,wj)
  use mpi_header
  use mhd_vars
  implicit none

  complex(kind=pr),intent(inout):: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent (inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr) :: beta
  integer :: ix,iy,iz,i
  real(kind=pr) :: x,y,z  

  beta=0.8d0

  do ix=ra(1),rb(1)
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        do iz=ra(3),rb(3)
           z=zl*(dble(iz)/dble(nz) -0.5d0)

           wj(ix,iy,iz,1)=-2.d0*dsin(y)
           wj(ix,iy,iz,2)=2.d0*dsin(x)
           wj(ix,iy,iz,3)=0.d0

           wj(ix,iy,iz,4)=beta*(-2.d0*dsin(2.d0*y) + dsin(z))
           wj(ix,iy,iz,5)=beta*(2.d0*dsin(x) + dsin(z))
           wj(ix,iy,iz,6)=beta*(dsin(x) + dsin(y))
        enddo
     enddo
  enddo
  
  do i=1,nd
     call fft(ubk(:,:,:,i),wj(:,:,:,i))
  enddo
end subroutine init_orszagtang


! Constant initial conditions for MHD
subroutine init_const(ubk,wj)
  use mpi_header
  use mhd_vars
  implicit none

  complex(kind=pr),intent(inout):: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent (inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  integer :: i

  wj=1.d0
  
  do i=1,nd
     call fft(ubk(:,:,:,i),wj(:,:,:,i))
  enddo
end subroutine init_const


! The Sean-Montgomery-Chen initial conditions
subroutine init_smc(ubk,ub)
  use mpi_header
  use mhd_vars
  implicit none

  complex(kind=pr),intent(inout):: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent (inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  integer :: mpicode
  integer :: i, ix,iy,iz
  real(kind=pr) :: pekin,pfluidarea,fluidarea,ekin
  real(kind=pr) :: ux,uy,uz
  real(kind=pr) :: enorm
  real(kind=pr) :: x,y,r
  real (kind=pr) :: a,b,c,d,k1,k2,h

  ! Create random initial conditions for the velocity:
  call randgen3d(ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3),ub(:,:,:,1))
  
  do i=1,3
     call ifft(ub(:,:,:,i),ubk(:,:,:,i))
  enddo

  ! Compute energy and fluid area (for normalization later).
  pekin=0.d0
  pfluidarea=0.d0
  do ix=ra(1),rb(1)
     do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
           if(mask(ix,iy,iz) == 0.d0) then
              ux=ub(ix,iy,iz,1)
              uy=ub(ix,iy,iz,2)
              uz=ub(ix,iy,iz,3)

              pekin = pekin + ux*ux + uy*uy + uz*uz
              pfluidarea = pfluidarea + 1.d0
           endif
        enddo
     enddo
  enddo

  pekin = 0.5*pekin
  
  call MPI_REDUCE(pekin,ekin,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(pfluidarea,fluidarea,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       0,MPI_COMM_WORLD,mpicode)
  
  if(mpirank == 0) then
     ekin = ekin/fluidarea
  endif
  call MPI_BCAST(ekin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)

  ! normalize energy to that used in Jorge's initialization file.
  enorm=dsqrt(3.1017126d-07/ekin)
  do i=1,3
     ub(:,:,:,i)=ub(:,:,:,i)*enorm
     call fft(ubk(:,:,:,i),ub(:,:,:,i))
  enddo
  
  ! Set up the magnetic field

  ub(:,:,:,4)=0.d0
  ub(:,:,:,5)=0.d0

  ub(:,:,:,6)=B0 ! bz set to B0 everywhere

  k1=Bc*r2/r1
  k2=Bc/r1
  A=(2.d0*k1 -k2*(r2-r3))/(r3*r3*r3 -3.d0*r2*r3*r3 +3.d0*r2*r2*r3 -r2*r2*r2)
  B=(k2 -3.d0*A*(r2*r2 -r3*r3))/(2.d0*r2 -2.d0*r3)
  C=-3.d0*A*r3*r3 -2.d0*B*r3
  D=2.d0*A*r3*r3*r3 +B*r3*r3

  do ix=ra(1),rb(1)
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        
        r=dsqrt(x*x + y*y)

        if(r < r2) then
           do iz=ra(3),rb(3)
              ub(ix,iy,iz,4)=-y*Bc/R1
              ub(ix,iy,iz,5)= x*Bc/R1
           enddo
        endif

        if(r >= r2 .and. r <= r3) then
           h=(A*r*r*r +B*r*r +C*r +D)
           do iz=ra(3),rb(3)
              ub(ix,iy,iz,4)= h*y/r
              ub(ix,iy,iz,5)=-h*x/r
           enddo
        endif

     enddo
  enddo

  do i=3,nd
     call fft(ubk(:,:,:,i),ub(:,:,:,i))
  enddo
end subroutine init_smc


! Create a random initial condition based with zero divergence in fk1,
! fk2, fk3. w is a 3D real work array.
! FIXME: please add more documentation.
subroutine randgen3d(fk1,fk2,fk3,w)
  use mpi_header
  use vars
  implicit none

  complex(kind=pr),intent(inout):: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  real(kind=pr),intent(inout) :: w(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr) :: ek
  real(kind=pr) :: e1,e1x,e1y,e1z
  real(kind=pr) :: e2,e2x,e2y,e2z
  real(kind=pr) :: k0
  real(kind=pr) :: rand1,rand2
  real(kind=pr) :: kx,ky,kz,k
  real(kind=pr) :: psi1,beta
  integer :: ix,iy,iz

  fk1=0.d0
  fk2=0.d0
  fk3=0.d0

  k0=3.d0*sqrt(2.d0)*pi/4.d0

  !-- generate 3D gaussian field 
  do iz=ca(1),cb(1)
     kz=scalez*dble(modulo(iz+nz/2,nz) -nz/2)
     do ix=ca(2),cb(2)
        kx=scalex*dble(ix)
        do iy=ca(3),cb(3)
           ky=scaley*dble(modulo(iy+ny/2,ny) -ny/2)
           k=dsqrt(kx*kx +ky*ky +kz*kz)

           !-- ek = sqrt(Ek/4pik^2), Ek imposed spectrum
           if(k /= 0.d0) then
              ek=dsqrt(1.d0/(4.d0*pi*k*k))
           else
              ek=0.d0
           endif

           e1=dsqrt(kx*kx+ky*ky)
           e1x=ky/e1
           e1y=-kx/e1
           e1z=0.d0
           e2=k*dsqrt(kx*kx+ky*ky)
           e2x=kx*kz/e2
           e2y=ky*kz/e2
           e2z=-(kx*kx+ky*ky)/e2

           if(kx.eq.0 .and. ky.eq.0) then
              fk1(iz,ix,iy)=dcmplx(0.5d0*ek,0.5d0*ek)
              fk2(iz,ix,iy)=dcmplx(0.5d0*ek,0.5d0*ek)
              fk3(iz,ix,iy)=dcmplx(0.d0,0.d0)
              if(kz.eq.0) then
                 fk1(iz,ix,iy)=dcmplx(0.d0,0.d0)
                 fk2(iz,ix,iy)=dcmplx(0.d0,0.d0)
                 fk3(iz,ix,iy)=dcmplx(0.d0,0.d0)
              endif
           else
              call random_number(rand1)
              call random_number(rand2)
              psi1=2.*pi*rand1
              beta=2.*pi*rand2
              fk1(iz,ix,iy)=dcmplx(&
                   ek*dcos(psi1)*(dcos(beta)*e1x+dsin(beta)*e2x),&
                   ek*dsin(psi1)*(dcos(beta)*e1x+dsin(beta)*e2x)&
                   )
              fk2(iz,ix,iy)=dcmplx(&
                   ek*dcos(psi1)*(dcos(beta)*e1y+dsin(beta)*e2y),&
                   ek*dsin(psi1)*(dcos(beta)*e1y+dsin(beta)*e2y)&
                   )
              fk3(iz,ix,iy)=dcmplx(&
                   ek*dcos(psi1)*(dsin(beta)*e2z),&
                   ek*dsin(psi1)*(dsin(beta)*e2z)&
                   )
           endif
        enddo
     enddo
  enddo

  ! Enforce Hermitian symmetry the lazy way:
  call ifft(w,fk1)
  call fft(fk1,w)

  call ifft(w,fk2)
  call fft(fk2,w)

  call ifft(w,fk3)
  call fft(fk3,w)
end subroutine randgen3d
