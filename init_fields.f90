! Wrapper for init_fields
subroutine init_fields(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  use mpi_header
  use vars
  implicit none

  integer,intent (inout) :: n1,it
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex (kind=pr), intent(inout) :: &
       uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex (kind=pr),intent(inout) :: &
       work_nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1)
  real (kind=pr), intent(inout) ::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)

  if (mpirank == 0) write(*,*) "Initializating fields:"
 
  select case(method)
  case("fsi") 
     call init_fields_fsi(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  case("mhd")
     call init_fields_mhd(n1,time,it,dt0,dt1,uk,work_nlk,vort,explin)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in init_fields."
     call abort
  end select
end subroutine init_fields


! Create a randomized, divergence-free field in fk1,fk3,fk3 which is
! normalized to have the given energy in the fluid domain.
subroutine perturbation(fk1,fk2,fk3,f1,f2,f3,energy)
  use mpi_header
  use vars
  implicit none
  
  complex(kind=pr),intent(inout):: fk1(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: fk2(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  complex(kind=pr),intent(inout):: fk3(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

  real(kind=pr),intent (inout) :: f1(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent (inout) :: f2(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent (inout) :: f3(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  real(kind=pr),intent (in) :: energy

  integer :: mpicode
  integer :: ix,iy,iz
  real(kind=pr) :: ux,uy,uz
  real(kind=pr) :: pekin,pfluidarea,fluidarea,ekin
  real(kind=pr) :: enorm

  ! Create random initial conditions for the velocity:
  call randgen3d(fk1,fk2,fk3,f1)

  ! Compute energy and fluid area (for normalization later).
  pekin=0.d0
  pfluidarea=0.d0
  do ix=ra(1),rb(1)
     do iy=ra(2),rb(2)
        do iz=ra(3),rb(3)
           if(mask(ix,iy,iz) == 0.d0) then
              ux=f1(ix,iy,iz)
              uy=f2(ix,iy,iz)
              uz=f3(ix,iy,iz)

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

  ! Normalize the field to the desired energy:
  enorm=dsqrt(energy/ekin)
  f1=f1*enorm
  f2=f2*enorm  
  f3=f3*enorm

  call fft(fk1,f1)
  call fft(fk2,f2)
  call fft(fk3,f3)
end subroutine perturbation


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


! Initialize the first three fields of ubk to the steady-state
! solution for Taylor-Couette flow.
subroutine init_taylorcouette_u(ubk,ub)
  use mpi_header
  use vars
  implicit none

  complex(kind=pr),intent(inout):: ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr),intent (inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  integer :: ix,iy,iz,i
  real(kind=pr) :: x,y,r
  real(kind=pr) :: A,B,F

    if(mpirank == 0) then
     if(r1 >= r2) then
        write (*,*) "r1 >= r2 is not allowed in Taylor-Coette flow; stopping."
        stop
     endif
  endif
  
  A=-omega1*r1*r1/(r2*r2 - r1*r1)
  B=omega1*r1*r1*r2*r2/(r2*r2 - r1*r1)

  do ix=ra(1),rb(1)
     x=xl*(dble(ix)/dble(nx) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)

        r=dsqrt(x*x +y*y)

        if(r <= R1) then
           do iz=ra(3),rb(3)
              ub(ix,iy,iz,1)=-omega1*y
              ub(ix,iy,iz,2)=omega1*x
           enddo
        endif

        if(r > R1 .and. r < R2) then
           do iz=ra(3),rb(3)
              F=A*r+B/r
              ub(ix,iy,iz,1)=-F*y/r
              ub(ix,iy,iz,2)=F*x/r
           enddo
        endif
        
        if(r >= R2) then
           do iz=ra(3),rb(3)
              ! NB: We assume that the outer wall is not moving.
              us(ix,iy,iz,1)=0.d0
              us(ix,iy,iz,2)=0.d0
           enddo
        endif

     enddo
  enddo

  us(:,:,:,3)=0.d0
  
  ! Transform velocity field to Fourier space:
  do i=1,3
     call fft(ubk(:,:,:,i),ub(:,:,:,i))
  enddo
  
end subroutine init_taylorcouette_u
