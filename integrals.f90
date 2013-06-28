! Wrapper for writing integral quantities to file
subroutine write_integrals(time,uk,u,vort,nlk)
  use mpi_header
  use vars
  implicit none

  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr), intent(in) :: time

  select case(method(1:3))
  case("fsi")
     call write_integrals_fsi(time,uk,u,vort,nlk)
  case("mhd")
     call write_integrals_mhd(time,uk,u,vort,nlk)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in write_integrals"
     call abort
  end select
end subroutine write_integrals


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time,uk,u,vort,nlk)
  use mpi_header
  use fsi_vars
  implicit none

  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr), intent(in) :: time
  
  ! FIXME: compute integral quantities

  if(mpirank == 0) then
     open(14,file='drag_data',status='unknown',position='append')
     write(14,'(7(es12.4,1x))')  time,GlobalIntegrals%Ekin,&
          GlobalIntegrals%Dissip, GlobalIntegrals%Force(1),&
          GlobalIntegrals%Force(2),GlobalIntegrals%Force(3),&
       GlobalIntegrals%Volume
     close(14)
  endif
end subroutine write_integrals_fsi


! mhd version of writing integral quantities to disk
subroutine write_integrals_mhd(time,ubk,ub,wj,nlk)
  use mpi_header
  use mhd_vars
  implicit none

  complex (kind=pr),intent(inout)::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(inout) ::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real(kind=pr), intent(in) :: time
  integer :: ix,iy,iz,i,mpicode
  type(MHDIntegrals) LI ! local copy of integrals for mhd
  ! Local loop variables
  real(kind=pr) :: w1,w2,w3,j1,j2,j3
  real(kind=pr) :: u1,u2,u3,b1,b2,b3
  
  ! NB: integral quantities are defined in mhd_vars
  LI%Ekin=0.d0
  LI%Emag=0.d0

  LI%jmax=0.d0
  LI%jx=0.d0
  LI%jy=0.d0
  LI%jz=0.d0
  
  ! Compute u and B to physical space
  do i=1,nd
     call ifft(ub(:,:,:,i),ubk(:,:,:,i))
  enddo
  
  ! Compute the vorticity and store the result in the first three 3D
  ! arrays of nlk.
  call curl(nlk(:,:,:,1),nlk(:,:,:,2),nlk(:,:,:,3),&
       ubk(:,:,:,1),ubk(:,:,:,2),ubk(:,:,:,3))

  ! Compute the current density and store the result in the last three
  ! 3D arrays of nlk.
  call curl(nlk(:,:,:,4),nlk(:,:,:,5),nlk(:,:,:,6),&
       ubk(:,:,:,4),ubk(:,:,:,5),ubk(:,:,:,6))

  ! Transform vorcitity and current density to physical space, store
  ! in wj
  do i=1,nd
     call ifft(wj(:,:,:,i),nlk(:,:,:,i))
  enddo

  ! FIXME: TODO: compute more integral quantities

  do iy=ra(3),rb(3)
     do ix=ra(2),rb(2)
        do iz=ra(1),rb(1)
           ! Loop-local variables for velocity and magnetic field:
           u1=ub(iz,ix,iy,1)
           u2=ub(iz,ix,iy,2)
           u3=ub(iz,ix,iy,3)
           b1=ub(iz,ix,iy,4)
           b2=ub(iz,ix,iy,5)
           b3=ub(iz,ix,iy,6)

           ! Loop-local variables for vorticity and current density:
           w1=wj(iz,ix,iy,1)
           w2=wj(iz,ix,iy,2)
           w3=wj(iz,ix,iy,3)
           j1=wj(iz,ix,iy,4)
           j2=wj(iz,ix,iy,5)
           j3=wj(iz,ix,iy,6)
           
           ! Ekin, Emag
           LI%Ekin=LI%Ekin + u1*u1 + u2*u2 + u3*u3
           LI%Emag=LI%Emag + b1*b1 + b2*b2 + b2*b2

           ! jmax, jx, jy, jz
           LI%jmax=max(LI%jmax,j1+j1 + j2*j2 + j3*j3)
           LI%jx=LI%jx + j1
           LI%jx=LI%jx + j2
           LI%jx=LI%jx + j3
        enddo
     enddo
  enddo
  LI%Ekin=LI%Ekin*dx*dy*dz
  LI%Emag=LI%Emag*dx*dy*dz

  ! Ekin, Emag
  call MPI_REDUCE(LI%Ekin,integrals%Ekin,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(LI%Emag,integrals%Emag,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  
  ! jmax, jx, jy, jz
  call MPI_REDUCE(LI%jmax,integrals%jmax,&
       1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
       MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(LI%jx,integrals%jx,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(LI%jy,integrals%jy,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(LI%jz,integrals%jz,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)

  if(mpirank == 0) then
     open(14,file='evt',status='unknown',position='append')
     write(14,97) time,integrals%Ekin,integrals%Emag
     close(14)
   open(14,file='jvt',status='unknown',position='append')
     write(14,97) time,integrals%jmax,integrals%jx,integrals%jy,integrals%jz
     close(14)
97   format(1X,9(E14.7,' ')) ! Why must Fortran require this nonsense?
  endif
end subroutine write_integrals_mhd
