! Wrapper for writing integral quantities to file
subroutine write_integrals(time,uk,u,vort,work)
  use mpi_header
  use vars
  implicit none

  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real (kind=pr),intent(inout) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr), intent(in) :: time

  select case(method(1:3))
  case("fsi")
     call write_integrals_fsi(time,uk,u,vort,work)
  case("mhd")
     call write_integrals_mhd(time,uk,u,vort,work)
  case default
     if (mpirank == 0) write(*,*) "Error! Unkonwn method in write_integrals"
     call abort
  end select
end subroutine write_integrals


! fsi version of writing integral quantities to disk
subroutine write_integrals_fsi(time,uk,u,vort,work)
  use mpi_header
  use fsi_vars
  implicit none

  complex (kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
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
subroutine write_integrals_mhd(time,ubk,ub,wj,work)
  use mpi_header
  use mhd_vars
  implicit none

  complex (kind=pr),intent(inout)::ubk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd)
  real (kind=pr),intent(inout) :: ub(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),intent(inout) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr), intent(in) :: time
  integer :: ix,iy,iz,i,mpicode
  real(kind=pr) :: kx,ky,kz,k2
  type(MHDIntegrals) local
  
  ! NB: integral quantities are defined in mhd_vars
  local%Ekin=0.d0
  local%Bkin=0.d0
  
  ! Compute u and B to physical space
  do i=1,nd
     call ifft(ub(:,:,:,i),ubk(:,:,:,i))
  enddo
  
  ! FIXME: TODO: compute more integral quantities  

  do iy=ra(3),rb(3)
     do ix=ra(2),rb(2)
        do iz=ra(1),rb(1)
           local%Ekin=local%Ekin &
                + ub(iz,ix,iy,1)*ub(iz,ix,iy,1)&
                + ub(iz,ix,iy,2)*ub(iz,ix,iy,2)&
                + ub(iz,ix,iy,3)*ub(iz,ix,iy,3)
           local%Bkin=local%Bkin &
                + ub(iz,ix,iy,4)*ub(iz,ix,iy,4)&
                + ub(iz,ix,iy,5)*ub(iz,ix,iy,5)&
                + ub(iz,ix,iy,6)*ub(iz,ix,iy,6)
        enddo
     enddo
  enddo
  local%Ekin=local%Ekin*dx*dy*dz
  local%Bkin=local%Bkin*dx*dy*dz
  
  call MPI_REDUCE(local%Ekin,integrals%Ekin,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  call MPI_REDUCE(local%Bkin,integrals%Bkin,&
       1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
       MPI_COMM_WORLD,mpicode)
  
  if(mpirank == 0) then
     open(14,file='evt',status='unknown',position='append')
     write(14,97) time,integrals%Ekin,integrals%Bkin
     close(14)
97   format(1X,9(E14.7,' ')) ! Why must Fortran require this nonsense?
  endif
end subroutine write_integrals_mhd
