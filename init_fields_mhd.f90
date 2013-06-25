! Initialize fields for mhd simulations
subroutine init_fields_mhd(n1,time,it,dt0,dt1,ubk,nlk,wj,explin)
  use mpi_header
  use fsi_vars
  implicit none

  integer,intent (inout) :: n1,it
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex (kind=pr),dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd),&
       intent (out) :: ubk
  complex (kind=pr),dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd,0:1),&
       intent (out) :: nlk
  real(kind=pr),intent (inout) :: wj(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real (kind=pr),dimension (ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf),&
       intent(inout) :: explin
  
  ! Assign zero values
  time=0.0d0
  dt1=0.0d0
  ubk=dcmplx(0.0d0,0.0d0)
  nlk=dcmplx(0.0d0,0.0d0)
  explin=0.0
  it=0
  wj=0.0d0

  ! TODO: add more initial conditions
  select case(inicond)
  case("quiescent")
     ! Set the velocity and magnetic fields to zero.
     !--------------------------------------------------  
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     ubk=dcmplx(0.0d0,0.0d0)

  case("orszagtang")
     call init_orszagtang(ubk,wj)
  case default
     if(inicond(1:8) == "backup::") then
        !--------------------------------------------------
        ! read from backup
        !--------------------------------------------------  
        if (mpirank==0) write (*,*) "*** inicond: retaking backup " // &
             inicond(9:len(inicond))
        call Read_Runtime_Backup(inicond(9:len(inicond)),&
             time,dt0,dt1,n1,it,ubk,nlk,explin,wj(:,:,:,1))
     else
        if (mpirank==0) write (*,*) inicond
        if (mpirank==0) write (*,*) '??? ERROR: Invalid initial condition'
        call abort
     endif
  end select

end subroutine init_fields_mhd

! The Orszag-Tang initial conditions for MHD.
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

  do iz=ra(3),rb(3)
     z=zl*(dble(iz)/dble(nz) -0.5d0)
     do iy=ra(2),rb(2)
        y=yl*(dble(iy)/dble(ny) -0.5d0)
        do ix=ra(1),rb(1)
           x=xl*(dble(ix)/dble(nx) -0.5d0)
           wj(ix,iy,iz,1)=-2.d0*dsin(y)
           wj(ix,iy,iz,2)=2.d0*dsin(x)
           wj(ix,iy,iz,3)=0.d0
           wj(ix,iy,iz,4)=beta*(-2.*dsin(2.*y) + dsin(z))
           wj(ix,iy,iz,5)=beta*(2.*dsin(x) + dsin(z))
           wj(ix,iy,iz,6)=beta*(dsin(x) + dsin(y))
        enddo
     enddo
  enddo
  
  do i=1,nd
     call fft(ubk(:,:,:,i),wj(:,:,:,i))
  enddo
     
end subroutine init_orszagtang

