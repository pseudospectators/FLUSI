subroutine lagrangian_acc(uk,u,acc)
!subroutine lagrangian_acc(uk,u,vort,acc,lap)    
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  implicit none

  real(kind=pr),intent(in)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  !real(kind=pr),intent(in)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) 
  complex(kind=pr),intent(in)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  
  real(kind=pr),intent(out)::acc(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  !real(kind=pr),intent(out)::lap(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd) 
 
  select case(method)
  case("fsi") 
     call lagrangian_acc_flusi(uk,u,acc)
  !case("mhd")
     !call lagrangian_acc_mhd(time,it,dt0,dt1,n0,n1,uk,nlk,vort,explin)   
  case default
     call abort(123789,"Error! Unkonwn method in lagrangian acceleration.")
  end select

end subroutine lagrangian_acc


subroutine lagrangian_acc_flusi(uk,u,acc) 
  use mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use penalization ! mask array etc
  use module_insects
  implicit none

  real(kind=pr),intent(in) :: u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  complex(kind=pr),intent(in) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nd) 

  real(kind=pr),intent(out):: acc(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  !real(kind=pr),intent(out):: lap(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  
  real(kind=pr) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr) :: work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))

  complex(kind=pr) ::acclk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr) ::pk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))   
  complex(kind=pr) ::pkyz(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr) :: nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  complex(kind=pr) ::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  
  integer :: ix,iy,iz, i
  real(kind=pr) :: kx,ky,kz,k2
  complex(kind=pr) :: imag

  type(diptera) :: Insect

  imag = dcmplx(0.d0,1.d0)
 
  !---------------------------------------------------------------------------
  !------- Compute the Lagrangian acceleration for hasegawa wakatani ---------
  !--- acc_lag= -\nabla P + nu \nabla^2 u -chi(u-us) ---
  !---------------------------------------------------------------------------    
  
  ! Compute \nabla P in Fourier Space
  call cal_nlk_fsi (0.d0,0,nlk,uk,u,vort,work,workc, Insect)
  call pressure( nlk,pk )
  ! total pressure in x-space
  call ifft( ink=pk, outx=work )
  ! get actuall pressure (we're in the rotational formulation)
  work = work - 0.5d0*( u(:,:,:,1)**2 + u(:,:,:,2)**2 + u(:,:,:,3)**2 )
  call fft(inx=work,outk=workc(:,:,:,1) )

  do iz=ca(1),cb(1)
     kz=wave_z(iz)
     do iy=ca(2),cb(2)
        ky=wave_y(iy)
        do ix=ca(3),cb(3)
          kx=wave_x(ix)
          k2= kx*kx + ky*ky + kz*kz
          pkyz(iz,iy,ix,1)= kx*imag*workc(iz,iy,ix,1)
          pkyz(iz,iy,ix,2)= ky*imag*workc(iz,iy,ix,1)
          pkyz(iz,iy,ix,3)= kz*imag*workc(iz,iy,ix,1)

          ! d2x(ux), d2y(uy) and d2z(uz)
          workc(iz,iy,ix,1) = -uk(iz,iy,ix,1)*k2
          workc(iz,iy,ix,2) = -uk(iz,iy,ix,2)*k2
          workc(iz,iy,ix,3) = -uk(iz,iy,ix,3)*k2

          acclk(iz,iy,ix,1) = -pkyz(iz,iy,ix,1) + nu*workc(iz,iy,ix,1)
          acclk(iz,iy,ix,2) = -pkyz(iz,iy,ix,2) + nu*workc(iz,iy,ix,2)
          acclk(iz,iy,ix,3) = -pkyz(iz,iy,ix,3) + nu*workc(iz,iy,ix,3)
      enddo
    enddo
  enddo
   
  ! Lagrangian acceleration in physical space
  call ifft3 ( ink=acclk(:,:,:,1:3), outx=acc(:,:,:,1:3) )
  
  acc(:,:,:,1)=acc(:,:,:,1) + mask*(u(:,:,:,1)-us(:,:,:,1))
  acc(:,:,:,2)=acc(:,:,:,2) + mask*(u(:,:,:,2)-us(:,:,:,2))
  acc(:,:,:,3)=acc(:,:,:,3) + mask*(u(:,:,:,3)-us(:,:,:,3))
  !! Lagrangian acceleration in physical space
  !call ifft(lap(:,:,:,1),workc(:,:,:,1))
  !call ifft(lap(:,:,:,2),workc(:,:,:,2))  
  !
end subroutine lagrangian_acc_flusi
