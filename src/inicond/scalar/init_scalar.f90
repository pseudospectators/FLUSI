!-------------------------------------------------------------------------------
! Passive scalar initial conditions, called from init_fields_fsi
! when resuming a backup file, this is automatically skipped (the scalar is then
! loaded from the runtime backup in init_fields_fsi-->Read_Runtime_Backup)
!
! INPUT/OUTPUT
!     uk: 3D array that will hold the passive scalar theta in F-space on output
!     work: three real valued work arrays
!     workc: one complex valued work array
!-------------------------------------------------------------------------------
subroutine init_passive_scalar(uk,work,workc,Insect,beams)
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  implicit none

  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))
  ! the workc array is not always allocated, ensure allocation before using
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3)) 
  
  
  
  
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw) !PASS AUFFFFFFF
 
 
 
 
 
 integer :: ix,iy,iz
  real (kind=pr) :: x,y,z
  type(solid),dimension(1:nBeams), intent(inout) :: beams
  type(diptera),intent(inout)::Insect 

!   select case(inicond_scalar)
!   case("cosine") 
!     !-- initialize scalar zero everywhere
!     work(:,:,:,1) = 0.d0
!     !-- set half the domain to one
!     do iz=ra(3), rb(3)
!       z = dble(iz)*dz
!       do iy=ra(2), rb(2)
!          y = dble(iy)*dy
!          do ix=ra(1), rb(1)
!            x = dble(ix)*dx - 0.5d0*xl
!            work(ix,iy,iz,1) = dcos(pi*x) + dcos(4.d0*pi*x) 
!          enddo
!       enddo
!     enddo
!   
!     call fft ( inx=work(:,:,:,1), outk=uk )
!   
!   case("right_left_discontinuous")
!     !---------------------------------------------------------------------------
!     ! one half of the domain is 1, the other is 0, divided along y-axis
!     ! no scalar inside mask (we build this here since FLUSI first loads inicond,
!     ! then creates the mask!)
!     !---------------------------------------------------------------------------
!     call create_mask( 0.d0, Insect,beams )
!     
!     !-- initialize scalar zero everywhere
!     work(:,:,:,1) = 0.d0
!     !-- set half the domain to one
!     do iz=ra(3), rb(3)
!       z = dble(iz)*dz
!       do iy=ra(2), rb(2)
!          y = dble(iy)*dy
!          do ix=ra(1), rb(1)
!            x = dble(ix)*dx
!            if (y>0.5d0*yl) then
!              work(ix,iy,iz,1) = 1.d0
!            endif
!          enddo
!       enddo
!     enddo
!     
!     !-- kill scalar inside obstacle
!     work(:,:,:,1) = work(:,:,:,1) * (1.0-mask*eps)
!     !-- transform this to F-space
!     call fft ( inx=work(:,:,:,1), outk=uk )
!     
!   case ("right_left_smooth")
!     !---------------------------------------------------------------------------
!     ! smoothed heaviside function, covering approx. half the domain
!     !---------------------------------------------------------------------------
!     call create_mask( 0.d0, Insect,beams )
!     
!     !-- initialize scalar zero everywhere
!     work(:,:,:,1) = 0.d0
!     !-- set half the domain to one
!     do iz=ra(3), rb(3)
!       z = dble(iz)*dz
!       do iy=ra(2), rb(2)
!          y = dble(iy)*dy
!          do ix=ra(1), rb(1)
!            x = dble(ix)*dx
!            call smoothstep( work(ix,iy,iz,1), abs(y-0.5*yl), 0.25*yl, 4.0*dy)
!          enddo
!       enddo
!     enddo
!     
!     !-- kill scalar inside obstacle
!     work(:,:,:,1) = work(:,:,:,1) * (1.0-mask*eps)
!     !-- transform this to F-space
!     call fft ( inx=work(:,:,:,1), outk=uk )  
!     
!   case ("right_left_smooth2")
!     !---------------------------------------------------------------------------
!     ! smoothed heaviside function, covering approx. half the domain
!     !---------------------------------------------------------------------------
!     call create_mask( 0.d0, Insect,beams )
!     
!     !-- initialize scalar zero everywhere
!     work(:,:,:,1) = 0.d0
!     !-- set half the domain to one
!     do iz=ra(3), rb(3)
!       z = dble(iz)*dz
!       do iy=ra(2), rb(2)
!          y = dble(iy)*dy
!          do ix=ra(1), rb(1)
!            x = dble(ix)*dx
!            call smoothstep( work(ix,iy,iz,1), abs(z-0.5*zl), 0.25*zl, 2.0*dz)
!          enddo
!       enddo
!     enddo
!     
!     !-- kill scalar inside obstacle
!     work(:,:,:,1) = work(:,:,:,1) * (1.0-mask*eps)
!     !-- transform this to F-space
!     call fft ( inx=work(:,:,:,1), outk=uk )  
! 
!   case ("empty")
!     !---------------------------------------------------------------------------
!     !-- initialize scalar zero everywhere
!     work(:,:,:,1) = 0.d0
!     !-- transform this to F-space
!     call fft ( inx=work(:,:,:,1), outk=uk )  
! 
!     
!   case default
!     if (mpirank==0) then
!       write(*,*) "init_passive_scalar:: unknown inicond_scalar="//inicond_scalar
!       call abort()
!     endif
!   end select
  
end subroutine init_passive_scalar