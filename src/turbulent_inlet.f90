module turbulent_inlet_module
  use fsi_vars
  implicit none
  
  ! inflow field
  complex(kind=pr),allocatable,save::uk_turb(:,:,:,:)
  
  !!!!!!!!
  contains
  !!!!!!!!
  
  
!-----------------------------------------------------------------------------
! Initialize the turbulent inflow
!-----------------------------------------------------------------------------
subroutine init_turbulent_inlet ( work )
  use p3dfft_wrapper
  implicit none
  
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr)::u1,u2,u3
  
  if (mpirank==0) write(*,*) "Initializing turbulent inlet..."
  
  !--allocate 
  if (.not.allocated(uk_turb)) then
    allocate(uk_turb(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq))
  endif
  
  !-- read files
  call Read_Single_File ( "ux_turb.h5", work )
  call fft( inx=work, outk=uk_turb(:,:,:,1) )
  
  call Read_Single_File ( "uy_turb.h5", work )
  call fft( inx=work, outk=uk_turb(:,:,:,2) )
  
  call Read_Single_File ( "uz_turb.h5", work )
  call fft( inx=work, outk=uk_turb(:,:,:,3) )
  
  if (mpirank==0) write(*,*) "DONE init turbulent inlet."
  
  
  call get_mean_flow(uk_turb,u1,u2,u3)
  if (mpirank==0) write(*,'("inlet mean flow= ",3(es15.8,1x))') u1,u2,u3
end subroutine
  

!-----------------------------------------------------------------------------
! create the mask function for turbulent inlet
!-----------------------------------------------------------------------------
subroutine turbulent_inlet( time, workc )
  use p3dfft_wrapper
  use penalization
  implicit none
  
  real(kind=pr),intent(in) :: time
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw) 
  complex(kind=pr)::expc
  integer::ix,iy,iz
    
  do iz=ca(1),cb(1)
    do iy=ca(2),cb(2)
        do ix=ca(3),cb(3)
          expc = exp(-dcmplx(0.d0,1.d0)*uxmean*time*wave_x(ix))
          workc(iz,iy,ix,1) = expc*uk_turb(iz,iy,ix,1)
          workc(iz,iy,ix,2) = expc*uk_turb(iz,iy,ix,2)
          workc(iz,iy,ix,3) = expc*uk_turb(iz,iy,ix,3)
      enddo
    enddo
  enddo
  
  
  call ifft( ink=workc(:,:,:,1), outx=us(:,:,:,1) )
  call ifft( ink=workc(:,:,:,2), outx=us(:,:,:,2) )
  call ifft( ink=workc(:,:,:,3), outx=us(:,:,:,3) )
  
  do ix=ra(1),rb(1)
    do iy=ra(2),rb(2)
      do iz=ra(3),rb(3)
        if (ix<=10) then
          mask(ix,iy,iz)=1.d0
          mask_color(ix,iy,iz)=0
        endif
      enddo
    enddo
  enddo  
  
  ! delete flow field outside inlet
  us(:,:,:,1) = us(:,:,:,1)*mask !+ uxmean
  us(:,:,:,2) = us(:,:,:,2)*mask !+ uymean
  us(:,:,:,3) = us(:,:,:,3)*mask !+ uzmean
  
!   us(:,:,:,1) = us(:,:,:,1) + uxmean
!   us(:,:,:,2) = us(:,:,:,2) + uymean
!   us(:,:,:,3) = us(:,:,:,3) + uzmean
!   
end subroutine
  
  
end module turbulent_inlet_module