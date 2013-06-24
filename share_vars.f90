! Variables for pseudospectral simnulations
module vars
  use mpi_header
  implicit none

  character(len=3),save:: method
  integer,save :: nf = 3 ! number of fields (3 for NS, 6 for MHD)

  integer,parameter :: pr = 8 ! precision of doubles
  integer,save :: nx,ny,nz,nt ! resolution and time-stepping

  ! Local array bounds
  integer,dimension (1:3),save :: ra,rb,rs,ca,cb,cs

  ! MPI and p3dfft variables and parameters
  integer,save :: mpisize,mpirank,mpicommcart
  integer,parameter :: mpiinteger=MPI_INTEGER
  integer,parameter :: mpireal=MPI_DOUBLE_PRECISION
  integer,parameter :: mpicomplex=MPI_DOUBLE_COMPLEX
  integer,dimension(2),save :: mpidims,mpicoords,mpicommslab
  integer,dimension (:,:),allocatable,save :: ra_table,rb_table

  ! used in params.f90
  integer,parameter :: nlines=2048 ! maximum number of lines in PARAMS-file

  ! Variables set via the parameters file
  integer,save :: iDealias,itdrag
  integer,save :: iDrag,iKinDiss
  integer,save :: iSaveVelocity,iSaveVorticity,iSavePress,iSaveMask
  integer,save :: iMoving,iPenalization,iMeanFlow,iDoBackup,iSaveSolidVelocity
  real (kind=pr),save :: tmax,cfl,nu,eps,pi,scalex,scaley,scalez,length
  real (kind=pr),save :: xl,yl,zl,dx,dy,dz
  real (kind=pr),save :: Ux,Uy,Uz,tstart,tsave
  real (kind=pr),save :: dt_fixed
  character (len=80),save :: iMask,iTimeMethodFluid,inicond
 
  ! The mask array.  TODO: move out of shave_vars?
  real (kind=pr),dimension (:,:,:),allocatable,save :: mask
  ! Velocity field inside the solid.  TODO: move out of shave_vars?
  real (kind=pr),allocatable,save :: us(:,:,:,:) 

  ! Vabiables timing statistics.  Global to simplify syntax.
  real (kind=pr),save :: time_fft,time_ifft,time_vis,time_mask,time_fft2
  real (kind=pr),save :: time_vor,time_curl,time_p,time_nlk,time_u, time_ifft2
  real (kind=pr),save :: time_bckp,time_save,time_total,time_fluid,time_nlk_fft
end module vars


! Variables for fsi simulations
module fsi_vars
  use vars
  implicit none

  real (kind=pr),save :: x0,y0,z0 ! Parameters for logical centre of obstacle

  ! The derived quantities for fluid-structure interactions.
  type Integrals
     real(kind=pr) :: time
     real(kind=pr) :: E_Kin
     real(kind=pr) :: Dissip
     real(kind=pr) :: Divergence
     real(kind=pr) :: Volume
     real(kind=pr),dimension(1:3) :: Force
  end type Integrals

  type(Integrals),save :: GlobalIntegrals
end module fsi_vars


! Variables for mhd simulations
module mhd_vars
  use vars
  implicit none

  ! The derived quantities for mhd.  
  type Integrals
     real(kind=pr) :: time
     real(kind=pr) :: E_Kin
     real(kind=pr) :: B_Kin
     real(kind=pr) :: Dissip
     real(kind=pr) :: Divergence
     real(kind=pr) :: Volume
     real(kind=pr),dimension(1:3) :: Force
  end type Integrals

  type(Integrals),save :: GlobalIntegrals
end module mhd_vars


! Compute the FFT of the real-valued 3D array inx and save the output
! in the complex-valued 3D array outk.
subroutine fft(outk,inx)
    use mpi_header
    use vars ! For precision specficiation and array sizes
    
    real(kind=pr),intent(in)::inx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    complex(kind=pr),intent(out)::outk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

    call coftxyz(inx,outk)
end subroutine fft


! Compute the inverse FFT of the complex-valued 3D array ink and save the
! output in the real-valued 3D array outx.
subroutine ifft(outx,ink)
    use mpi_header
    use vars ! For precision specficiation and array sizes
    
    real(kind=pr),intent(out)::outx(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
    complex(kind=pr),intent(in)::ink(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3))

    call cofitxyz(ink,outx)
end subroutine ifft
