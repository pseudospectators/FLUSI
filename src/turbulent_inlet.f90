module turbulent_inlet_module
  use vars
  implicit none

  ! inflow field
  real(kind=pr),allocatable,save :: u_turb(:,:,:,:)

  integer :: nx_turb, ny_turb, nz_turb
  real(kind=pr) :: xl_turb, yl_turb, zl_turb

  !!!!!!!!
  contains
  !!!!!!!!


!-----------------------------------------------------------------------------
! Initialize the turbulent inflow
!-----------------------------------------------------------------------------
subroutine init_turbulent_inlet ( )
  use p3dfft_wrapper
  use basic_operators
  implicit none

  real(kind=pr) :: u1, u2, u3
  integer :: nx_org, L, N, iy, iz, q
  real(Kind=pr),allocatable :: f(:,:), fk(:,:)
  real(kind=pr) :: time, umax,viscosity_dummy

  complex (kind=pr), dimension (:), allocatable :: cworkx, expx0, expx ! work array and shift in x

  if (mpirank==0) write(*,*) "Initializing turbulent inlet..."

  if (iMoving==0) then
    call abort(100, "Using turbulent inlet requires setting iMoving=1 since the us field is not updated otherwise. ")
  endif

  ! The resolution in the inlet file and the resolution of our simulation do not
  ! nessesarily match in the x-direction (they DO match in y-z direction). We
  ! first determine how many points in x-direction the file has, then allocate
  ! the field u_turb accordingly
  call fetch_attributes( "ux_turb.h5", nx_turb,ny_turb,nz_turb,&
       xl_turb,yl_turb,zl_turb,time,viscosity_dummy, origin)


  if (mpirank==0) write(*,*) "inlet file resolution is ",nx_turb, ny_turb, nz_turb
  !HACK
  yl_turb = yl
  zl_turb = zl
  xl_turb = yl*dble(nx_turb)/dble(ny_turb)
  !HACK
  if (mpirank==0) write(*,*) "inlet file domain is ",xl_turb, yl_turb, zl_turb

  ! initialize global constants. note we assume that the file resolution matches
  ! in y,z directions the one of our simulation
  nx_org = nx
  ra(1) = 0
  rb(1) = nx_turb-1
  nx = nx_turb
  if ((ny_turb.ne.ny).or.(nz_turb.ne.nz)) then
    call abort(101, "the turbulent inlet file resolution must match ny,nz and it doesn't: stop.")
  endif

  ! allocate inlet velocity in physical space with the dimensions gathered from
  ! the file
  allocate( u_turb(0:nx_turb-1,ra(2):rb(2),ra(3):rb(3),1:3) )

  ! read files into allocated array
  call Read_Single_File ( "ux_turb.h5", u_turb(:,:,:,1) )
  call Read_Single_File ( "uy_turb.h5", u_turb(:,:,:,2) )
  call Read_Single_File ( "uz_turb.h5", u_turb(:,:,:,3) )

  ! determine the maximum velocity of the perturbation field
  umax = field_max_magnitude(u_turb)
  if (mpirank==0) write (*,*) " turbulence intensity in infile", umax

  u_turb = u_turb * rescale

  umax = field_max_magnitude(u_turb)
  if (mpirank==0) write (*,*) " turbulence intensity, in use", umax

  ! reset dimensions to old value (for the rest of the program)
  nx = nx_org
  rb(1) = nx-1

  if (mpirank==0) write(*,*) "DONE init turbulent inlet."
end subroutine


!-----------------------------------------------------------------------------
! create the mask function for turbulent inlet
!-----------------------------------------------------------------------------
subroutine turbulent_inlet( time )
  use p3dfft_wrapper
  use penalization
  implicit none

  real(kind=pr),intent(in) :: time
  real(kind=pr) :: x_turb, x, dx_turb, C1, C2
  integer :: ix,iy,iz
  integer :: n, i1,i2

  ! thickness of inflow velocity sponge
  n = inlet_thickness
  dx_turb = xl_turb / dble(nx_turb)

  ! note x-direction is contiguous among MPI procs, i.e. it is NOT split
  do ix = 0,n-1
    x = dble(ix)*dx
    x_turb = -uxmean*time + x

    if (x_turb >= xl_turb) x_turb = x_turb - xl_turb*dble(ceiling(abs(x_turb/xl_turb)))
    if (x_turb < 0.d0 )    x_turb = x_turb + xl_turb*dble(ceiling(abs(x_turb/xl_turb)))

    i1 = floor(x_turb/dx_turb)
    i2 = i1 + 1

    ! weights for linear interpolation
    C1 = (dble(i2)*dx_turb - x_turb)  /  dx_turb
    C2 = (x_turb - dble(i1)*dx_turb)  /  dx_turb

    ! ensure periodic indices
    i1 = per(i1,nx_turb)
    i2 = per(i2,nx_turb)

    if ((i1<0 .or. i1>nx_turb-1).or.(i2<0 .or. i2>nx_turb-1)) then
      if(mpirank==0) write(*,*) "indices for inlet fail: ", i1,i2
      call abort(999991, "indices for inlet fail....")
    endif

    mask(ix,:,:) = 1.d0
    mask_color(ix,:,:) = 0
    us(ix,:,:,1) = u_turb(i1,:,:,1)*C1 + u_turb(i2,:,:,1)*C2 + uxmean
    us(ix,:,:,2) = u_turb(i1,:,:,2)*C1 + u_turb(i2,:,:,2)*C2!+ uymean
    us(ix,:,:,3) = u_turb(i1,:,:,3)*C1 + u_turb(i2,:,:,3)*C2!+ uzmean
  enddo

end subroutine


end module turbulent_inlet_module
