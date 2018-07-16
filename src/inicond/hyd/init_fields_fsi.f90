! Set initial conditions for fsi code.
subroutine init_fields_fsi(time,it,dt0,dt1,n0,n1,uk,nlk,vort,explin,workc,&
           press,scalars,scalars_rhs,Insect,beams, work, u)
  use module_ini_files_parser_mpi
  use vars
  use p3dfft_wrapper
  use solid_model
  use insect_module
  use basic_operators
  implicit none

  integer,intent (inout) :: n1,it,n0
  real (kind=pr),intent (inout) :: time,dt1,dt0
  complex(kind=pr),intent(inout)::uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq)
  ! the workc array is not always allocated, ensure allocation before using
  complex(kind=pr),intent(inout)::workc(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:ncw)
  complex(kind=pr),intent(inout)::nlk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:neq,0:nrhs-1)
  real(kind=pr),intent(inout)::vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)
  real(kind=pr),intent(inout)::explin(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:nf)
  real(kind=pr),intent(inout)::press(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::scalars(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars)
  real(kind=pr),intent(inout)::scalars_rhs(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:n_scalars,0:nrhs-1)
  real(kind=pr),intent(inout)::work(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nrw)
  real(kind=pr),intent(inout)::u(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:nd)

  real(kind=pr),dimension(:,:,:),allocatable::tmp
  type(solid),dimension(1:nBeams), intent(inout) :: beams
  type(diptera),intent(inout)::Insect
  integer :: ix,iy,iz, nxs,nys,nzs, nxb,nyb,nzb,k
  real (kind=pr) :: x,y,z,r,a,b,gamma0,x00,r00,omega,viscosity_dummy
  real (kind=pr) :: uu,Ek,E,Ex,Ey,Ez,kx,ky,kz,theta1,theta2,phi,kabs,kh,kp,maxdiv
  complex(kind=pr) :: alpha,beta
  real(kind=pr), dimension(0:nx-1) :: S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin, kvec
  real(kind=pr), dimension(:,:), allocatable :: spec_array

  ! Assign zero values
  time = 0.0d0
  dt1  = tsave
  it = 0

  uk = dcmplx(0.0d0,0.0d0)
  nlk = dcmplx(0.0d0,0.0d0)
  explin = 0.0d0
  vort = 0.0d0

  select case(inicond)
  case ("couette")
    !--------------------------------------------------
    ! couette flow
    !--------------------------------------------------
    R1=0.4d0
    R2=1.0d0
    omega=1.25d0

    a = omega*(-R1**2 / (R2**2 - R1**2))
    b = omega*(R1**2 * R2**2) / (R2**2 - R1**2)

    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy - 0.5d0*yl
        z = dble(iz)*dz - 0.5d0*zl
        R = dsqrt(y**2 + z**2)

        if ((R>R1).and.(R<R2)) then
          ! fluid domain
          uu = a*R + b/R
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) =+uu*z/R
          vort(:,iy,iz,3) =-uu*y/R
        elseif (R>=R2) then
          ! outer cylinder
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) = 0.d0
          vort(:,iy,iz,3) = 0.d0
        elseif (R<=R1) then
          ! inner cylinder
          vort(:,iy,iz,1) = 0.d0
          vort(:,iy,iz,2) = +omega*z
          vort(:,iy,iz,3) = -omega*y
        endif
      enddo
    enddo
    call fft3 ( uk,vort )

  case("random-given-spectrum")
    if (root) write (*,*) "*** inicond: random field with given spectrum"

    ! step 1
    ! create a random vorticity field in phys space. there is no need to make it
    ! div-free or regular or anything. even the energy does not matter - just some
    ! noise.
    call random_seed()
    do iz=ra(3), rb(3)
      do iy=ra(2), rb(2)
        do ix=ra(1), rb(1)
          vort(ix,iy,iz,1) = rand_nbr()
          vort(ix,iy,iz,2) = rand_nbr()
          vort(ix,iy,iz,3) = rand_nbr()
        end do
      end do
    end do

    ! step 2
    ! bring this field to k-space an project it on the incompressible manifold.
    ! we end up with a velocity field which is random and div-free, but other than
    ! that has no remarkable property-
    call fft3(inx=vort, outk=nlk(:,:,:,:,0) )
    call Vorticity2Velocity(nlk(:,:,:,1:3,0), uk(:,:,:,1:3))

    ! step 3
    ! compute the spectrum of the new field, and define the spectrum that you want
    ! to have.
    call compute_spectrum( time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

    ! read dsired spectrum from file.
    call count_lines_in_ascii_file_mpi(inicond_spectrum_file, k, n_header=0)
    ! read the array from file, two columns, one is k second is E(k)
    allocate(spec_array(0:k-1,1:2))
    call read_array_from_ascii_file_mpi(inicond_spectrum_file, spec_array, n_header=0)

    ! step 4
    ! in wavenumber space, multiply all fourier coeffs with the ratio of desired to
    ! current spectrum (theres a sqrt in it). Now the field will still be div-free
    ! but have the spectrum that you impose.
    do iz=ca(1),cb(1)
      kz = wave_z(iz)
      do iy=ca(2),cb(2)
        ky = wave_y(iy)
        do ix=ca(3),cb(3)
          kx = wave_x(ix)
          ! compute magnitude of wavenumber vector
          kabs = dsqrt(kx**2 + ky**2 + kz**2)
          ! index of spectrum (just the rounded wavenumber, if domain is 2*pi ^3. if not
          ! the routine will fail.)
          k = nint(kabs)
          ! scale spectrum. (avoid dividing by zero...)
          if (S_Ekin(k) > 1.0d-13) then
            uk(iz,iy,ix,1:3) = uk(iz,iy,ix,1:3) * dsqrt(spec_array(k,2) / S_Ekin(k))
          else
            uk(iz,iy,ix,1:3) = 0.0d0
          endif
        enddo
      enddo
    enddo

    ! step 5
    ! done. no some checks:
    deallocate( spec_array )

    ! check what mean energy we end up with
    call compute_spectrum( time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )

    ! check if our initial condition is indeed divergence-free
    call divergence( ink=uk, outk=workc(:,:,:,1) )
    call ifft( ink=workc(:,:,:,1), outx=vort(:,:,:,1) )
    ! vort(:,:,:,1) is now div in phys space
    maxdiv = fieldmax(vort(:,:,:,1))

    if(root) then
      write(*,*) "Mean energy of field=", sum(S_Ekin)
      write(*,*) "Please note that output of integrals is E=",sum(S_Ekin)*(2.d0*pi)**3
      write(*,*) "because it is the volume integral and not the mean"
      write(*,*) "Maximum divergence in field=", maxdiv
      write(*,'(80("-"))')
    endif

  case("turbulence_rogallo")
    !---------------------------------------------------------------------------
    ! randomized initial condition with given spectrum k^4*exp(-k^2 /2)
    ! the field is normalized to have the given energy "omega1"=E, where omega1
    ! is set in parameter file
    ! NOTE: I actually did not figure out what happens if xl=yl=zl/=2*pi
    ! which is a rare case in all isotropic turbulence situtations, and neither
    ! of the corresponding routines have been tested for that case.
    ! NOTE: This technique is validated; just give it a spectrum and see
    ! if that is what you get in return
    !---------------------------------------------------------------------------

    if (root) then
      write(*,'(80("-"))')
      write(*,*) "Initial condition: turbulence_rogallo"
      write(*,*) "randomized divergence-free field with given spectum"
      write(*,*) "and total energy ", omega1
    endif

    ! initialize the random number generator
    call random_seed()

    ! loop over wavenumbers
    do iz=ca(1),cb(1)
      kz = wave_z(iz)
      do iy=ca(2),cb(2)
        ky = wave_y(iy)
        do ix=ca(3),cb(3)
          kx = wave_x(ix)
          ! The spectrum E(k) is defined for "binned" wavenumbers, i.e.
          ! for example 0.5 <= k <= 1.5 is K=1 (this is a spherical wavenumber
          ! shell in k-space). this operation can easily be achieved by rounding
          ! kabs to the nearest integer and then re-converting it to double precision
          kabs = dble( nint(dsqrt(kx**2 + ky**2 + kz**2)) )
          kh   = dsqrt( kx**2 + ky**2)
          ! the wavenumber kp has the highest energy in the field. Mostly, you will
          ! want to set either 1.d0 or 2.d0
          kp   = 2.d0

          ! random numbers for rogallo's initial conditions. Note we use the
          ! wrapper rand_nbr which is in vars.f90
          theta1 = 2.d0*pi*rand_nbr()
          theta2 = 2.d0*pi*rand_nbr()
          phi    = 2.d0*pi*rand_nbr()

          ! this is the given spectrum of the field, which depends on the (radially
          ! rounded) absolute wavenumber: E(k). You can set whatever you like
          ! here, but most people use k^4 * exp(-k^2) (normalized to a given energy)
          Ek = 1.0d2 * (kabs**4 / kp**5) * exp(-2.d0*(kabs/kp)**2)

          ! the coefficients alpha and beta follow exactly Rogallo's notation
          if (kabs /= 0.d0) then
            alpha = 2.d0 * sqrt(Ek/(4.d0*pi*kabs**2)) * exp(cmplx(0.d0,theta1)) * cos(phi)
            beta  = 2.d0 * sqrt(Ek/(4.d0*pi*kabs**2)) * exp(cmplx(0.d0,theta2)) * sin(phi)
          else
            alpha = cmplx(0.d0,0.d0)
            beta  = cmplx(0.d0,0.d0)
          endif

          ! These formulaes come out of Rogallo's NASA 1981 paper, page 53 directly
          ! but note that Rogallo made a typo in the third component (the one
          ! here is correct)
          if (kabs*kh /= 0.d0) then
            uk(iz,iy,ix,1) = ( alpha*kabs*ky + beta*kx*kz) / (kabs*kh)
            uk(iz,iy,ix,2) = (-alpha*kabs*kx + beta*kz*ky) / (kabs*kh)
            uk(iz,iy,ix,3) = (-beta*kh**2) / (kabs*kh)
          else
            uk(iz,iy,ix,1) = cmplx(0.d0,0.d0)
            uk(iz,iy,ix,2) = cmplx(0.d0,0.d0)
            uk(iz,iy,ix,3) = cmplx(0.d0,0.d0)
          endif
       enddo
      enddo
    enddo

    ! enforce hermitian symmetry the lazy way
    call ifft3(ink=uk,outx=vort)
    call fft3(inx=vort,outk=uk)

    ! we now renomalize the velocity, such that it has the given energy omega1
    ! which is set in the parameter file
    call compute_spectrum( time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )
    Ek = sum(S_Ekin)
    if (root) write(*,*) "MEAN Energy before normalization=" ,ek

    ! actual normalization
    uk = uk * dsqrt(omega1/Ek)

    ! check if this really worked
    call compute_spectrum( time,kvec,uk(:,:,:,1:3),S_Ekinx,S_Ekiny,S_Ekinz,S_Ekin )
    if (root) write(*,*) "MEAN Energy after normalization =" ,sum(S_Ekin)
    if (root) write(*,*) "Please note that output of integrals is E=",sum(S_Ekin)*(2.d0*pi)**3
    if (root) write(*,*) "because it is the volume integral and not the mean"


    ! check if our initial condition is indeed divergence-free
    call divergence( ink=uk, outk=workc(:,:,:,1) )
    call ifft( ink=workc(:,:,:,1), outx=vort(:,:,:,1) )
    ! vort(:,:,:,1) is now div in phys space
    maxdiv = fieldmax(vort(:,:,:,1))
    if(mpirank == 0) write(*,*) "Maximum divergence in field=", maxdiv
    if(mpirank == 0) write(*,'(80("-"))')



  case ("taylor_green_2d")
    !--------------------------------------------------
    ! taylor green vortices
    !--------------------------------------------------
    do iz=ra(3),rb(3)
      do iy=ra(2),rb(2)
        y = dble(iy)*dy
        z = dble(iz)*dz
        vort(:,iy,iz,2) = dsin( y ) * dcos( z )
        vort(:,iy,iz,3) =-dcos( y ) * dsin( z )
      enddo
    enddo
    call fft3 ( uk,vort )

  case("infile")
     !--------------------------------------------------
     ! read HDF5 files
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: reading infiles"
     ! fetch time from one file
     call fetch_attributes( file_ux, nxs, nys, nzs, x, y, z, time, viscosity_dummy, origin )
     ! read files
     call Read_Single_File ( file_ux, vort(:,:,:,1) )
     call Read_Single_File ( file_uy, vort(:,:,:,2) )
     call Read_Single_File ( file_uz, vort(:,:,:,3) )
     call fft3 ( uk(:,:,:,1:3), vort(:,:,:,1:3) )

     if (equation=='artificial-compressibility') then
       call Read_Single_File ( file_p, vort(:,:,:,1) )
       call fft( inx=vort(:,:,:,1), outk=uk(:,:,:,4) )
     endif

     if (mpirank==0) write (*,*) "*** done reading infiles"

  case("infile_inlet")
    ! read fields from file, but repeat it in the x-direction, since it is too short
    ! note that the resulting field is non-periodic (but the boundary conditions
    ! will regularize this, i.e. with sponges)

    nxb=nx
    nyb=ny
    nzb=nz

    if (mpirank==0) write(*,*) "Inicond infile_inlet "//file_ux
    if (mpirank==0) write(*,*) "Inicond infile_inlet "//file_uy
    if (mpirank==0) write(*,*) "Inicond infile_inlet "//file_uz

    call fetch_attributes( file_ux, nxs, nys, nzs, x, y, z, time, viscosity_dummy, origin )
    time = 0.d0
    ra(1) = 0
    rb(1) = nxs-1
    nx=nxs
    if(mpirank==0) write(*,*) "Reduced to ", nxs, ra, rb
    allocate (tmp(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3)) )

    call read_single_file( file_ux, tmp )
    vort(0:nxs-1,:,:,1) = tmp
    vort(nxs:nxb-1,:,:,1) = tmp(0:nxb-1-nxs+1,:,:)


    call read_single_file( file_uy, tmp )
    vort(0:nxs-1,:,:,2) = tmp
    vort(nxs:nxb-1,:,:,2) = tmp(0:nxb-1-nxs+1,:,:)

    call read_single_file( file_uz, tmp )
    vort(0:nxs-1,:,:,3) = tmp
    vort(nxs:nxb-1,:,:,3) = tmp(0:nxb-1-nxs+1,:,:)

    deallocate (tmp)

    nx=nxb
    ra(1)=0
    rb(1)=nx-1
    if(mpirank==0) write(*,*) "Back to ", nx, ra, rb

    call fft3(inx=vort,outk=uk)

  case("infile_perturbed")
    ! read velocity field from file, as is done in inicond="infile", but add a
    ! small perturbation to the vorticity, as done in "half_HIT"
    ! useful to read in a turbulence field, add a perturbation, to trigger a different decay
    if (mpirank==0) write (*,*) "*** inicond: reading infiles, and perturbing them"
    call Read_Single_File ( file_ux, vort(:,:,:,1) )
    call Read_Single_File ( file_uy, vort(:,:,:,2) )
    call Read_Single_File ( file_uz, vort(:,:,:,3) )
    if (mpirank==0) write (*,*) "*** done reading infiles"

    call fft3(inx=vort,outk=uk)
    call curl3_inplace(uk)
    call ifft3(ink=uk,outx=vort)
    call compute_energies(vort,E,Ex,Ey,Ez)
    omega1 = 0.02d0 *  sqrt( 2.d0*E / (xl*yl*zl) )

    if(mpirank==0) write(*,*) "Perturbation is ", omega1

    do iz=ra(3), rb(3)
      do iy=ra(2), rb(2)
          do ix=ra(1), rb(1)
            vort(ix,iy,iz,1) = vort(ix,iy,iz,1)+omega1*(2.0d0*rand_nbr() - 1.d0)
            vort(ix,iy,iz,2) = vort(ix,iy,iz,2)+omega1*(2.0d0*rand_nbr() - 1.d0)
            vort(ix,iy,iz,3) = vort(ix,iy,iz,3)+omega1*(2.0d0*rand_nbr() - 1.d0)
          end do
      end do
    end do

    call Vorticity2Velocity_old(uk(:,:,:,1:3), nlk(:,:,:,1:3,0), vort(:,:,:,1:3))

  case("VortexRing")
     !--------------------------------------------------
     ! Vortex ring intial condition
     ! The vortex ring is a poloidal vorticty field, as for example defined in
     ! Jause-Labert et al., "Numerical validation of the volume penalization method in three-dimensional
     ! pseudo-spectral simulations" Computers & Fluids 2012
     ! This ring has the strength omega1, its center is at (x0,y0,z0), and its radius is "length"
     !--------------------------------------------------
     ! radius of vortexring (set in params file)
     r00 = length
     a  = 0.4131d0 * r00
     a  = 0.82d0 * r00
     ! strength (set in params file)
     gamma0 = omega1

     if (mpirank==0) then
       write(*,'(80("*"))')
       write (*,*) "*** inicond: vortex ring initial condition"
       write (*,'("r00=",g12.4," a=",g12.4," gamma0=",g12.4," RE=",g13.4)') r00,a,gamma0,gamma0/nu
       write(*,'("center=",3(g12.4,1x))') x0,y0,z0
       write(*,'(80("*"))')
     endif

     ! define vorticity in phy space
     do iz=ra(3), rb(3)
       do iy=ra(2), rb(2)
         do ix=ra(1), rb(1)
           x = dble(ix)*dx
           y = dble(iy)*dy
           z = dble(iz)*dz
           ! radius in cylinder coordinates. note we do not need the polar
           ! angle, luckily, as we can directly set their cosine/sine with
           ! z/r and y/r. This is much cheaper.
           r = dsqrt( (y-y0)**2 + (z-z0)**2 )

           ! angular vorticity component. The vorticity vector in cylinder
           ! coordinates is (r,theta,z) = (0,omega,0), see eg
           ! Jause-Labert et al., "Numerical validation of the volume penalization method in three-dimensional
           ! pseudo-spectral simulations" Computers & Fluids 2012
           omega = (gamma0/(pi*a**2))*dexp(-((x-x0)**2 + (r-r00)**2)/(a**2) )

           ! avoid division by zero if R is very small
           if ( dabs(r)> 1.0d-12) then
             vort(ix,iy,iz,2) =-omega * ( (z-z0)/r) ! this is sin(theta)
             vort(ix,iy,iz,3) = omega * ( (y-y0)/r) ! this is cos(theta)
           endif
         enddo
       enddo
     enddo
     ! go to Fourier space
     call fft3(inx=vort, outk=nlk(:,:,:,:,0))
     ! invert curl (Biot-Savart Operator)
     call Vorticity2Velocity(nlk(:,:,:,1:3,0),uk(:,:,:,1:3))
     ! add mean flow (note: any inverted curl always has zero mean flow)
     call set_mean_flow(uk,time)

 case ("three-vortices")
     if (nx /= 1) call abort(13, "three-vortices is a 2d flow! set nx=1")

     ! define vorticity in phy space
     do iz=ra(3), rb(3)
         do iy=ra(2), rb(2)
             y = dble(iy)*dy
             z = dble(iz)*dz

             gamma0 = +1.0_pr
             a = 1.0_pr/pi
             y0 = 0.75_pr*pi ! in farge 1997, this is X
             z0 = 1.00_pr*pi ! in farge 1997, this is Y
             omega = (gamma0/(pi*a**2))*dexp(-((y-y0)**2 + (z-z0)**2)/(a**2) )

             gamma0 = +1.0_pr
             a = 1.0_pr/pi
             y0 = 1.25_pr*pi ! in farge 1997, this is X
             z0 = 1.00_pr*pi ! in farge 1997, this is Y
             omega = omega + (gamma0/(pi*a**2))*dexp(-((y-y0)**2 + (z-z0)**2)/(a**2) )

             gamma0 = -0.5_pr
             a = 1.0_pr/pi
             y0 = 1.25_pr*pi ! in farge 1997, this is X
             z0 = 1.00_pr*pi + pi/(2.0_pr*sqrt(2.0_pr)) ! in farge 1997, this is Y
             omega = omega + (gamma0/(pi*a**2))*dexp(-((y-y0)**2 + (z-z0)**2)/(a**2) )

             vort(:,iy,iz,1) = omega
             vort(:,iy,iz,2:3) = 0.0_pr
         enddo
     enddo


     gamma0 = fieldmean(vort(:,:,:,1))
     write(*,*) "mean vort", gamma0
     vort(:,:,:,1) = vort(:,:,:,1) - gamma0

     call Vorticity2Velocity_old(uk(:,:,:,1:3), nlk(:,:,:,1:3,0), vort(:,:,:,1:3))

  case ("vortex")
     if (mpirank==0) write (*,*) "*** inicond: vortex ring initial condition"
     r00=yl/8.d0
     a  =0.4131d0 * r00
     a  =0.82d0 * r00
     gamma0=12.0d0
     x00=0.5d0 * xl

     ! define vorticity in phy space
     do iz=ra(3), rb(3)
        z=dble(iz) * zl / dble(nz)
        do iy=ra(2), rb(2)
           y=dble(iy) * yl / dble(ny)
           do ix=ra(1), rb(1)
              x=dble(ix) * xl / dble(nx)

              r=dsqrt( (x-0.5d0*xl)**2+ (y-0.5d0*yl)**2 + (z-0.5d0*zl)**2 )

              omega = 10.d0*exp( -(r / (0.1*zl))**2 )

              vort (ix,iy,iz,3)=omega
           enddo
        enddo
     enddo

     call Vorticity2Velocity_old(uk(:,:,:,1:3), nlk(:,:,:,1:3,0), vort(:,:,:,1:3))

     call set_mean_flow(uk,time)

  case("turbulence")
     !--------------------------------------------------
     ! random vorticity
     ! This case has two parameters read from ini file:
     ! the maximum vorticity in each direction and the
     ! smoothing parameter
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: turbulence (random vorticity) initial condition"
     call random_seed()
     if (nx>1) then
       do iz=ra(3), rb(3)
          do iy=ra(2), rb(2)
             do ix=ra(1), rb(1)
                vort(ix,iy,iz,1)=omega1*(2.0d0*rand_nbr() - 1.d0)
                vort(ix,iy,iz,2)=omega1*(2.0d0*rand_nbr() - 1.d0)
                vort(ix,iy,iz,3)=omega1*(2.0d0*rand_nbr() - 1.d0)
             end do
          end do
       end do
     else
       do iz=ra(3), rb(3)
          do iy=ra(2), rb(2)
              vort(:,iy,iz,1)=omega1*(2.0d0*rand_nbr() - 1.d0)
              vort(:,iy,iz,2)=0.0d0
              vort(:,iy,iz,3)=0.0d0
          end do
       end do
     end if

     call cal_vis( nu_smoothing/nu, explin(:,:,:,1))
     call fft3( inx=vort, outk=nlk(:,:,:,:,0) )
     nlk(:,:,:,1,0)=nlk(:,:,:,1,0)*explin(:,:,:,1)
     nlk(:,:,:,2,0)=nlk(:,:,:,2,0)*explin(:,:,:,1)
     nlk(:,:,:,3,0)=nlk(:,:,:,3,0)*explin(:,:,:,1)
     call ifft3( ink=nlk(:,:,:,:,0), outx=vort )

     call Vorticity2Velocity_old(uk(:,:,:,1:3), nlk(:,:,:,1:3,0), vort(:,:,:,1:3))

  case("half_HIT")
    ! this is a very specialized case. it reads a field from files, but the field
    ! is only half as long in the x-direction. it is then padded by itself (we
    ! take the same field twice), we compute the curl and add a 0.5% white noise
    ! to the vorticity, recompute the velocity and go from there
    allocate (tmp(0:nx/2-1,ra(2):rb(2),ra(3):rb(3)) )

    if (mpirank==0) write(*,*) "Inicond HALF_HIT "//file_ux
    if (mpirank==0) write(*,*) "Inicond HALF_HIT "//file_uy
    if (mpirank==0) write(*,*) "Inicond HALF_HIT "//file_uz

    nx=nx/2
    rb(1)=nx-1
    if(mpirank==0) write(*,*) "Reduced to ", nx, ra, rb



    call read_single_file( file_ux, tmp )
    vort(0:nx-1,:,:,1) = tmp
    vort(nx:2*nx-1,:,:,1) = tmp

    call read_single_file( file_uy, tmp )
    vort(0:nx-1,:,:,2) = tmp
    vort(nx:2*nx-1,:,:,2) = tmp

    call read_single_file( file_uz, tmp )
    vort(0:nx-1,:,:,3) = tmp
    vort(nx:2*nx-1,:,:,3) = tmp

    deallocate (tmp)

    nx = nx*2
    rb(1)=nx-1
    if(mpirank==0) write(*,*) "Back to ", nx, ra, rb

    call fft3(inx=vort,outk=uk)
    call curl3_inplace(uk)
    call ifft3(ink=uk,outx=vort)

    call compute_energies(vort,E,Ex,Ey,Ez)

    omega1 = 0.02d0 *  sqrt( 2.d0*E / (xl*yl*zl) )

    if(mpirank==0) write(*,*) "Perturbation is ", omega1

    do iz=ra(3), rb(3)
      do iy=ra(2), rb(2)
          do ix=ra(1), rb(1)
            vort(ix,iy,iz,1) = vort(ix,iy,iz,1)+omega1*(2.0d0*rand_nbr() - 1.d0)
            vort(ix,iy,iz,2) = vort(ix,iy,iz,2)+omega1*(2.0d0*rand_nbr() - 1.d0)
            vort(ix,iy,iz,3) = vort(ix,iy,iz,3)+omega1*(2.0d0*rand_nbr() - 1.d0)
          end do
      end do
    end do

    call Vorticity2Velocity_old(uk(:,:,:,1:3), nlk(:,:,:,1:3,0), vort(:,:,:,1:3))

  case("MeanFlow")
     !--------------------------------------------------
     ! mean flow only
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: mean flow"
     uk=dcmplx(0.0d0,0.0d0)
     ! note this inicond also works without meanflow forcing, it is then
     ! really just an inicond

     ! forcing = zeroth Fourier mode only
     if ( (ca(1) == 0) .and. (ca(2) == 0) .and. (ca(3) == 0) ) then
       uk(0, 0, 0,1) = Uxmean
       uk(0, 0, 0,2) = Uymean
       uk(0, 0, 0,3) = Uzmean
     endif

  case("quiescent")
     !--------------------------------------------------
     ! fluid at rest
     !--------------------------------------------------
     if (mpirank==0) write (*,*) "*** inicond: fluid at rest"
     uk=dcmplx(0.0d0,0.0d0)

  case default
     if(inicond(1:8) == "backup::") then
        !--------------------------------------------------
        ! read from backup
        !--------------------------------------------------
        if (mpirank==0) write (*,*) "*** inicond: retaking backup " // &
             inicond(9:len(inicond))
        call Read_Runtime_Backup(inicond(9:len(inicond)),time,dt0,dt1,n1,it,uk,&
             nlk,explin,vort(:,:,:,1),scalars,scalars_rhs)
     else
        !--------------------------------------------------
        ! unknown inicond : error
        !--------------------------------------------------
        write (*,*) '??? ERROR: Invalid initial condition' // inicond
        call abort(55523, '??? ERROR: Invalid initial condition')
     endif
  end select


  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the solid solver
  !-----------------------------------------------------------------------------
  if (use_solid_model=="yes") then
    if(mpirank==0) write(*,*) "Initializing solid solver and testing..."
    call init_beams( beams )
    call surface_interpolation_testing( time, beams(1), press )
    call init_beams( beams )
  endif

  !-----------------------------------------------------------------------------
  ! If module is in use, initialize also the passive scalar(s)
  !-----------------------------------------------------------------------------
  if ((use_passive_scalar==1).and.(index(inicond,"backup::")==0)) then
    ! only if not resuming a backup
    call init_passive_scalar(scalars,scalars_rhs,Insect,beams)
  endif

  !-----------------------------------------------------------------------------
  ! when computing running time avg, initialize (note that if we're resuming
  ! a backup, it is read from that file)
  !-----------------------------------------------------------------------------
  if ((time_avg=="yes").and.(index(inicond,"backup::")==0)) then
    if (vel_avg=="yes") uk_avg(:,:,:,1:3) = uk(:,:,:,1:3)
    if (ekin_avg=="yes") e_avg=0.d0
  endif
end subroutine init_fields_fsi


! Computes the divergence-free velocity in Fourier space u given vort
! in physical space.  work is a work array
subroutine Vorticity2Velocity_old(uk,work,vort)
  use mpi
  use vars
  use p3dfft_wrapper
  implicit none

  complex (kind=pr),intent(out) :: uk(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  complex (kind=pr),intent(inout)::work(ca(1):cb(1),ca(2):cb(2),ca(3):cb(3),1:3)
  real (kind=pr), intent (in) :: vort(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:3)
  integer :: ix, iy, iz, i
  real (kind=pr) :: kx,ky,kz,kx2,ky2,kz2,k_abs_2
  complex (kind=pr) :: im
  ! imaginary unit
  im=dcmplx(0.d0,1.d0)
  !-------------------------------------------------
  ! Compute vorticity in Fourier space
  !-------------------------------------------------
  call fft3(inx=vort,outk=work)

  !-------------------------------------------------
  ! Compute streamfunction in Fourier space
  ! work(:,:,:,1:3, 1) will contain the three components of
  ! streamfunction
  !-------------------------------------------------
   do ix=ca(3),cb(3)
    kx=wave_x(ix)
    kx2=kx*kx
    do iy=ca(2),cb(2)
      ky=wave_y(iy)
      ky2=ky*ky
      do iz=ca(1),cb(1)
        kz=wave_z(iz)
        kz2=kz*kz

        k_abs_2=kx2+ky2+kz2
        if (abs(k_abs_2) .ne. 0.0) then
          work(iz,iy,ix,1)=+work(iz,iy,ix,1) / k_abs_2
          work(iz,iy,ix,2)=+work(iz,iy,ix,2) / k_abs_2
          work(iz,iy,ix,3)=+work(iz,iy,ix,3) / k_abs_2
        else
          work(iz,iy,ix,1)=dcmplx(0.d0,0.d0)
          work(iz,iy,ix,2)=dcmplx(0.d0,0.d0)
          work(iz,iy,ix,3)=dcmplx(0.d0,0.d0)
        endif
      enddo
    enddo
  enddo

  !-----------------------------------------------
  !-- compute velocity as curl of streamfunction
  !-----------------------------------------------
  do ix=ca(3),cb(3)
    !-- wavenumber in x-direction
    kx = wave_x(ix)
    do iy=ca(2),cb(2)
      !-- wavenumber in y-direction
      ky = wave_y(iy)
      do iz=ca(1),cb(1)
        !-- wavenumber in z-direction
        kz = wave_z(iz)

        uk(iz,iy,ix,1)=im*(ky*work(iz,iy,ix,3) - kz*work(iz,iy,ix,2))
        uk(iz,iy,ix,2)=im*(kz*work(iz,iy,ix,1) - kx*work(iz,iy,ix,3))
        uk(iz,iy,ix,3)=im*(kx*work(iz,iy,ix,2) - ky*work(iz,iy,ix,1))
      enddo
    enddo
  enddo
end subroutine Vorticity2Velocity_old
