!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper of flexible wings
!-------------------------------------------------------------------------------
subroutine flexible_wing_motions ( time, wing, kine)
  implicit none

  real(kind=pr),intent(in) :: time
  type(wingkinematics), intent(inout) :: kine
  type(flexible_wing), intent (inout) :: wing
  integer :: i

  select case(wing%Motion)
  case ("simple_harmonic")
    call simple_harmonic_motion (time, wing)
  case ("stationary")
    continue
  case ("from_file")
    call read_kine_from_file (time, wing, kine)
  case ("revolving_wing")
    call revolving_wing(time,wing)
  case ("harmonic_ocsillation")
    wing%vt0(1) = 0.d0
    wing%vt0(2) = 0.25*(1*pi)*cos(1*pi*time)
    wing%vt0(3) = 0.d0

    wing%at0(1) = 0.d0
    wing%at0(2) = - 0.25*(1*pi)**2*sin(1*pi*time)
    wing%at0(3) = 0.d0
  end select

end subroutine

subroutine revolving_wing (time, wing)

    implicit none

    real(kind=pr),intent(in) :: time
    type(flexible_wing), intent (inout) :: wing
    integer :: j
    real(kind=pr) :: phi_y, ttau
    real(kind=pr) :: phi, phi_dt, alpha, alpha_dt, theta, theta_dt
    real(kind=pr),dimension(1:3,1:3) :: mat_Rx
    real(kind=pr),dimension(1:3) :: u

    ! revolving wing kinematics, pre-defined set. We fix alpha to 45deg and increase
    ! phi linearily with a short startup conditioner as suggested in [1]. The startup
    ! time is fixed to 0.4, which gives phi=31.35deg at the end of that interval
    ! [3] D. Kolomenskiy, Y. Elimelech and K. Schneider. Leading-edge vortex shedding from rotating wings. Fluid Dyn. Res., 46, 031421, 2014.
    ttau = 0.4
    ! position angle (is directly given in radiant)
    ! we use PHI_DOT = 1 as normalization as well (since we have no frequency in this case)
    phi = 1.d0*( ttau*dexp(-time/ttau) + time) - ttau
    phi_dt = 1.d0*(1.d0-dexp(-time/ttau))
    ! feathering angle is constant
    alpha = deg2rad(-45.d0)
    alpha_dt = 0.d0
    ! elevation angle is always zero
    theta = 0.d0
    theta_dt = 0.d0

    wing%phi = phi
    wing%alpha = alpha
    wing%theta = theta

    wing%phi_dt = phi_dt
    wing%alpha_dt = alpha_dt
    wing%theta_dt = theta_dt

    !wing%ar0 = 0.d0
    !wing%ar0(3) = 1.0d0/ttau*dexp(-time/ttau)


    !call rotate_vector_into_global_system(wing,wing%ar0)
    !wing%vr0(1) = phi_dt
    !wing%vr0(2) = alpha_dt
    !wing%vr0(3) = theta_dt

    !This angular velocity is defined in the wing system, we need to put it back into
    !the coordinate system of the MSM solver
    !wing%vr0 = matmul(wing%M_solver,wing%vr0)

end subroutine

subroutine read_kine_from_file (time, wing, kine)

    implicit none

    real(kind=pr),intent(in) :: time
    type(wingkinematics), intent(inout) :: kine
    type(flexible_wing), intent (inout) :: wing
    real(kind=pr) :: phi, phi_dt, alpha, alpha_dt, theta, theta_dt
    type(inifile) :: kinefile

    !---------------------------------------------------------------------------
    ! Load kinematics for one stroke from file. This can be applied for example
    ! in hovering. if the wing motion varies appreciably between strokes,
    ! the kinematic loader is the method of choice. The file is specified
    ! in the params file and stored in Insect%infile. An example file
    ! (kinematics_fourier_example.ini) is in the git-repository
    !---------------------------------------------------------------------------
    if (index( kine%infile,".ini")==0) then
      call abort(2030,"you're trying to load an old kinematics file,please convert it to &
      &a new *.ini file (see src/insects/kinematics_example.ini&
      & the insects-tools repository can help you do that!")
    endif

    !---------------------------------------------------------------------------
    ! this block is excecuted only once
    !---------------------------------------------------------------------------
    if (.not.kine%initialized) then
        if (root) then
          write(*,'(80("<"))')
          write(*,*) "Initializing wing kinematics!"
          write(*,*) "*.ini file is: "//trim(adjustl(kine%infile))
          write(*,'(80("<"))')
        endif
      ! parse ini file
      call read_ini_file_mpi(kinefile, kine%infile, .true.)

      ! how to interpret numbers: Fourier or Hermite?
      call read_param_mpi(kinefile,"kinematics","type",kine%infile_type,"none")
      ! what units are given, degree or radiant?
      call read_param_mpi(kinefile,"kinematics","units",kine%infile_units,"degree")
      ! what convention/definition does the data follow?
      call read_param_mpi(kinefile,"kinematics","convention",kine%infile_convention,"flusi")

      ! inform about your interpretation
      select case (kine%infile_type)
      case ("Fourier","fourier","FOURIER")
        if (root) write(*,*) "The input file is interpreted as FOURIER coefficients"
      case ("Hermite","hermite","HERMITE")
        if (root) write(*,*) "The input file is interpreted as HERMITE coefficients"
      case default
        call abort(77771, "kinematics file does not appear to be valid, set type=fourier or type=hermite")
      end select

      ! how many coefficients will be read
      call read_param_mpi(kinefile,"kinematics","nfft_phi",kine%nfft_phi,0)
      call read_param_mpi(kinefile,"kinematics","nfft_alpha",kine%nfft_alpha,0)
      call read_param_mpi(kinefile,"kinematics","nfft_theta",kine%nfft_theta,0)

      ! read coefficients
      call read_param_mpi(kinefile,"kinematics","a0_phi",kine%a0_phi,0.d0)
      call read_param_mpi(kinefile,"kinematics","a0_alpha",kine%a0_alpha,0.d0)
      call read_param_mpi(kinefile,"kinematics","a0_theta",kine%a0_theta,0.d0)

      call read_param_mpi(kinefile,"kinematics","ai_phi",kine%ai_phi(1:kine%nfft_phi))
      call read_param_mpi(kinefile,"kinematics","bi_phi",kine%bi_phi(1:kine%nfft_phi))
      call read_param_mpi(kinefile,"kinematics","ai_alpha",kine%ai_alpha(1:kine%nfft_alpha))

      call read_param_mpi(kinefile,"kinematics","bi_alpha",kine%bi_alpha(1:kine%nfft_alpha))
      call read_param_mpi(kinefile,"kinematics","ai_theta",kine%ai_theta(1:kine%nfft_theta))
      call read_param_mpi(kinefile,"kinematics","bi_theta",kine%bi_theta(1:kine%nfft_theta))
      kine%initialized = .true.
      call clean_ini_file_mpi( kinefile )

      if (root) write(*,'(80(">"))')
    endif

    !---------------------------------------------------------------------------
    ! get actual kinematics from the coefficients
    !---------------------------------------------------------------------------
    ! this block is executed every time. it is called a few times only per time step
    ! so don't worry about performance, we can use the string comparison
    select case (kine%infile_type)
    case ("Fourier","fourier","FOURIER")
      ! evaluate fourier series
      call fseries_eval(time,phi,phi_dt, kine%a0_phi, kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
      call fseries_eval(time,alpha,alpha_dt, kine%a0_alpha, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
      call fseries_eval(time,theta,theta_dt, kine%a0_theta, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    case ("Hermite","hermite","HERMITE")
      ! evaluate hermite interpolation
      call hermite_eval(time,phi,phi_dt    , kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
      call hermite_eval(time,alpha,alpha_dt, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
      call hermite_eval(time,theta,theta_dt, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    case default
      call abort(1717,"kinematics file does not appear to be valid, set type=fourier or type=hermite")
    end select

    !---------------------------------------------------------------------------
    ! make sure the output is in the right units (it HAS to be radiants!)
    !---------------------------------------------------------------------------
    select case (kine%infile_units)
    case ("degree","DEGREE","Degree")
      ! the rest of the code gets radiants, so convert here
      phi    = deg2rad(phi)
      alpha  = deg2rad(alpha)
      theta  = deg2rad(theta)
      phi_dt = deg2rad(phi_dt)
      alpha_dt = deg2rad(alpha_dt)
      theta_dt = deg2rad(theta_dt)
    case ("radiant","RADIANT","Radiant")
      ! if the file is already in radiants, do nothing and be happy!
    case default
      call abort(1718,"kinematics file does not appear to be valid, set units=degree or units=radiant")
    end select

    !---------------------------------------------------------------------------
    ! make sure convention / definition of angles is respected
    !---------------------------------------------------------------------------
    select case (kine%infile_convention)
    case ("FLUSI","flusi","Flusi","FluSI")
      ! defintion as used in flusi. all angles are positive in the right hand rule
      ! that means especially that positive deviation puts wing downwards in a
      ! horizontal stroke plane
    case ("SISC","sisc","siam-sisc")
      ! the same as FLUSI, except for the deviation angle. The sign has been changed
      ! to be more in agreement with other people's definitions. however, flusi
      ! always internally works with the right hand rule
      theta    = -theta
      theta_dt = -theta_dt
    case default
      call abort(1719,"kinematics file does not appear to be valid, set convention=flusi or convention=sisc")
    end select

    ! Transger the kinematics into Wing variable
    wing%phi = phi
    wing%alpha = alpha
    wing%theta = theta

    wing%phi_dt = phi_dt
    wing%alpha_dt = alpha_dt
    wing%theta_dt = theta_dt

end subroutine

subroutine simple_harmonic_motion (time, wing)

  implicit none

  real(kind=pr),intent(in) :: time
  type(flexible_wing), intent (inout) :: wing
  integer :: j

    do j=1,nVeins_BC
      wing%z_BC(-1,j) = wing%z0_BC(-1,j) + 0.25*sin(1*pi*time)
      wing%z_BC(0,j) = wing%z0_BC(0,j) + 0.25*sin(1*pi*time)
    enddo

end subroutine

subroutine translation_acceleration_of_wing_plane (time,dt0,dt1,it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), dimension (1:nWings), intent(inout) :: wings
integer :: i

do i=1,nWings
    wings(i)%at0(1) = 0.d0
    wings(i)%at0(2) = - 0.25*(1*pi)**2*sin(1*pi*time)
    wings(i)%at0(3) = 0.d0
enddo

end subroutine

subroutine moving_noninertial_frame_in_reference_frame(time,dt0,dt1, it,wings)

implicit none

real(kind=pr),intent(in) :: time,dt0,dt1
integer,intent(in) :: it
type(flexible_wing), intent(inout) :: wings
real(kind=pr) :: c1, c2, c3, r
integer :: np

np = wings%np

!Calculate scheme coefficients for variable-step BDF2 scheme
r = dt1/dt0
c1=(1+r)**2/(1+2*r)
c2=r**2/(1+2*r)
c3=(1+r)/(1+2*r)

if (it==0) then

  wings%x(1:np) = wings%x(1:np) + dt1*wings%vt0(1) !dt1**2*wings%at0(1)
  wings%y(1:np) = wings%y(1:np) + dt1*wings%vt0(2)
  wings%z(1:np) = wings%z(1:np) + dt1*wings%vt0(3)

else

  wings%x(1:np) = wings%x(1:np) + dt1*wings%vt0(1) !dt1**2*wings%at0(1)
  wings%y(1:np) = wings%y(1:np) + dt1*wings%vt0(2)
  wings%z(1:np) = wings%z(1:np) + dt1*wings%vt0(3)

endif

end subroutine
