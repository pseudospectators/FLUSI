!-------------------------------------------------------------------------------
! WRAPPER Motion protocoll wrapper left wing
!-------------------------------------------------------------------------------
subroutine FlappingMotion_left ( time, Insect )
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  type(diptera) :: Insect
  
  call FlappingMotion ( time, Insect, Insect%FlappingMotion_left, &
       Insect%phi_l, Insect%alpha_l, Insect%theta_l, Insect%phi_dt_l,&
       Insect%alpha_dt_l, Insect%theta_dt_l )
end subroutine FlappingMotion_left


!-------------------------------------------------------------------------------
! WARPPER Motion protocoll wrapper right wing
!-------------------------------------------------------------------------------
subroutine FlappingMotion_right ( time, Insect )
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  type(diptera) :: Insect
  
  call FlappingMotion ( time, Insect, Insect%FlappingMotion_right, &
       Insect%phi_r, Insect%alpha_r, Insect%theta_r, Insect%phi_dt_r, &
       Insect%alpha_dt_r, Insect%theta_dt_r )
end subroutine FlappingMotion_right



!-------------------------------------------------------------------------------
! Flapping wing motion protocoll, different choices.
! Input: 
!       time (self explanatory)
!       protocoll: string containing what motion you want to have (may be 
!                  different for both wings)
! Output:
!       phi: flapping/positonal angle
!       alpha: feathering angle / angle of attack
!       theta: deviation angle / Out-of-stroke-plane
!       phi_dt: flapping angular velocity
!       alpha_dt: feathering angular velocity
!       theta_dt: deviatory angular velocity
! The actual motion depends on the choices in the parameter file, namely
! Insect%WingMotion, and sub-parameters that may further precise a given motion 
! protocoll. Note we allow both wings to follow a differen motion, but they both 
! call this routine here.
!-------------------------------------------------------------------------------
subroutine FlappingMotion(time, Insect, protocoll, phi, alpha, theta, phi_dt, alpha_dt, theta_dt)
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=pr), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  character (len=*), intent(in) :: protocoll
  real(kind=pr) :: phi_max,alpha_max, phase,f
  real(kind=pr) :: bi_alpha_flapper(1:29) ! For comparison with Sane&Dickinson
  real(kind=pr) :: ai_phi_flapper(1:31) ! For comparison with Sane&Dickinson
  real(kind=pr) :: tadv,t0adv ! For comparison with Dickinson, Ramamurti
  real(kind=pr) :: posi,elev,feth,posi_dt,elev_dt,feth_dt,angles ! Comp. w. Maeda
  real(kind=pr) :: dangle_posi,dangle_elev,dangle_feth ! Comp. w. Maeda
  real(kind=pr) :: dangle_posi_dt,dangle_elev_dt,dangle_feth_dt ! Comp. w. Maeda
  real(kind=pr) :: a_posi(1:4),b_posi(1:4),a_elev(1:4),b_elev(1:4),a_feth(1:4),b_feth(1:4)
  real(kind=pr) :: s,c,a0_phi
  real(kind=pr) :: phicdeg
  real(kind=pr) :: alphacdeg
  integer :: i,mpicode
  character(len=strlen) :: dummy
  
  select case ( protocoll )
  case ("from_file")
    !---------------------------------------------------------------------------
    ! Load kinematics for one stroke from file. This can be applied for example
    ! in hovering. if the wing motion varies appreciably between strokes,
    ! the kinematic loader is the method of choice. The file is specified
    ! in the params file and stored in Insect%infile. Input Fourier coeffients
    ! are in DEGREE, output of this function is IN RADIANT!  An example file 
    ! (kinematics_fourier_example.in) is in the git-repository
    !---------------------------------------------------------------------------
    if (.not.allocated(ai_phi)) then
      ! this block is excecuted only once
      call check_file_exists( Insect%infile )
      ! learn how many Fourier coefficients to expect
      if (mpirank==0) then
        open(37, file=Insect%infile, form='formatted', status='old')
        read(37,*) dummy ! skip lines added for documentation
        read(37,*) Insect%nfft_phi
        read(37,*) Insect%nfft_alpha
        read(37,*) Insect%nfft_theta
        write(*,'("Reading kinematics: nfft_phi=",i2," nfft_alpha=",i2," nfft_theta=",i2)')&
          Insect%nfft_phi, Insect%nfft_alpha, Insect%nfft_theta
      endif
      ! BCAST nfft to all procs
      call MPI_BCAST( Insect%nfft_phi,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%nfft_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%nfft_theta,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      
      ! allocate fourier coefficent arrays
      allocate ( ai_phi(1:Insect%nfft_phi) )
      allocate ( bi_phi(1:Insect%nfft_phi) ) 
      allocate ( ai_alpha(1:Insect%nfft_alpha) )
      allocate ( bi_alpha(1:Insect%nfft_alpha) )
      allocate ( ai_theta(1:Insect%nfft_theta) )
      allocate ( bi_theta(1:Insect%nfft_theta) )
      
      ! read coefficients
      if (mpirank==0) then
        read(37,*) dummy ! skip lines added for documentation
        read(37,*) dummy
        read(37,*) Insect%a0_phi
        read(37,*) dummy
        read(37,*) ai_phi
        read(37,*) dummy
        read(37,*) bi_phi

        read(37,*) dummy
        read(37,*) dummy
        read(37,*) Insect%a0_alpha
        read(37,*) dummy
        read(37,*) ai_alpha
        read(37,*) dummy
        read(37,*) bi_alpha

        read(37,*) dummy
        read(37,*) dummy
        read(37,*) Insect%a0_theta
        read(37,*) dummy
        read(37,*) ai_theta
        read(37,*) dummy
        read(37,*) bi_theta

        close(37)
      endif
      
      call MPI_BCAST( Insect%a0_phi,1,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( ai_phi,Insect%nfft_phi,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_phi,Insect%nfft_phi,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%a0_alpha,1,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( ai_alpha,Insect%nfft_alpha,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_alpha,Insect%nfft_alpha,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%a0_theta,1,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( ai_theta,Insect%nfft_theta,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_theta,Insect%nfft_theta,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
    endif
      
    ! this block is executed every time
    call fseries_eval(time,phi,phi_dt,Insect%a0_phi,ai_phi,bi_phi)
    call fseries_eval(time,alpha,alpha_dt,Insect%a0_alpha,ai_alpha,bi_alpha)
    call fseries_eval(time,theta,theta_dt,Insect%a0_theta,ai_theta,bi_theta)
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)
    
  case ("from_file_hermite")
    !---------------------------------------------------------------------------
    ! motion protocol from hermite coeffients (instead of Fourier coefficients)
    ! it is designed closely to "from_file", Insect%ai... is function values and
    ! Insect%bi is derivatives. Insect%a0 is unused.
    ! Input values are assumed IN DEGREE (°), output of this function is 
    ! IN RADIANTS! An example file (kinematics_hermite_example.in) is in the
    ! git-repository
    !---------------------------------------------------------------------------
    if (.not.allocated(ai_phi)) then
      ! this block is excecuted only once
      call check_file_exists( Insect%infile )
      ! learn how many Fourier coefficients to expect
      if (mpirank==0) then
        open(37, file=Insect%infile, form='formatted', status='old')
        read(37,*) dummy ! skip lines added for documentation
        read(37,*) Insect%nfft_phi
        read(37,*) Insect%nfft_alpha
        read(37,*) Insect%nfft_theta
        write(*,'("Reading kinematics: n_phi=",i2," n_alpha=",i2," n_theta=",i2)')&
          Insect%nfft_phi, Insect%nfft_alpha, Insect%nfft_theta
        write(*,*) "Using Hermite interpolation...."
      endif
      ! BCAST nfft to all procs
      call MPI_BCAST( Insect%nfft_phi,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%nfft_alpha,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( Insect%nfft_theta,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
      
      ! allocate fourier coefficent arrays
      allocate ( ai_phi(1:Insect%nfft_phi) )
      allocate ( bi_phi(1:Insect%nfft_phi) ) 
      allocate ( ai_alpha(1:Insect%nfft_alpha) )
      allocate ( bi_alpha(1:Insect%nfft_alpha) )
      allocate ( ai_theta(1:Insect%nfft_theta) )
      allocate ( bi_theta(1:Insect%nfft_theta) )
      
      ! read coefficients
      if (mpirank==0) then
        read(37,*) dummy ! skip lines added for documentation
        read(37,*) dummy
        read(37,*) ai_phi
        read(37,*) dummy
        read(37,*) bi_phi
        read(37,*) dummy
        read(37,*) dummy
        read(37,*) ai_alpha
        read(37,*) dummy
        read(37,*) bi_alpha
        read(37,*) dummy
        read(37,*) dummy
        read(37,*) ai_theta
        read(37,*) dummy
        read(37,*) bi_theta

        close(37)
      endif
      
      call MPI_BCAST( ai_phi,Insect%nfft_phi,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_phi,Insect%nfft_phi,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( ai_alpha,Insect%nfft_alpha,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_alpha,Insect%nfft_alpha,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( ai_theta,Insect%nfft_theta,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
      call MPI_BCAST( bi_theta,Insect%nfft_theta,MPI_DOUBLE_PRECISION,0,&
           MPI_COMM_WORLD,mpicode )
    endif
      
    ! this block is executed every time
    call hermite_eval(time,phi,phi_dt,ai_phi,bi_phi)
    call hermite_eval(time,alpha,alpha_dt,ai_alpha,bi_alpha)
    call hermite_eval(time,theta,theta_dt,ai_theta,bi_theta)
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)
  
  case ("Drosophila_hovering_fry")
    !---------------------------------------------------------------------------
    ! motion protocoll digitalized from Fry et al JEB 208, 2303-2318 (2005)
    !
    ! fourier coefficients analyzed with matlab
    !---------------------------------------------------------------------------
    if (.not.allocated(ai_alpha)) then
      Insect%nfft_alpha = 10
      Insect%nfft_theta = 10
      Insect%nfft_phi   = 10
      allocate ( ai_phi(1:Insect%nfft_phi), bi_phi(1:Insect%nfft_phi) ) 
      allocate ( ai_alpha(1:Insect%nfft_alpha), bi_alpha(1:Insect%nfft_alpha) )
      allocate ( ai_theta(1:Insect%nfft_theta), bi_theta(1:Insect%nfft_theta) )
      
      Insect%a0_phi   =25.4649398
      Insect%a0_alpha =-0.3056968
      Insect%a0_theta =-17.8244658  ! - sign (Dmitry, 10 Nov 2013)
      
      ai_phi   =(/71.1061858,2.1685448,-0.1986978,0.6095268,-0.0311298,&
                -0.1255648,-0.0867778,0.0543518,0.0,0.0/)
      bi_phi   =(/5.4547058,-3.5461688,0.6260698,0.1573728,-0.0360498,-0.0205348,&
                -0.0083818,-0.0076848,0.0,0.0/)
      ai_alpha =(/3.3288788,0.6303878,-10.9780518,2.1123398,-3.2301198,&
                -1.4473158,0.6141758,-0.3071608,0.1458498,0.0848308/)
      bi_alpha =(/67.5430838,0.6566888,9.9226018,3.9183988,-2.6882828,0.6433518,&
                  -0.8792398,-0.4817838,0.0300078,-0.1015118/)
      ai_theta =(/-3.9750378,-8.2808998,0.0611208,0.3906598,-0.4488778,0.120087,&
                0.0717048,-0.0699578,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
      bi_theta =(/-2.2839398,-3.5213068,1.9296668,-1.0832488,-0.3011748,0.1786648,&
                -0.1228608,0.0004808,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
    endif
    
    call fseries_eval(time,phi,phi_dt,Insect%a0_phi,ai_phi,bi_phi)
    call fseries_eval(time,alpha,alpha_dt,Insect%a0_alpha,ai_alpha,bi_alpha)
    call fseries_eval(time,theta,theta_dt,Insect%a0_theta,ai_theta,bi_theta)
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)

  case ("Drosophila_hovering_sun")
    !---------------------------------------------------------------------------
    ! 6th wingbeat from Chen & Sun (2014)
    ! Fourier coefficients analyzed with matlab
    ! Note that it begins with UPSTROKE, unlike Fry's kinematics
    !---------------------------------------------------------------------------
    if (.not.allocated(ai_alpha)) then
      Insect%nfft_alpha = 10
      Insect%nfft_theta = 10
      Insect%nfft_phi   = 10
      allocate ( ai_phi(1:Insect%nfft_phi), bi_phi(1:Insect%nfft_phi) ) 
      allocate ( ai_alpha(1:Insect%nfft_alpha), bi_alpha(1:Insect%nfft_alpha) )
      allocate ( ai_theta(1:Insect%nfft_theta), bi_theta(1:Insect%nfft_theta) )
      
      Insect%a0_phi   =38.2280144124915
      Insect%a0_alpha =1.09156750841542
      Insect%a0_theta =-17.0396438317138
      ai_phi   =(/-60.4766838884452,-5.34194400577534,-2.79100982711466,&
                  0.334561891664093,-0.0202655263017256,-0.323801616724358,&
                  -0.474435962657283,-0.111655938200879,0.00958151551130752,&
                  0.119666224142596/)
      bi_phi   =(/3.00359559936894,-7.55184535084148,-1.32520563461209,&
                  -0.297351445375239,-0.213108013812305,-0.0328282543472566,&
                  0.0146299151981855,-0.0385423658663155,-0.512411386850196,&
                  0.0785978606299901/)
      ai_alpha =(/-7.73228268249956,8.09409174482393,3.98349294858406,&
                  6.54460609657175,4.20944598804824,0.138380341939039,&
                  1.38813149742271,0.625930107014395,0.607953761451392,&
                  -0.688049862096416/)
      bi_alpha =(/-52.2064112743351,-3.83568699799253,-14.8200023306913,&
                  4.57428035662431,1.01656074807431,-0.113387395332322,&
                  0.614350733080735,0.637197524906845,0.878128257565861,&
                  0.271333646229075/)
      ai_theta =(/-0.0876540586926952,-6.32825811106895,0.461710021840119,&
                  -0.196290365124456,0.266829219535534,0.191278837298358,&
                  0.0651359677824226,0.132751714936873,0.104342707251547,&
                  0.0251049936194829/)  
      bi_theta =(/5.05944308805575,-0.677547202362310,-1.30945644385840,&
                  -0.741962147111828,-0.351209215835472,-0.0374224119537382,&
                  -0.0940160961803885,-0.0563224030429001,0.0533369476976694,&
                0.0507212428142968/)
    endif
    
    call fseries_eval(time,phi,phi_dt,Insect%a0_phi,ai_phi,bi_phi)
    call fseries_eval(time,alpha,alpha_dt,Insect%a0_alpha,ai_alpha,bi_alpha)
    call fseries_eval(time,theta,theta_dt,Insect%a0_theta,ai_theta,bi_theta)
    
    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    
    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)
  case ("Drosophila_hovering_maeda")
    !---------------------------------------------------------------------------
    ! Drosophila hovering kinematics protocol 
    !
    ! Fourier coefficients provided by Maeda
    ! Diditized from Fry et al.
    !---------------------------------------------------------------------------
    a_posi = (/  0.22700d0,  1.24020d0,  0.03610d0, -0.00360d0/)
    b_posi = (/  0.00000d0,  0.08880d0, -0.07000d0,  0.01250d0/)
    a_elev = (/  0.16125d0,  0.06750d0,  0.14500d0,  0.00540d0/)
    b_elev = (/  0.00000d0,  0.03670d0,  0.06840d0, -0.03390d0/)
    a_feth = (/ -0.00864d0, -0.04890d0, -0.02056d0,  0.19649d0/)
    b_feth = (/  0.00000d0, -1.17586d0, -0.01216d0, -0.17590d0/)

    ! Initialize angles and velocities
    posi = 0.0d0
    elev = 0.0d0
    feth = 0.0d0
    posi_dt = 0.0d0
    elev_dt = 0.0d0
    feth_dt = 0.0d0

    do i=0,3 !! Fourier series
      !! time dependent angle
      angles  = 2.0d0*dble(i)*pi*time

      selectcase( i )
      case( 0 ) ! Fourier 0th order
        ! mean
        dangle_posi = a_posi(1) ! +shift_mean_posi_
        dangle_elev = a_elev(1) ! +shift_mean_elev_
        dangle_feth = a_feth(1) ! +shift_mean_feth_

        dangle_posi_dt = 0.0d0
        dangle_elev_dt = 0.0d0
        dangle_feth_dt = 0.0d0
  
      case default !! Fourier n-th orders
  
        call get_dangle( &
            & angles, &                !! intent(in)
            & i, &                     !! intent(in)
            & a_posi(i+1), & !! intent(in)
            & b_posi(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_posi, &      !! intent(out)
            & dangle_posi_dt &   !! intent(out)
        & )
  
        call get_dangle( &
            & angles, &                !! intent(in
            & i, &                     !! intent(in)
            & a_elev(i+1), & !! intent(in)
            & b_elev(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_elev, &      !! intent(out)
            & dangle_elev_dt &   !! intent(out)
        & )
  
        call get_dangle( &
            & angles, &                !! intent(in
            & i, &                     !! intent(in)
            & a_feth(i+1), & !! intent(in)
            & b_feth(i+1), & !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & 0.0d0, &                 !! intent(in)
            & dangle_feth, &      !! intent(out)
            & dangle_feth_dt &   !! intent(out)
        & )
  
      endselect
  
      posi = posi +dangle_posi
      elev = elev +dangle_elev
      feth = feth +dangle_feth

      posi_dt = posi_dt +dangle_posi_dt
      elev_dt = elev_dt +dangle_elev_dt
      feth_dt = feth_dt +dangle_feth_dt
    enddo
 
    ! Convert to FLUSI's variables
    phi = posi
    alpha = -feth
    theta = -elev
    
    phi_dt = posi_dt
    alpha_dt = -feth_dt
    theta_dt = -elev_dt
    
  case ("flapper_sane")
    !---------------------------------------------------------------------------
    ! motion protocol from Sane and Dickinson, JEB 204, 2607-2626 (2001)
    !
    ! feathering: fourier coefficients analyzed with matlab, 2nd order
    !             Butterworth filter with cutoff at k=10
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 2 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in JEB 204, p. 2613 
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.807554373967804d0,&
     0.0d0,11.14661083909663d0,0.0d0,2.242734216805251d0,&
     0.0d0,-0.6141899985692184d0,0.0d0,-0.7426551158681146d0,&
     0.0d0,-0.2329560587573768d0,0.0d0,0.038749678276091284d0,&
     0.0d0,0.07083462320831221d0,0.0d0,0.028982501947490313d0,&
     0.0d0,-0.0025202918494477244d0,0.0d0,-0.010221019942802941d0,&
     0.0d0,-0.005614021318470698d0,0.0d0,1.1958884364596903d-6,&
     0.0d0,0.002186832241254999d0,0.0d0,0.0015347995090793172d0/)

    alpha = 0.0
    alpha_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of JEB 204 (eg alphedeg=90-50 for Fig 3D)
    alphacdeg = 90.0d0 - 00.0d0
    alpha = alphacdeg/40.0d0 * alpha
    alpha_dt = alphacdeg/40.0d0 * alpha_dt
   
    ! convert in radians 
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)
    
    ! *** II. position ***
    ai_phi_flapper =(/72.96795908179631d0,&
     0.0d0,8.064401876272864d0,0.0d0,2.769062401215844d0,&
     0.0d0,1.2200252377066352d0,0.0d0,0.5584689705779989d0,&
     0.0d0,0.2545617536476344d0,0.0d0,0.11829515180579572d0,&
     0.0d0,0.05754453975774996d0,0.0d0,0.02964141751269772d0,&
     0.0d0,0.016177705089515895d0,0.0d0,0.009315101869467001d0,&
     0.0d0,0.005625663922446026d0,0.0d0,0.0035424425357352385d0,&
     0.0d0,0.0023130422432356247d0,0.0d0,0.001558278163264511d0,&
     0.0d0,0.001078213692334021d0/)

    phi = 0.0
    phi_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo
    
    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0d0
    phi = phicdeg/180.0d0 * phi
    phi_dt = phicdeg/180.0d0 * phi_dt
   
    ! convert in radians 
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)
    
    ! *** III. elevation ***
    theta = 0.0d0
    theta_dt = 0.0d0

  case ("flapper_dickinson")
    !---------------------------------------------------------------------------
    ! motion protocol from Dickinson, Lehmann and Sane, Science (1999)
    !
    ! feathering: fourier coefficients analyzed with matlab,
    !             Gaussian filter that fits fig 3D
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 5 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in Science
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.23094285611071d0,&
      0.0d0,10.224154661301371d0,0.0d0,2.1623763046726396d0,&
      0.0d0,0.05049394424178093d0,0.0d0,-0.17550942623071494d0,&
      0.0d0,-0.06634193748204852d0,0.0d0,-0.008925020495896451d0,&
      0.0d0,0.0011292567942149407d0,0.0d0,6.471071566666472d-4,&
      0.0d0,1.0018757795834964d-4,0.0d0,3.0105550216312524d-6,&
      0.0d0,-1.237567150768195d-6,0.0d0,-1.988004402010933d-7,&
      0.0d0,-1.10165545174181d-8,0.0d0,2.4135650975460306d-10/)

    ! Advanced rotation (+ sign) or delayed rotation (- sign)
    tadv = 0.08
    !tadv = - 0.08

    alpha = 0.0
    alpha_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+tadv)) 
      c = dcos(f*dble(i)*(time+tadv))
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of Science (eg alphedeg=90-40 for Fig 3)
    alphacdeg = 90.0d0 - 40.0d0
    alpha = alphacdeg/40.0d0 * alpha
    alpha_dt = alphacdeg/40.0d0 * alpha_dt
   
    ! convert in radians 
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)
    
    ! *** II. position ***
    ai_phi_flapper =(/63.24528806534019d0,&
      0.0d0,5.753991800610726d0,0.0d0,1.3887974015525626d0,&
      0.0d0,0.3889856512386744d0,0.0d0,0.10577402496901325d0,&
      0.0d0,0.026061339604144987d0,0.0d0,0.005623376646981709d0,&
      0.0d0,0.001042285996467963d0,0.0d0,1.639611509380189d-4,&
      0.0d0,2.1716252827442023d-5,0.0d0,2.408190194815521d-6,&
      0.0d0,2.2268710288534648d-7,0.0d0,1.7118916093759426d-8,&
      0.0d0,1.0914870312823793d-9,0.0d0,5.76135101855556d-11,&
      0.0d0,2.513944479978149d-12/)

    phi = 0.0
    phi_dt = 0.0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo
    
    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0d0
    phi = phicdeg/180.0d0 * phi
    phi_dt = phicdeg/180.0d0 * phi_dt
   
    ! convert in radians 
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0d0
    theta_dt = 0.0d0

  case ("flapper_ramamurti")
    !---------------------------------------------------------------------------
    ! motion protocol from Ramamurti and Sandberg, JEB (2002)
    !
    ! feathering: fourier coefficients analyzed with matlab,
    !             Gaussian filter that fits fig 2
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 10 Jun 2013
    !---------------------------------------------------------------------------

    ! Time shift as in JEB
    t0adv = 0.28175d0

    ! *** I. feathering motion ***
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
!    bi_alpha_flapper =(/48.23094285611071d0,&
!      0.0d0,10.224154661301371d0,0.0d0,2.1623763046726396d0,&
!      0.0d0,0.05049394424178093d0,0.0d0,-0.17550942623071494d0,&
!      0.0d0,-0.06634193748204852d0,0.0d0,-0.008925020495896451d0,&
!      0.0d0,0.0011292567942149407d0,0.0d0,6.471071566666472d-4,&
!      0.0d0,1.0018757795834964d-4,0.0d0,3.0105550216312524d-6,&
!      0.0d0,-1.237567150768195d-6,0.0d0,-1.988004402010933d-7,&
!      0.0d0,-1.10165545174181d-8,0.0d0,2.4135650975460306d-10/)

    bi_alpha_flapper =(/58.5622945117485d0,&
      0.0d0,9.70856389020196d0,0.0d0,1.06772463979698d0,&
      0.0d0,-0.127455709998572d0,0.0d0,-0.0553559839123380d0,&
      0.0d0,-0.00500430568361157d0,0.0d0,0.000181637429802491d0,&
      0.0d0,5.24470001944736d-05,0.0d0,2.39488192125763d-06,&
      0.0d0,-1.62462115890220d-08,0.0d0,-3.60732923972996d-09,&
      0.0d0,-7.58943924434529d-11,0.0d0,-3.17275659750908d-16,&
      0.0d0,1.50381167957580d-14,0.0d0,1.41146073296091d-16/)

    ! Advanced rotation (+ sign) or delayed rotation (- sign)
    tadv = 0.08d0
    !tadv = - 0.08d0

    alpha = 0.0d0
    alpha_dt = 0.0d0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+tadv+t0adv)) 
      c = dcos(f*dble(i)*(time+tadv+t0adv))
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! convert in radians 
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)
    
    ! *** II. position ***
!    ai_phi_flapper =(/63.24528806534019d0,&
!      0.0d0,5.753991800610726d0,0.0d0,1.3887974015525626d0,&
!      0.0d0,0.3889856512386744d0,0.0d0,0.10577402496901325d0,&
!      0.0d0,0.026061339604144987d0,0.0d0,0.005623376646981709d0,&
!      0.0d0,0.001042285996467963d0,0.0d0,1.639611509380189d-4,&
!      0.0d0,2.1716252827442023d-5,0.0d0,2.408190194815521d-6,&
!      0.0d0,2.2268710288534648d-7,0.0d0,1.7118916093759426d-8,&
!      0.0d0,1.0914870312823793d-9,0.0d0,5.76135101855556d-11,&
!      0.0d0,2.513944479978149d-12/)

    a0_phi = 20.0d0
    ai_phi_flapper =(/71.1023524748246d0,&
      0.0d0,6.43355277659369d0,0.0d0,1.53592591540135d0,&
      0.0d0,0.423188982190934d0,0.0d0,0.112580047175950d0,&
      0.0d0,0.0269873997935609d0,0.0d0,0.00563410946130459d0,&
      0.0d0,0.00100469248408058d0,0.0d0,0.000151188714225520d0,&
      0.0d0,1.90437108218227d-05,0.0d0,1.99627349027768d-06,&
      0.0d0,1.73399283714958d-07,0.0d0,1.24381425530736d-08,&
      0.0d0,7.34711932879597d-10,0.0d0,3.56493134534579d-11,&
      0.0d0,1.41757814786810d-12/)

    phi = 0.5d0 * a0_phi
    phi_dt = 0.0d0
    
    ! frequency factor
    f = 2.d0*pi
    
    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+t0adv)) 
      c = dcos(f*dble(i)*(time+t0adv))
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo
    
    ! convert in radians 
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0d0
    theta_dt = 0.0d0
   
  case ("takeoff")
    !--------------------------------------------------
    ! Fontaine et al. 
    !--------------------------------------------------
    if (Insect%KineFromFile/="no") then
      call wing_kine_interp(time,phi,alpha,theta,phi_dt,alpha_dt,theta_dt)
      ! position angle
      phi = deg2rad(phi)
      phi_dt = deg2rad(phi_dt)
      ! feathering angle
      alpha = deg2rad(alpha)         
      alpha_dt = deg2rad(alpha_dt)    
      ! elevation angle in flusi coordinates
      theta = -theta
      theta_dt = - theta_dt
      theta = deg2rad(theta)
      theta_dt = deg2rad(theta_dt)
    endif 
 
  case ("simplified")
    !---------------------------------------------------------------------------
    ! simplified motion protocoll
    !
    ! J. comput. Phys. 231 (2012) 1822-1847 "A fluid-structure interaction 
    ! model of insect flight with flexible wings"
    !
    ! the pase shift "phase" was my idea
    !---------------------------------------------------------------------------
    phi_max     = 60.d0*pi/180.d0  ! phi is up/down angle (flapping)
    alpha_max   = 0.d0!45.d0*pi/180.d0  ! alpha is tethering
    phase       = 0.d0! 10.d0*pi/180.d0  ! phase shift between flapping and tethering
    f = 1.d0*2.0*pi
    
    phi      = phi_max  *dcos(f*time)
    alpha    = alpha_max*dsin(f*(time+phase))
    theta    = 0.0
    
    phi_dt   =-phi_max *f *dsin(f*time)
    alpha_dt = alpha_max*f*dcos(f*(time+phase))
    theta_dt = 0.0
  case ("debug")
    phi      = deg2rad(45.d0)   
    alpha    = deg2rad(0.d0)
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case ("none")
    phi      = 0.0
    alpha    = 0.0
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case default
    if (mpirank==0) then
    write (*,*) "insects.f90::FlappingMotion: motion case (protocoll) undefined"
    call abort()
    endif    
  end select
  
end subroutine FlappingMotion

!-------------------------------------------------------------------------------
! evaluate a fourier series given by the coefficents a0,ai,bi
! at the time "time", return the function value "u" and its
! time derivative "u_dt". Uses assumed-shaped arrays, requires an interface.
!-------------------------------------------------------------------------------
subroutine fseries_eval(time,u,u_dt,a0,ai,bi)
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: a0, time
  real(kind=pr), intent(in), dimension(:) :: ai,bi
  real(kind=pr), intent(out) :: u, u_dt
  real(kind=pr) :: c,s,f
  integer :: nfft, i
  
  nfft=size(ai)
  
  ! frequency factor
  f = 2.d0*pi

  u = 0.5d0*a0  
  u_dt = 0.d0
  
  do i=1,nfft
      s = dsin(f*dble(i)*time) 
      c = dcos(f*dble(i)*time)
      ! function value
      u    = u + ai(i)*c + bi(i)*s
      ! derivative (in time)
      u_dt = u_dt + f*dble(i)*(-ai(i)*s + bi(i)*c)
  enddo
end subroutine fseries_eval


!-------------------------------------------------------------------------------
! evaluate hermite series, given by coefficients ai (function values)
! and bi (derivative values) at the locations x. Note that x is assumed periodic;
! do not include x=1.0.
! a valid example is x=(0:N-1)/N
!-------------------------------------------------------------------------------
subroutine hermite_eval(time,u,u_dt,ai,bi)
  use vars
  implicit none
  
  real(kind=pr), intent(in) :: time
  real(kind=pr), intent(in), dimension(:) :: ai,bi
  real(kind=pr), intent(out) :: u, u_dt
  real(kind=pr) :: dt,h00,h10,h01,h11,t
  integer :: n, j1,j2
  
  n=size(ai)
  
  dt = 1.d0 / dble(n)
  j1 = floor(time/dt) + 1
  j2 = j1 + 1
  ! periodization
  if (j2 > n) j2=1
  ! normalized time (between two data points)
  t = (time-dble(j1-1)*dt) /dt
  
  ! values of hermite interpolant
  h00 = (1.d0+2.d0*t)*((1.d0-t)**2)
  h10 = t*((1.d0-t)**2)
  h01 = (t**2)*(3.d0-2.d0*t)
  h11 = (t**2)*(t-1.d0)
  
  ! function value
  u = h00*ai(j1) + h10*dt*bi(j1) &
    + h01*ai(j2) + h11*dt*bi(j2)
  
  ! derivative values of basis functions
  h00 = 6.d0*t**2 - 6.d0*t
  h10 = 3.d0*t**2 - 4.d0*t + 1.d0
  h01 =-6.d0*t**2 + 6.d0*t 
  h11 = 3.d0*t**2 - 2.d0*t
  
  ! function derivative value
  u_dt = (h00*ai(j1) + h10*dt*bi(j1) &
       + h01*ai(j2) + h11*dt*bi(j2) ) / dt
end subroutine hermite_eval
