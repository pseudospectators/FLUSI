!------------------------------------------------------------------------------
! Draws an insect's body, several options available.
! the body is, in the local coordinate system, always aligned with the
! x-axis. also, we currently use only rotational symmetric bodies.
subroutine DrawBody(ix,iy,iz,x_body)
  use fsi_vars
  use mpi
  implicit none
  real(kind=pr) :: a_body, R, R0, steps, x, y, z, s, s1, x1, x_tmp, R_tmp
  real(kind=pr) :: rbc, x0bc, z0bc, thbc1, thbc2, xcs, zcs
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x_body(1:3)
  
  select case (Insect%BodyType)
  case ('ellipsoid')  
    ! ------------------------------------
    ! ellipsoid body (jerry)
    ! ------------------------------------
    a_body = Insect%L_body / 2.d0
    
    ! check if inside the surrounding box (save comput. time)
    if ( dabs(x_body(2)) <= Insect%b_body + Insect%safety ) then
    if ( dabs(x_body(3)) <= Insect%b_body + Insect%safety ) then
    
    ! check for length inside ellipsoid:
    if ( dabs(x_body(1) ) < Insect%L_body/2 + Insect%safety ) then
        R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
        ! this gives the R(x) shape
        if ( (x_body(1)/a_body)**2 <= 1.d0) then
        R0 = dsqrt( Insect%b_body**2 *(1.d0- (x_body(1)/a_body)**2 ) )

        if ( R < R0 + Insect%safety ) then
          mask(ix,iy,iz)= max(steps(R,R0),mask(ix,iy,iz))
          ! body has the color "2"
          mask_color(ix,iy,iz) = 2
        endif
        endif
    endif
    endif
    endif
    
    
  case ('drosophila')
    ! ------------------------------------
    ! two b-splines body (abdomen+thorax)
    ! ------------------------------------    
    x = x_body(1) + 0.8067 ! centers the thickest part of the thorax at the origin
    
    ! check if inside body bounds (in x-direction)
    if ( (x>=-Insect%safety) .and. (x<=1.2+Insect%safety) ) then    
      R0=0.0
      ! compute radius as a function of x (counting from the tail on)
      if (x < 0.6333) then
        ! we're in the ABDOMEN
        R0 = max( -1.2990*x**2 + 0.9490*x + 0.0267, 0.d0)
      elseif ((x >= 0.6333) .and. (x <=1.0 )) then
        ! we're in the THORAX 
        R0 = max( -2.1667*x**2 + 3.4661*x - 1.2194, 0.d0)
      elseif ((x >= 1.0) .and. (x <=1.2 )) then
        ! we're in the HEAD
        R0 = max( -12.68*x**2 + 27.4960*x - 14.7360, 0.d0)
      endif
    
      ! radius at this point
      R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
      
      ! smoothing in x-direction
      if (x<Insect%safety) then  ! xs is chordlength coordinate
        x_tmp = steps(-x, Insect%smooth)
      else
        x_tmp = steps( x,1.2-Insect%smooth)
      endif     
      
      
      if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
        R_tmp = steps(R,R0)        
        mask(ix,iy,iz)= max( R_tmp*x_tmp , mask(ix,iy,iz) )
        ! body has the color "2"
        mask_color(ix,iy,iz) = 2
      endif      
    
    endif
    
  
  
  case ('drosophila_maeda')
    ! ------------------------------------
    ! approximation to mesh from Maeda
    ! similar to Aono et al.
    ! ------------------------------------    
    x = x_body(1)
    y = x_body(2)
    z = x_body(3)

    ! symmetry plane is xz
    ! +x direction is forward
    ! body centerline is an arc with center at x0bc,y0bc
    ! radius rbc and angles th counted from negative z
    rbc = 0.9464435146443515
    thbc1 = 112.0d0 *pi/180.0d0
    thbc2 = 53.0d0 *pi/180.0d0
    x0bc = -0.24476987447698745d0
    z0bc = -0.9301255230125524
  
    ! chordwise dimensionless coordinate, from head to abdomen
    s = (atan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1) 
    ! body center coordinates at s
    xcs = x0bc + (x-x0bc)*rbc/sqrt((x-x0bc)**2+(z-z0bc)**2)
    zcs = z0bc + (z-z0bc)*rbc/sqrt((x-x0bc)**2+(z-z0bc)**2)
    ! distance to the body center at s
    R = sqrt( (x-xcs)**2 + y**2 + (z-zcs)**2 )

    ! check if inside body bounds (in s-direction)
    if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then    
      R0 = 0.0d0
      ! distortion of s
      s1 = 1.0d0 - ( s + 0.08d0*tanh(30.0d0*s) ) / (1.0d0+0.08d0*tanh(30.0d0))
      s1 = ( s1 + 0.04d0*tanh(60.0d0*s1) ) / (1.0d0+0.04d0*tanh(60.0d0))
      s1 = (sin(1.2d0*s1)/sin(1.2d0))**1.25d0
      x1 = 1.075d0 * s1 
      ! compute radius as a function of x1 (counting from the tail on)
      ! same shape as 'drosophila'
      if (x1 < 0.6333d0) then
        ! we're in the ABDOMEN
        R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
      elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
        ! we're in the THORAX 
        R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
      endif
      ! distortion of R0
      R0 = 0.8158996d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0

      ! smoothing
      if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
        R_tmp = steps(R,R0)        
        mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
        ! body has the color "2"
        mask_color(ix,iy,iz) = 2
      endif      
    
    endif

    
  case ('nobody')
    ! doesn't do anything
  case default
    if (mpirank==0) then
    write (*,*) "In DrawBody: unkown body type."
    stop
    endif
  end select
  
end subroutine




!------------------------------------------------------------------------------
! Draws a sphere with radius R, as we need for the head and the eyes, if present
subroutine DrawSphere(ix,iy,iz,x,R0)
  use fsi_vars
  use mpi
  implicit none
  real(kind=pr), intent(in) :: R0
  real(kind=pr) :: R, steps
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  
  if (abs(x(1))<R0+Insect%safety) then
  if (abs(x(2))<R0+Insect%safety) then
  if (abs(x(3))<R0+Insect%safety) then
      R = sqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
      if ( R <= R0+Insect%safety ) then
        mask(ix,iy,iz) = max(steps(R,R0),mask(ix,iy,iz))
        ! body parts have the color "2"
        mask_color(ix,iy,iz) = 2
      endif
  endif
  endif
  endif  
end subroutine

! as long as we have only spherical eyes, this is just a wrapper
subroutine DrawEye(ix,iy,iz,x)
  use fsi_vars
  use mpi
  implicit none
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  if (Insect%HasEye=="yes") then
  call DrawSphere(ix,iy,iz,x,Insect%R_eye)
  endif
end subroutine

! as long as we have only spherical heads, this is just a wrapper
subroutine DrawHead(ix,iy,iz,x)
  use fsi_vars
  use mpi
  implicit none
  integer, intent(in) :: ix,iy,iz
  real(kind=pr),intent(in) :: x(1:3)
  real(kind=pr) :: x_head,z_head,dx_head,dz_head,R,R0,steps
  if (Insect%HasHead=="yes") then
  
     select case (Insect%BodyType)
     case ('ellipsoid')  
      ! an ellipsoid body goes with a spherical head
      call DrawSphere(ix,iy,iz,x,Insect%R_head)

     case ('drosophila')
      ! drosophilae have different heads.

     case ('drosophila_maeda')  
      ! ellipsoid head, assumes xc_head=0 in .ini file
      x_head = 0.17d0
      z_head = -0.1d0
      dx_head = 0.5d0*  0.185d0
      dz_head = 0.5d0*  0.27d0
      ! check if inside the surrounding box (save comput. time)
      if ( dabs(x(2)) <= dz_head + Insect%safety ) then
      if ( dabs(x(3)-z_head) <= dz_head + Insect%safety ) then
      ! check for length inside ellipsoid:
      if ( dabs(x(1)-x_head) < dx_head + Insect%safety ) then
        R  = dsqrt ( x(2)**2 + (x(3)-z_head)**2 )
        ! this gives the R(x) shape
        if ( ((x(1)-x_head)/dx_head)**2 <= 1.d0) then
        R0 = dz_head*dsqrt(1.d0- ((x(1)-x_head)/dx_head)**2 )
        if ( R < R0 + Insect%safety ) then
          mask(ix,iy,iz)= max(steps(R,R0),mask(ix,iy,iz))
          ! body parts have the color "2"
          mask_color(ix,iy,iz) = 2
        endif
        endif
      endif
      endif
      endif

     case default
      ! do nothing
     end select 
  endif
end subroutine