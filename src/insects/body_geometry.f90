subroutine draw_body( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none
  
  type(diptera),intent(inout) :: Insect  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)
  
  real(kind=pr) :: t1
  t1 = MPI_wtime()
  
  select case (Insect%BodyType)
  case ("nobody")
    return
  case ("jerry")
    call draw_body_jerry( mask, mask_color, us, Insect, color_body, M_body)
  case ("drosophila")
    call draw_body_drosophila( mask, mask_color, us, Insect, color_body, M_body)
  case ("drosophila_maeda","drosophila_slim")
    call draw_body_drosophila_maeda( mask, mask_color, us, Insect, color_body, M_body)
  case default
    if(mpirank==0) write(*,*) "Insect::draw_body::Insect%BodyType unknown..."
    if(mpirank==0) write(*,*) Insect%bodytype
    call abort()
  end select

  time_insect_body = time_insect_body + MPI_wtime() - t1
end subroutine


!-------------------------------------------------------------------------------

! Body adapted from Maeda & Liu. It assumes Insect%x_head=0.0
subroutine draw_body_drosophila_maeda( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none
  
  type(diptera),intent(inout) :: Insect  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  integer :: ix,iy,iz
  real(kind=pr) :: x,y,z,s,s1, a_body, R,R0,R_tmp,x1
  real(kind=pr) :: x_glob(1:3),x_body(1:3),x_head(1:3)
  real(kind=pr) :: rbc,thbc1,thbc2,x0bc,z0bc,xcs,zcs
  real(kind=pr) :: xx_head,zz_head,dx_head,dz_head,a_head
  
  
  !-----------------------------------------------------------------------------
  ! Body
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
     do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
           ! x_glob is in the global coordinate system
           x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
           ! x_body is in the body coordinate system
           x_body = matmul(M_body,x_glob-Insect%xc_body)
           
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
            rbc = 0.9464435146443515d0
            thbc1 = 112.0d0 *pi/180.0d0
            thbc2 = 53.0d0 *pi/180.0d0
            x0bc = -0.24476987447698745d0
            z0bc = -0.9301255230125524d0
          
            ! chordwise dimensionless coordinate, from head to abdomen
            s = (datan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1) 
            ! body center coordinates at s
            xcs = x0bc + (x-x0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)
            zcs = z0bc + (z-z0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)

            ! check if inside body bounds (in s-direction)
            if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then    
              R0 = 0.0d0
              ! round section by default
              if (Insect%BodyType == 'drosophila_slim') then 
                a_body = 1.09d0
              else
                a_body = 1.0d0
              endif
              ! distortion of s
              s1 = 1.0d0 - ( s + 0.08d0*dtanh(30.0d0*s) ) / (1.0d0+0.08d0*dtanh(30.0d0))
              s1 = ( s1 + 0.04d0*dtanh(60.0d0*s1) ) / (1.0d0+0.04d0*dtanh(60.0d0))
        !       s1 = ( max(dsin(1.2d0*s1)/dsin(1.2d0), 0.d0) )**1.25  
              s1 = ( dsin(1.2d0*s1)/dsin(1.2d0) )**1.25
              
              x1 = 1.075d0 * s1 
              ! compute radius as a function of x1 (counting from the tail on)
              ! same shape as 'drosophila'
              if (x1 < 0.6333d0) then
                ! we're in the ABDOMEN
                R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
              elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
                ! we're in the THORAX 
                R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
                ! slim body
                if (Insect%BodyType == 'drosophila_slim') &
                  a_body = 1.09d0-0.19d0*(x1-0.6333d0)*(x1-1.075d0)/0.0488d0
              endif
              ! distortion of R0
              R0 = 0.8158996d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0
              ! distance to the body center at s
              R = dsqrt( (x-xcs)**2 + (a_body*y)**2 + (z-zcs)**2 )

              ! smoothing
              if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
                R_tmp = steps(R,R0)        
                mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
                mask_color(ix,iy,iz) = color_body
              endif      
            
            endif
        enddo
     enddo
  enddo  
  
  !-----------------------------------------------------------------------------
  ! Head
  !-----------------------------------------------------------------------------
  if (Insect%BodyType == 'drosophila_slim') then
    a_head = 1.09d0
  else
    a_head = 1.0d0
  endif
  
  ! ellipsoid head, assumes xc_head=0 in .ini file
  xx_head = 0.17d0
  zz_head = -0.1d0
  dx_head = 0.5d0 * 0.185d0
  dz_head = 0.5d0 * 0.27d0

  do iz = ra(3), rb(3)
     do iy = ra(2), rb(2)
        do ix = ra(1), rb(1)
          !-- define the head coordinate systems we are going to use
          x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
          x_body   = matmul(M_body,x_glob-Insect%xc_body)
          x_head   = x_body - Insect%x_head
          
          ! check if inside the surrounding box (save comput. time)
          if ( dabs(x_head(2)) <= dz_head + Insect%safety ) then
            if ( dabs(x_head(3)-zz_head) <= dz_head + Insect%safety ) then
              ! check for length inside ellipsoid:
              if ( dabs(x_head(1)-xx_head) < dx_head + Insect%safety ) then
              
                R  = dsqrt ( (a_head*x_head(2))**2 + (x_head(3)-zz_head)**2 )
                ! this gives the R(x) shape
                if ( ((x_head(1)-xx_head)/dx_head)**2 <= 1.d0) then
                  R0 = dz_head*dsqrt(1.d0- ((x_head(1)-xx_head)/dx_head)**2 )
                  if ( R < R0 + Insect%safety ) then
                    mask(ix,iy,iz)= max(steps(R,R0),mask(ix,iy,iz))
                    mask_color(ix,iy,iz) = color_body
                  endif
                endif
              endif
            endif
          endif
           
           
        enddo
     enddo
  enddo  
end subroutine draw_body_drosophila_maeda

!-------------------------------------------------------------------------------

! frist attempt drosophila body ("doro"). It is a straight body (center line
! undistorted), and its radius is a single function, defined in thorax
! abdomen and head
subroutine draw_body_drosophila( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none
  
  type(diptera),intent(inout) :: Insect  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)
  
  real(kind=pr) :: x,R0,R,R_tmp,x_tmp
  real(kind=pr) :: x_body(1:3), x_glob(1:3)
  integer :: ix,iy,iz
  
  
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob-Insect%xc_body)
        
        ! ------------------------------------
        ! two b-splines body (abdomen+thorax)
        ! ------------------------------------    
        x = x_body(1) + 0.8067d0 ! centers the thickest part of the thorax at the origin
        
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
            mask_color(ix,iy,iz) = color_body
          endif      
        
        endif  
      enddo
    enddo
  enddo
  
end subroutine draw_body_drosophila

!-------------------------------------------------------------------------------

subroutine draw_body_jerry( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none
  
  type(diptera),intent(inout) :: Insect  
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)
  
  real(kind=pr) :: x,R0,R,R_tmp,x_tmp,a_body
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
  integer :: ix,iy,iz
  
  ! overwrite
  Insect%R_head = 0.125d0
  Insect%R_eye = 0.0625d0
  Insect%x_pivot_r =(/ 0.05d0, -0.2165d0, 0.d0 /)
  Insect%x_pivot_l =(/ 0.05d0, +0.2165d0, 0.d0 /)
  Insect%b_body = 0.1
  Insect%L_body = 1.0
  Insect%x_head = (/0.5*Insect%L_body,0.d0,0.d0 /)
  Insect%x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
                   *0.8d0*(/1.d0,+1.d0,1.d0/)
  Insect%x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
                   *0.8d0*(/1.d0,-1.d0,1.d0/)
  
  !-----------------------------------------------------------------------------
  ! ellipsoid body
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob-Insect%xc_body)
        
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
                    mask_color(ix,iy,iz) = color_body
                  endif
                endif
            endif
          endif
        endif
      enddo
    enddo
  enddo
  
  !-----------------------------------------------------------------------------
  ! spherical head
  !-----------------------------------------------------------------------------
  x_head = Insect%xc_body + matmul(transpose(M_body),Insect%x_head)
  call drawsphere( x_head,Insect%R_head,mask,mask_color,us,Insect,color_body )
  
  x_eye = Insect%xc_body + matmul(transpose(M_body),Insect%x_eye_l)
  call drawsphere( x_eye,Insect%R_eye,mask,mask_color,us,Insect,color_body )
  
  x_eye = Insect%xc_body + matmul(transpose(M_body),Insect%x_eye_r)
  call drawsphere( x_eye,Insect%R_eye,mask,mask_color,us,Insect,color_body )
end subroutine draw_body_jerry

!-------------------------------------------------------------------------------

! draw a sphere with radius R0 centered at xc. as part of the body, it has no
! velocity vector (this is added to the entire insect later)
subroutine drawsphere( xc,R0,mask,mask_color,us,Insect,icolor )
  use vars
  implicit none
  
  real(kind=pr),intent(in)::xc(1:3),R0
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: icolor
  type(diptera),intent(inout) :: Insect
  
  integer :: ix,iy,iz
  real(kind=pr)::x(1:3),R
  
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x is in the global coordinate system
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x is now centered in the sphere's center point
        x = x - xc
        
        if (dabs(x(1)) <= R0+Insect%safety) then
          if (dabs(x(2)) <= R0+Insect%safety) then
            if (dabs(x(3)) <= R0+Insect%safety) then
              R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
              if ( R <= R0+Insect%safety ) then
                mask(ix,iy,iz) = max(steps(R,R0),mask(ix,iy,iz))
                mask_color(ix,iy,iz) = icolor
              endif
            endif
          endif
        endif 
        
      enddo
    enddo
  enddo  
  
end subroutine 