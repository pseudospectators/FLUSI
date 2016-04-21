subroutine draw_body( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: t1
  t1 = MPI_wtime()


  select case (Insect%BodyType)
  case ("nobody")
    return
  case ("suzuki_thin_rod")
    call draw_suzuki_thin_rod( mask, mask_color, us, Insect, color_body, M_body)
  case ("jerry","Jerry")
    call draw_body_jerry( mask, mask_color, us, Insect, color_body, M_body)
  case ("hawkmoth","Hawkmoth")
    call draw_body_hawkmoth( mask, mask_color, us, Insect, color_body, M_body)
  case ("particle")
    call draw_body_particle( mask, mask_color, us, Insect, color_body, M_body)
  case ("platicle")
    call draw_body_platicle( mask, mask_color, us, Insect, color_body, M_body)
  case ("coin")
    call draw_body_coin( mask, mask_color, us, Insect, color_body, M_body)
  case ("sphere","SPHERE","Sphere")
    call draw_body_sphere( mask, mask_color, us, Insect, color_body, M_body)
  case ("drosophila")
    call draw_body_drosophila( mask, mask_color, us, Insect, color_body, M_body)
  case ("drosophila_maeda","drosophila_slim")
    call draw_body_drosophila_maeda( mask, mask_color, us, Insect, color_body, M_body)
  case ("bumblebee")
    call draw_body_bumblebee( mask, mask_color, us, Insect, color_body, M_body)
  case ("mosquito_iams")
    call draw_body_mosquito_iams( mask, mask_color, us, Insect, color_body, M_body)
  case default
    if(mpirank==0) write(*,*) "Insect::draw_body::Insect%BodyType unknown..."
    if(mpirank==0) write(*,*) Insect%bodytype
    call abort()
  end select

  time_insect_body = time_insect_body + MPI_wtime() - t1
end subroutine


!-------------------------------------------------------------------------------

! Bumblebee body, BB1 in Dudley & Ellington JEB 1990
subroutine draw_body_bumblebee( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout) :: mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout) :: us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout) :: mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  integer :: ix,iy,iz,j
  real(kind=pr) :: x,y,z,s,s1,a_body,R,R0,R_tmp,x1
  real(kind=pr) :: x_glob(1:3),x_body(1:3),x_head(1:3)
  real(kind=pr) :: rbc,thbc1,thbc2,x0bc,z0bc,xcs,zcs
  real(kind=pr) :: xx_head,zz_head,dx_head,dz_head,a_head
  real(kind=pr) :: xl1(5),yl1(5),zl1(5),rl1(4),xl2(5),yl2(5),zl2(5),rl2(4),&
  xl3(5),yl3(5),zl3(5),rl3(4),xf(2),yf(2),zf(2),rf,xan(2),yan(2),zan(2),ran,&
  xmin_bbox,xmax_bbox,ymin_bbox,ymax_bbox,zmin_bbox,zmax_bbox

  !-----------------------------------------------------------------------------
  ! Body
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob)

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
        rbc = 1.3d0
        thbc1 = 112.0d0 *pi/180.0d0
        thbc2 = 53.0d0 *pi/180.0d0
        x0bc = 0.0782301255230126d0
        z0bc = -1.26512552301255d0

        ! chordwise dimensionless coordinate, from head to abdomen
        s = (datan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1)
        ! body center coordinates at s
        xcs = x0bc + (x-x0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)
        zcs = z0bc + (z-z0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)

        ! check if inside body bounds (in s-direction)
        if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then
          R0 = 0.0d0
          ! round section by default
          a_body = 1.0d0
          ! distortion of s
          s1 = 1.0d0 - ( s + 0.08d0*dtanh(30.0d0*s) ) / (1.0d0+0.08d0*dtanh(30.0d0))
          s1 = ( s1 + 0.04d0*dtanh(60.0d0*s1) ) / (1.0d0+0.04d0*dtanh(60.0d0))
          s1 = ( dsin(1.2d0*s1)/dsin(1.2d0) )**1.25

          x1 = 1.075d0 * s1
          ! compute radius as a function of x1 (counting from the tail on)
          ! same shape as 'drosophila'
          if (x1 < 0.6333d0) then
            ! we're in the ABDOMEN
            R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
            ! flatten abdomen
            a_body = 1.0d0-0.07d0*(x1-0.6333d0)*x1/0.0488d0
          elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
            ! we're in the THORAX
            R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
          endif
          ! distortion of R0
          R0 = 1.2d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0
          ! distance to the body center at s
          R = dsqrt( (x-xcs)**2 + y**2 + (a_body*(z-zcs))**2 )

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
  a_head = 1.04d0

  ! ellipsoid head, assumes xc_head=0 in .ini file
  xx_head = 0.58125d0
  zz_head = -0.1d0
  dx_head = 0.5d0 * 0.2035d0
  dz_head = 0.5d0 * 0.297d0

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the head coordinate systems we are going to use
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        x_body   = matmul(M_body,x_glob)
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

  !-----------------------------------------------------------------------------
  ! Legs, antennae and proboscis
  !-----------------------------------------------------------------------------
  if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="legs")&
  .or.(Insect%HasDetails=="antennae_proboscis")) then
  ! Bounding box
  xmin_bbox = -0.74-Insect%safety
  xmax_bbox = 0.8+Insect%safety
  ymin_bbox = -0.35-Insect%safety
  ymax_bbox = 0.35+Insect%safety
  zmin_bbox = -0.38-Insect%safety
  zmax_bbox = 0.1+Insect%safety

  ! Parameters of legs, antennae and proboscis
  xl1 = (/-0.74,-0.63,-0.4,-0.1,0.1/)
  yl1 = (/0.32,0.32,0.31,0.3,0.12/)
  zl1 = (/-0.35,-0.37,-0.2,-0.1,-0.16/)
  rl1 = (/0.015,0.03,0.04,0.03/)*1.3
  xl2 = (/-0.24,-0.15,0.02,0.17,0.19/)
  yl2 = (/0.33,0.33,0.32,0.3,0.15/)
  zl2 = (/-0.29,-0.28,-0.2,-0.15,-0.19/)
  rl2 = (/0.015,0.03,0.04,0.03/)*1.3
  xl3 = (/0.28,0.35,0.45,0.4,0.35/)
  yl3 = (/0.31,0.30,0.28,0.2,0.15/)
  zl3 = (/-0.3,-0.28,-0.25,-0.18,-0.18/)
  rl3 = (/0.015,0.02,0.03,0.02/)*1.3
  xf = (/0.43,0.6/)
  yf = (/0.0,0.0/)
  zf = (/-0.28,-0.23/)
  rf = 0.017*1.3
  xan = (/0.63,0.8/)
  yan = (/0.05,0.27/)
  zan = (/-0.03,0.1/)
  ran = 0.015*1.3

  ! Assign values to mask pointwise
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        !-- define the head coordinate systems we are going to use
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        x_body   = matmul(M_body,x_glob)

        !-- check bounds
        if ((x_body(1)>=xmin_bbox).and.(x_body(1)<=xmax_bbox)&
        .and.(x_body(2)>=ymin_bbox).and.(x_body(2)<=ymax_bbox)&
        .and.(x_body(3)>=zmin_bbox).and.(x_body(3)<=zmax_bbox)) then

        !-- left and right legs, antennae and proboscis
        if (x_body(2)>0) then
          if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="legs")) then
            do j=1,4
              call draw_cylinder(x_body,xl1(j),yl1(j),zl1(j),&
              xl1(j+1),yl1(j+1),zl1(j+1),rl1(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
            do j=1,4
              call draw_cylinder(x_body,xl2(j),yl2(j),zl2(j),&
              xl2(j+1),yl2(j+1),zl2(j+1),rl2(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
            do j=1,4
              call draw_cylinder(x_body,xl3(j),yl3(j),zl3(j),&
              xl3(j+1),yl3(j+1),zl3(j+1),rl3(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
          endif
          if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="antennae_proboscis")) then
            call draw_cylinder(x_body,xan(1),yan(1),zan(1),&
            xan(2),yan(2),zan(2),ran,mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
          endif
        else
          if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="legs")) then
            do j=1,4
              call draw_cylinder(x_body,xl1(j),-yl1(j),zl1(j),&
              xl1(j+1),-yl1(j+1),zl1(j+1),rl1(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
            do j=1,4
              call draw_cylinder(x_body,xl2(j),-yl2(j),zl2(j),&
              xl2(j+1),-yl2(j+1),zl2(j+1),rl2(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
            do j=1,4
              call draw_cylinder(x_body,xl3(j),-yl3(j),zl3(j),&
              xl3(j+1),-yl3(j+1),zl3(j+1),rl3(j),mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
            enddo
          endif
          if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="antennae_proboscis")) then
            call draw_cylinder(x_body,xan(1),-yan(1),zan(1),&
            xan(2),-yan(2),zan(2),ran,mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
          endif
        endif
        if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="antennae_proboscis")) then
          call draw_cylinder(x_body,xf(1),yf(1),zf(1),&
          xf(2),yf(2),zf(2),rf,mask(ix,iy,iz),mask_color(ix,iy,iz),color_body,Insect%safety)
        endif
      endif
    enddo
  enddo
enddo
endif

end subroutine draw_body_bumblebee

!------------------------------------------------------------------------------

! Body adapted from Maeda & Liu. It assumes Insect%x_head=0.0
subroutine draw_body_drosophila_maeda( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
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
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob)

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
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob)
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
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
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
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x_glob)

        ! ------------------------------------
        ! two b-splines body (abdomen+thorax)
        ! ------------------------------------
        x = x_body(1) + 0.8067d0 ! centers the thickest part of the thorax at the origin

        ! check if inside body bounds (in x-direction)
        if ( (x>=-Insect%safety) .and. (x<=1.2+Insect%safety) ) then
          R0=0.0d0
          ! compute radius as a function of x (counting from the tail on)
          if (x < 0.6333d0) then
            ! we're in the ABDOMEN
            R0 = max( -1.2990d0*x**2 + 0.9490d0*x + 0.0267d0, 0.0d0)
          elseif ((x >= 0.6333d0) .and. (x <=1.0d0 )) then
            ! we're in the THORAX
            R0 = max( -2.1667d0*x**2 + 3.4661d0*x - 1.2194d0, 0.0d0)
          elseif ((x >= 1.0d0) .and. (x <=1.2d0 )) then
            ! we're in the HEAD
            R0 = max( -12.68d0*x**2 + 27.4960d0*x - 14.7360d0, 0.0d0)
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
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
  integer :: ix,iy,iz

  ! the following are coordinates of specific points on the insect's body, for
  ! example the position of the head, its size etc. In older versions, these
  ! parameters were set in the *.ini file, which proved to be too much flexibility.
  ! in practice, the insect is created once, while implementing it, and then
  ! no longer changed. For Jerry, we overwrite the coordinates with hard-coded
  ! values here. the advantage is that we can set   BodyType=jerry;  and voilà!
  Insect%R_head = 0.125d0
  Insect%R_eye = 0.0625d0
  Insect%x_pivot_r =(/ 0.05d0, -0.2165d0, 0.d0 /)
  Insect%x_pivot_l =(/ 0.05d0, +0.2165d0, 0.d0 /)
  Insect%b_body = 0.1d0
  Insect%L_body = 1.0d0
  Insect%x_head = (/0.5d0*Insect%L_body,0.d0,0.d0 /)
  Insect%x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
  *0.8d0*(/1.d0,+1.d0,1.d0/)
  Insect%x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
  *0.8d0*(/1.d0,-1.d0,1.d0/)

  a_body = Insect%L_body / 2.d0
  !-----------------------------------------------------------------------------
  ! Jerry's body is an ellipsoid
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system, which is centered at Insect%xc_body
        x_body = matmul( M_body, x_glob)
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
  ! Jerry's head and eyes are spheres
  !-----------------------------------------------------------------------------
  x_head = Insect%xc_body + matmul(transpose(M_body),Insect%x_head)
  call drawsphere( x_head,Insect%R_head,mask,mask_color,us,Insect,color_body )

  x_eye = Insect%xc_body + matmul(transpose(M_body),Insect%x_eye_l)
  call drawsphere( x_eye,Insect%R_eye,mask,mask_color,us,Insect,color_body )

  x_eye = Insect%xc_body + matmul(transpose(M_body),Insect%x_eye_r)
  call drawsphere( x_eye,Insect%R_eye,mask,mask_color,us,Insect,color_body )
end subroutine draw_body_jerry


! a body that is just a sphere of unit diameter. used for particles.
subroutine draw_body_sphere( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: x,R0,R,R_tmp,x_tmp,a_body
  real(kind=pr) :: corner
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)

  x_head = Insect%xc_body
  call drawsphere( x_head,0.50d0,mask,mask_color,us,Insect,color_body )

end subroutine draw_body_sphere


! draw a cylinder defined by points (x1,y1,z1), (x2,y2,z2) and radius R0
subroutine draw_cylinder( xp,x1,y1,z1,x2,y2,z2,R0,mask_val,color_val,icolor,safety )
  use vars
  implicit none

  real(kind=pr),intent(in)::xp(1:3),x1,x2,y1,y2,z1,z2,R0,safety
  real(kind=pr),intent(inout)::mask_val
  integer(kind=2),intent(in)::icolor
  integer(kind=2),intent(inout)::color_val

  real(kind=pr)::x(1:3),R,xab,yab,zab,xu,yu,zu,xvp,yvp,zvp,&
  cbx,cby,cbz,rbx,rby,rbz

  ! draw cylinder without endpoint treatment
  cbx = 0.5*(x1+x2) - xp(1)
  cby = 0.5*(y1+y2) - xp(2)
  cbz = 0.5*(z1+z2) - xp(3)
  rbx = x1-x2
  rby = y1-y2
  rbz = z1-z2
  if ( cbx*cbx+cby*cby+cbz*cbz < 0.25*(rbx*rbx+rby*rby+rbz*rbz) ) then
    xab = xp(1)-x1
    yab = xp(2)-y1
    zab = xp(3)-z1
    xu = x2-x1
    yu = y2-y1
    zu = z2-z1
    xvp = yab*zu-zab*yu
    yvp = zab*xu-xab*zu
    zvp = xab*yu-yab*xu
    R = sqrt((xvp*xvp+yvp*yvp+zvp*zvp)/(xu*xu+yu*yu+zu*zu))
    if ( R <= R0+safety ) then
      mask_val = max(steps(R,R0),mask_val)
      color_val = icolor
    endif
  endif

  ! spheres at endpoints
  x = xp - (/x1,y1,z1/)
  if (dabs(x(1)) <= R0+safety) then
    if (dabs(x(2)) <= R0+safety) then
      if (dabs(x(3)) <= R0+safety) then
        R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
        if ( R <= R0+safety ) then
          mask_val = max(steps(R,R0),mask_val)
          color_val = icolor
        endif
      endif
    endif
  endif
  x = xp - (/x2,y2,z2/)
  if (dabs(x(1)) <= R0+safety) then
    if (dabs(x(2)) <= R0+safety) then
      if (dabs(x(3)) <= R0+safety) then
        R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
        if ( R <= R0+safety ) then
          mask_val = max(steps(R,R0),mask_val)
          color_val = icolor
        endif
      endif
    endif
  endif

end subroutine



!-------------------------------------------------------------------------------
subroutine draw_body_particle( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body, projected_length
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_part(1:3), n_part(1:3)
  integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k

  !-----------------------------------------------------------------------------
  ! initialization phase, executed only once
  !-----------------------------------------------------------------------------
  if (.not.allocated(particle_points)) then
    ! initialization, this is the first call to the routine.
    if (mpirank==0) then
      open(37, file='particle.in', form='formatted', status='old')
      read(37,*) npoints
      write(*,'("reading particle with ",i7," points")') npoints
    endif
    call MPI_BCAST( npoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
    allocate ( particle_points(1:npoints,1:6) )

    if (mpirank==0) then
      do ix = 1,npoints
        read(37,*) particle_points(ix,:)
      enddo
      write(*,*) "done reading particle.in (that is excellent news!)"
      close(37)
    endif
    ! make sure all cpu know the particle well
    call MPI_BCAST( particle_points(:,1),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
    call MPI_BCAST( particle_points(:,2),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
    call MPI_BCAST( particle_points(:,3),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
    call MPI_BCAST( particle_points(:,4),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
    call MPI_BCAST( particle_points(:,5),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
    call MPI_BCAST( particle_points(:,6),npoints,MPI_DOUBLE_PRECISION,0,&
    MPI_COMM_WORLD,mpicode )
  endif

  ! initialize signed distance as very far away
  mask = 1.d8

  !-----------------------------------------------------------------------------
  ! now, we are sure that all CPU know all points on the particle, and we can
  ! proceed to draw it
  !-----------------------------------------------------------------------------
  npoints = size(particle_points,1)
  ! loop over the marker points
  do ip = 1, npoints
    ! coordinate of surface point (in body system)
    x_part = particle_points(ip,1:3)
    ! normal vector of surface point (in body system)
    n_part = particle_points(ip,4:6)

    ! go to laboratory frame:
    x_glob = matmul( transpose(M_body), x_part) + Insect%xc_body
    ! periodize:
    if (x_glob(1)<0.0) x_glob(1)=x_glob(1)+xl
    if (x_glob(2)<0.0) x_glob(2)=x_glob(2)+yl
    if (x_glob(3)<0.0) x_glob(3)=x_glob(3)+zl
    if (x_glob(1)>=xl) x_glob(1)=x_glob(1)-xl
    if (x_glob(2)>=yl) x_glob(2)=x_glob(2)-yl
    if (x_glob(3)>=zl) x_glob(3)=x_glob(3)-zl

    ! coordinate (global) in integer space:
    ijk = nint( x_glob/dx )
    ! size of box around markers
    box = 1

    ! loop over neigborhood of marker
    do iz = ijk(3)-box, ijk(3)+box
      do iy = ijk(2)-box, ijk(2)+box
        do ix = ijk(1)-box, ijk(1)+box
          ! check if this point is on my rank. note this also checks implicitly
          ! if we're in the domain at all
          if ( on_proc( (/ix,iy,iz/)  ) ) then
            ! x_glob is in the global coordinate system
            x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
            x_glob = periodize_coordinate(x_glob - Insect%xc_body)
            ! x_body is in the body coordinate system, which is centered at Insect%xc_body
            x_body = matmul( M_body, x_glob )

            ! unsigned distance to point
            R = dsqrt(  (x_body(1)-x_part(1))*(x_body(1)-x_part(1)) + &
                        (x_body(2)-x_part(2))*(x_body(2)-x_part(2)) + &
                        (x_body(3)-x_part(3))*(x_body(3)-x_part(3)) )

            if ( R<=abs(mask(ix,iy,iz)) ) then
              ! this is closer, so we overwrite
              mask(ix,iy,iz) = R
              mask_color(ix,iy,iz) = color_body
              ! compute scalar product of difference vector and outward pointing normal
              projected_length = (x_body(1)-x_part(1))*n_part(1) + &
                                 (x_body(2)-x_part(2))*n_part(2) + &
                                 (x_body(3)-x_part(3))*n_part(3)

              if ( projected_length <= 0.d0  ) then
                ! we're inside the particle
                mask(ix,iy,iz) = -R!mask(ix,iy,iz)
              endif
            endif

          endif ! on_proc
        enddo ! nx
      enddo ! ny
    enddo ! nz
  enddo ! np


  ! do iz = ra(3), rb(3)
  !   do iy = ra(2), rb(2)
  !     do ix = ra(1), rb(1)
  !       if (mask(ix,iy,iz) < 0.d0) then
  !         R=0.d0
  !         do k=iz-1,iz+1
  !           do j=iy-1,iy+1
  !             do i=ix-1,ix+1
  !               R=R+mask( per(i,nx),per(j,ny),per(k,nz) )
  !             enddo
  !           enddo
  !         enddo
  !         R=R-mask(ix,iy,iz)
  !         if (R>26.5e8) mask(ix,iy,iz) = 1.0d8
  !
  !       endif
  !     enddo
  !   enddo
  ! enddo

  !-----------------------------------------------------------------------------
  ! fill the interior of the particle with "-" signs (the above concentrates
  ! on the interface!)
  ! we exploit the fact that the x-direction is continuous in memory and not
  ! split among processes
  !-----------------------------------------------------------------------------
  ! we start at a point which is surely not inside
  ! the particle. the algorithm cannot start on interior
  ! points.
  ! start = per(nint( (Insect%xc_body(1)-0.5*xl)/dx), nx)
  ! if(root) write(*,*) "point is", Insect%xc_body(1)-0.5*xl

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = start, start+nx ! we run all points, still
        ! is the point not yet touched (i.e. large value)?
        if (mask(per(ix,nx),iy,iz) > 99.99e6) then
          ! is either of the neighbors inside (e.g. negative)
          if (mask(per(ix+1,nx),iy,iz)<0.d0.or.mask(per(ix-1,nx),iy,iz)<0.d0) then
            mask(per(ix,nx),iy,iz) = -mask(per(ix,nx),iy,iz)
            mask_color(per(ix,nx),iy,iz) = color_body
          endif
        endif
      enddo
    enddo
  enddo

  ! call save_field_hdf5(0.d0,'./mask_00',mask)
  ! stop

  !-----------------------------------------------------------------------------
  ! convert signed distance function to mask function chi
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        mask(ix,iy,iz) = steps( mask(ix,iy,iz),0.d0 )
      enddo
    enddo
  enddo


end subroutine draw_body_particle


!-------------------------------------------------------------------------------
! In our terminology, a macroscopic particle is an insect without wings and no
! flapping motion in free flight. Therefore, the insect module contains nowadays
! also body shapes that are not related to insects. This one is a flat plate of
! size
! Insect%L_span x Insect%L_chord x Insect%WingThickness
!-------------------------------------------------------------------------------
subroutine draw_body_platicle( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body, projected_length
  real(kind=pr) :: x_body(1:3), x(1:3), xc(1:3), n_part(1:3)
  integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x is in the global coordinate system
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x is now centered in the plate's center point
        x = periodize_coordinate(x - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x)

        ! bounding box checks
        if (dabs(x_body(1)) <= Insect%L_span+Insect%safety) then
          if (dabs(x_body(2)) <= Insect%L_chord+Insect%safety) then
            if (dabs(x_body(3)) <= Insect%WingThickness+Insect%safety) then
              ! signed distance:
              R = maxval( (/ dabs(x_body(3))-Insect%WingThickness/2.d0,&
                             dabs(x_body(2))-Insect%L_chord/2.d0,&
                             dabs(x_body(1))-Insect%L_span/2.d0 &
                          /) )
              mask(ix,iy,iz) = max(steps(R,0.d0),mask(ix,iy,iz))
              mask_color(ix,iy,iz) = color_body
            endif
          endif
        endif


      enddo
    enddo
  enddo

end subroutine draw_body_platicle


!-------------------------------------------------------------------------------
! In our terminology, a macroscopic particle is an insect without wings and no
! flapping motion in free flight. Therefore, the insect module contains nowadays
! also body shapes that are not related to insects. This one is a flat COIN (D=1)
!-------------------------------------------------------------------------------
subroutine draw_body_coin( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body, projected_length
  real(kind=pr) :: x_body(1:3), x(1:3), xc(1:3), n_part(1:3)
  integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x is in the global coordinate system
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x is now centered in the sphere's center point
        x = periodize_coordinate(x - Insect%xc_body)
        ! x_body is in the body coordinate system
        x_body = matmul(M_body,x)

        if (dabs(x_body(1)) <= 0.5d0+Insect%safety) then
          if (dabs(x_body(2)) <= 0.5d0+Insect%safety) then
            if (dabs(x_body(3)) <= Insect%WingThickness+Insect%safety) then
              ! signed distance:
              R = maxval( (/ dabs(x_body(3)) - Insect%WingThickness/2.d0,&
                             dsqrt(x_body(2)**2 + x_body(1)**2)-0.5d0 &
                          /) )
              mask(ix,iy,iz) = max(steps(R,0.d0),mask(ix,iy,iz))
              mask_color(ix,iy,iz) = color_body
            endif
          endif
        endif


      enddo
    enddo
  enddo

end subroutine draw_body_coin



!-------------------------------------------------------------------------------
! Thin rod-like body used in Suzuki et al. JFM 2015 to model a butterfly
!-------------------------------------------------------------------------------
subroutine draw_suzuki_thin_rod( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  real(kind=pr),intent(inout)::us(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ga(1):gb(1),ga(2):gb(2),ga(3):gb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a,RR0
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
  integer :: ix,iy,iz

  R0 = ( 0.5d0*Insect%WingThickness + Insect%Safety )**2
  RR0 = 0.5d0*Insect%WingThickness

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system, which is centered at Insect%xc_body
        x_body = matmul( M_body, x_glob)

        if ( dabs(x_body(1))<=0.5d0+Insect%safety) then
          R = x_body(2)**2 + x_body(3)**2
          if ( R < R0) then
            a = steps(dsqrt(R),RR0)
            if (mask(ix,iy,iz)<=a) then
              mask(ix,iy,iz) = a
              mask_color(ix,iy,iz) = color_body
            endif
           endif
         endif

      enddo
    enddo
  enddo

end subroutine draw_suzuki_thin_rod

!-------------------------------------------------------------------------------
subroutine draw_body_hawkmoth( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3), x_eye_r(1:3), x_eye_l(1:3)
  real(kind=pr), dimension(1:3) :: x1,x2
  integer :: ix,iy,iz

  Insect%R_head = 0.125d0
  Insect%R_eye = 0.0625d0
  ! Insect%x_pivot_r =(/ 0.05d0, -0.2165d0, 0.d0 /)
  ! Insect%x_pivot_l =(/ 0.05d0, +0.2165d0, 0.d0 /)
  Insect%b_body = 0.15d0
  Insect%L_body = 1.0d0
  Insect%x_head = (/0.5d0*Insect%L_body,0.d0,0.d0 /)
  Insect%x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
  *0.8d0*(/1.d0,+1.d0,1.d0/)
  Insect%x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
  *0.8d0*(/1.d0,-1.d0,1.d0/)

  a_body = Insect%L_body / 2.d0
  !-----------------------------------------------------------------------------
  ! The body is an ellipsoid
  !-----------------------------------------------------------------------------
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system, which is centered at Insect%xc_body
        x_body = matmul( M_body, x_glob)
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
  ! Head is a sphere, we add antennae, which are cylinders
  !-----------------------------------------------------------------------------
  x_head = Insect%xc_body + matmul(transpose(M_body),Insect%x_head)
  call drawsphere( x_head,Insect%R_head,mask,mask_color,us,Insect,color_body )

  ! these guys are in the body system:
  x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.d0,+1.d0,1.d0/)
  x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*1.8d0*(/1.d0,+1.d0,1.d0/)
  ! back to global system
  x1 = Insect%xc_body + matmul(transpose(M_body),x_eye_l)
  x2 = Insect%xc_body + matmul(transpose(M_body),x_eye_r)
  ! draw the cylinder (with spheres at the ends)
  call draw_cylinder_new( x1, x2, 0.015d0*1.3d0, mask, mask_color, us, Insect, color_body )


  ! these guys are in the body system:
  x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.d0,-1.d0,1.d0/)
  x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*1.8d0*(/1.d0,-1.d0,1.d0/)
  ! back to global system
  x1 = Insect%xc_body + matmul(transpose(M_body),x_eye_l)
  x2 = Insect%xc_body + matmul(transpose(M_body),x_eye_r)
  ! draw the cylinder (with spheres at the ends)
  call draw_cylinder_new( x1, x2, 0.015d0*1.3d0, mask, mask_color, us, Insect, color_body )
end subroutine draw_body_hawkmoth



!-------------------------------------------------------------------------------
! The mosquito is based on the simplified model presented in
! [1] Iams "Flight stability of mosquitos: A reduced model" SIAM J. Appl. Math. 74(5) 1535--1550 (2014)
!-------------------------------------------------------------------------------
subroutine draw_body_mosquito_iams( mask, mask_color, us, Insect, color_body, M_body)
  use vars
  implicit none

  type(diptera),intent(inout) :: Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_body
  real(kind=pr),intent(in)::M_body(1:3,1:3)

  real(kind=pr) :: R0,R,a_body, a,b,c, alpha, Ralpha(1:3,1:3)
  real(kind=pr) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3), x_eye_r(1:3), x_eye_l(1:3)
  real(kind=pr) :: x0_head(1:3), x0_abdomen(1:3), x0_thorax(1:3)
  real(kind=pr), dimension(1:3) :: x1,x2
  integer :: ix,iy,iz

  ! The mosquito consists of three parts: head, thorax and abdomen (sphere, ellipsoid, ellipsoid)
  ! positions are measured from fig. 1 in [1], we computed also the center of gravity
  ! for this mosquito, Insect%xc_body is thus the center of gravity
  x0_head = (/ 0.5652d0, 0.d0, -0.0434d0 /)
  x0_thorax = (/ 0.2579d0, 0.d0, 0.1267d0 /)
  x0_abdomen = (/-0.437d0, 0.d0, -0.2024d0 /)

  !-----------------------------------------------------------------------------
  ! head
  !-----------------------------------------------------------------------------
  ! the head is a simple sphere with radius 0.1154
  R0 = 0.1154d0
  x1 = x0_head + Insect%xc_body
  call drawsphere( x1, R0, mask,mask_color,us,Insect,color_body )

  !-----------------------------------------------------------------------------
  ! thorax
  !-----------------------------------------------------------------------------
  ! the thorax is a triaxial ellipsiod without rotation
  a = 0.2628d0
  b = 0.1603d0
  c = b ! HACK: for simplicity, assume b=c, otherwise it can be very tough to draw

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system, which is centered at Insect%xc_body
        x_body = matmul( M_body, x_glob)
        ! translate to origin of thorax
        x_body = x_body - x0_thorax

        ! check if inside the surrounding box (save comput. time)
        if ( dabs(x_body(1)) <= b + Insect%safety ) then
        if ( dabs(x_body(2)) <= b + Insect%safety ) then
        if ( dabs(x_body(3)) <= a + Insect%safety ) then
          ! the x-y plane are circles
          R  = dsqrt ( x_body(1)**2 + x_body(2)**2 )
          ! this gives the R(x) shape
          if ( x_body(3)/a <= 1.d0) then
            R0 = b * dsqrt( 1.d0 - (x_body(3)/a)**2 )
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
  ! abdomen
  !-----------------------------------------------------------------------------
  ! the abdomen is a axi-symmetric ellipsiod inclined by 30.44°
  a = 0.6026d0
  b = 0.1282d0
  ! angle by which the abdomen is tilted (measured from figure 1 in [1])
  alpha = deg2rad(-30.44d0)

  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        x_glob = periodize_coordinate(x_glob - Insect%xc_body)
        ! x_body is in the body coordinate system, which is centered at Insect%xc_body
        x_body = matmul(M_body, x_glob)
        ! translate to origin of abdomen
        x_body = x_body - x0_abdomen
        ! rotate into abdomens principal axis
        call Ry(Ralpha,alpha)
        x_body = matmul(Ralpha, x_body)

        ! check if inside the surrounding box (save comput. time)
        if ( dabs(x_body(1)) <= a + Insect%safety ) then
        if ( dabs(x_body(2)) <= b + Insect%safety ) then
        if ( dabs(x_body(3)) <= b + Insect%safety ) then
          ! the y-z plane are circles
          R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
          ! this gives the R(x) shape
          if ( x_body(1)/a <= 1.d0) then
            R0 = b * dsqrt( 1.d0 - (x_body(1)/a)**2 )
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



end subroutine draw_body_mosquito_iams





!-------------------------------------------------------------------------------
! draw a cylinder defined by GLOBALS points (x1,y1,z1), (x2,y2,z2) and radius R0
! At the start/end point, we add a sphere.
! The solid velocity field us is not touched -- we consider this routine for bodies
! therefore the solid velocity field (which is a solid body rotation around
! insect%xc) is added in the main insect drawing routine.
! The color of the new cylinder will be what you pass in color_val
!-------------------------------------------------------------------------------
subroutine draw_cylinder_new( x1, x2, R0, mask, mask_color, us, Insect, color_val)
  use vars
  implicit none

  real(kind=pr),dimension(1:3),intent(inout )::x1,x2
  real(kind=pr),intent(in)::R0
  type(diptera),intent(inout)::Insect
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: color_val

  real(kind=pr),dimension(1:3)::x_glob, e_x, tmp
  real(kind=pr)::ceta, R, clength, safety, t
  integer :: ix,iy,iz

  safety = Insect%safety

  ! unit vector in cylinder axis direction and cylinder length
  e_x = (x2-x1)
  clength = norm2(e_x)
  e_x = e_x / clength

  ! first we draw the cylinder, then the endpoint spheres
  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x_glob is in the global coordinate system
        x_glob = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)

        ! position on cylinder axis (projection of x-x1 on e_x)
        ceta = dot_product(x_glob-x1, e_x)

        ! is the height inside the bounding box?
        if ( ceta >= -safety .and. ceta <= clength+safety) then
          ! we're maybe inside the cylinder
          ! Radius is the distance to centreline
          tmp = x_glob - (x1 + ceta*e_x)
          R = (tmp(1)**2 + tmp(2)**2 + tmp(3)**2)

          if ( R <= (R0+safety)**2 ) then
            R = dsqrt(R)
            ! the new value we'd love to set here
            t = steps(R,R0)
            if (t>= mask(ix,iy,iz)) then
              mask(ix,iy,iz) = t
              mask_color(ix,iy,iz) = color_val
            endif
          endif
        endif
      enddo
    enddo
  enddo

  ! spheres at endpoints
  call drawsphere( x1, R0, mask,mask_color,us,Insect,color_val )
  call drawsphere( x2, R0, mask,mask_color,us,Insect,color_val )
end subroutine draw_cylinder_new



!-------------------------------------------------------------------------------
! draw a sphere with radius R0 centered at the point xc (GLOBAL SYSTEM)
! This routiner's intended use is for drawing the insect's body, for example
! the head and eyes of Jerry. The velocity field inside the body is added
! later, thus, the field us is untouched in this routines.
!-------------------------------------------------------------------------------
subroutine drawsphere( xc,R0,mask,mask_color,us,Insect,icolor )
  use vars
  implicit none

  real(kind=pr),intent(inout)::xc(1:3)
  real(kind=pr),intent(in)::R0
  real(kind=pr),intent(inout)::mask(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  real(kind=pr),intent(inout)::us(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3),1:neq)
  integer(kind=2),intent(inout)::mask_color(ra(1):rb(1),ra(2):rb(2),ra(3):rb(3))
  integer(kind=2),intent(in) :: icolor
  type(diptera),intent(inout) :: Insect

  integer :: ix,iy,iz
  real(kind=pr)::x(1:3),R,tmp

  ! periodization: if the center point is out of the domain, then correct that
  if (xc(1)<0.0) xc(1)=xc(1)+xl
  if (xc(2)<0.0) xc(2)=xc(2)+yl
  if (xc(3)<0.0) xc(3)=xc(3)+zl

  if (xc(1)>=xl) xc(1)=xc(1)-xl
  if (xc(2)>=yl) xc(2)=xc(2)-yl
  if (xc(3)>=zl) xc(3)=xc(3)-zl



  do iz = ra(3), rb(3)
    do iy = ra(2), rb(2)
      do ix = ra(1), rb(1)
        ! x is in the global coordinate system
        x = (/ dble(ix)*dx, dble(iy)*dy, dble(iz)*dz /)
        ! x is now centered in the sphere's center point
        x = periodize_coordinate(x - xc)

        ! bounding box check
        if (dabs(x(1)) <= R0+Insect%safety) then
          if (dabs(x(2)) <= R0+Insect%safety) then
            if (dabs(x(3)) <= R0+Insect%safety) then
              ! compute radius
              R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
              if ( R <= R0+Insect%safety ) then
                tmp = steps(R,R0)
                if (tmp>=mask(ix,iy,iz)) then
                  ! set new value
                  mask(ix,iy,iz) = tmp
                  mask_color(ix,iy,iz) = icolor
                endif
              endif
            endif
          endif
        endif


      enddo
    enddo
  enddo

end subroutine
