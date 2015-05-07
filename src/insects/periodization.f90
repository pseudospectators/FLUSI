!-------------------------------------------------------------------------------
subroutine init_periodic_insect(Insect)
  use vars
  type(diptera),intent(inout) :: Insect
  integer :: i,j,k,q
  real(kind=pr) :: perio(1:3), ghosts(1:3,1:27), R0

  ! if periodization is turned off, this routine does nothing
  if (.not.Insect%periodic) return

  ! First step: collect ALL virtual centers in the array ghosts. there are 27 of
  ! them, but we'll later figure out which ones we need. The ghosts are located
  ! at xc (+-xl, +-yl, +-zl), and we do that in loops for simplicity
  perio=(/-1.0,0.0,1.0/)

  q=1
  do i=1,3
    do j=1,3
      do k=1,3
        ghosts(1,q) = Insect%xc_body(1) + xl*perio(i)
        ghosts(2,q) = Insect%xc_body(2) + yl*perio(j)
        ghosts(3,q) = Insect%xc_body(3) + zl*perio(k)
        q=q+1
      enddo
    enddo
  enddo

  ! second step. Check if the box, centered around the ghost center, of size R0
  ! touches the compuational domain, ie if either of the eight corners of that
  ! box lies within the usual domain. If so, the ghost image is important, if
  ! not, we can drop it. The first task is to figure how how many ghost images
  ! are relevant, then we can allocate the memory to store them
  R0 = 1.2d0
  Insect%n_relevant_ghosts = 0
  ! note we also loop over the ORIGINAL image, which is the only point that actually
  ! lies inside the domain. The number of revelant ghost is thus always >=1 and never 0
  do q = 1,27
    if (  in_domain( ghosts(:,q)+(/-R0,-R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,+R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,-R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,+R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,-R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,+R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,-R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,+R0,-R0/) ) ) then
      ! the insect datatype holds the number of ghosts relevant
      Insect%n_relevant_ghosts = Insect%n_relevant_ghosts + 1
    endif
  enddo
  ! we now know how many of the 27 ghosts are in fact relevant to us

  ! third step, allocate the buffer and store the relevant ghost images
  if (allocated(relevant_ghosts)) deallocate (relevant_ghosts)
  allocate(relevant_ghosts(1:3,1:Insect%n_relevant_ghosts))
  ! the array dist will be used to figure out which is the closest ghost image
  ! when looping over the cartesian grid
  if (allocated(dists)) deallocate (dists)
  allocate(dists(1:Insect%n_relevant_ghosts))

  ! save relevant ghost images in the array "relevant_ghosts" which is a member
  ! of the insect module
  i=1
  do q = 1,27
    if (  in_domain( ghosts(:,q)+(/-R0,-R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,+R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,-R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/-R0,+R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,-R0,-R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,+R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,-R0,+R0/) ) .or. &
          in_domain( ghosts(:,q)+(/+R0,+R0,-R0/) ) ) then

      relevant_ghosts(1:3,i) = ghosts(:,q)
      i=i+1
    endif
  enddo

end subroutine

!-------------------------------------------------------------------------------
! free memory occupied by relevant ghosts, since in the next time step, we will
! count again and maybe that number has changed.
!-------------------------------------------------------------------------------
subroutine free_periodic_insect()
  ! deallocate(relevant_ghosts,dists)
end subroutine



!-------------------------------------------------------------------------------
! Get nearest center point.
!-------------------------------------------------------------------------------
! Due to periodicity, we do not have only one insect, but rather 27 identical
! ones surrounding it. In fact, it is even an infinite number, but let's not
! be too picky. The 27 neighbors lie in 27 boxes at each side and corner of the
! domain. Now, to draw a periodic insect, which leaves the domain at one side
! and re-enters on the other, we simply consider the closest of these 27
! insects at each Cartesian grid point. This routine returns the nearest center.
!-------------------------------------------------------------------------------
function get_nearest_center(Insect, x_glob)
  use vars
  type(diptera),intent(inout) :: Insect
  real(kind=pr),dimension(1:3) :: get_nearest_center
  real(kind=pr),dimension(1:3) :: x_glob
  integer :: q, qq(1)


  if (Insect%periodic) then
    if (Insect%n_relevant_ghosts == 1) then
      ! this is the usual case, the insect is far away from the boundary
      get_nearest_center = Insect%xc_body
      return
    endif
    ! look for nearest image
    do q = 1, Insect%n_relevant_ghosts
      dists(q) = dsqrt( (x_glob(1)-relevant_ghosts(1,q))*(x_glob(1)-relevant_ghosts(1,q)) &
                      + (x_glob(2)-relevant_ghosts(2,q))*(x_glob(2)-relevant_ghosts(2,q)) &
                      + (x_glob(3)-relevant_ghosts(3,q))*(x_glob(3)-relevant_ghosts(3,q)) )
    enddo
    qq = minloc( dists )
    get_nearest_center(1) = relevant_ghosts(1,qq(1))
    get_nearest_center(2) = relevant_ghosts(2,qq(1))
    get_nearest_center(3) = relevant_ghosts(3,qq(1))
  else
    ! just return body coordinate, if not periodized
    get_nearest_center = Insect%xc_body
  endif


end function
